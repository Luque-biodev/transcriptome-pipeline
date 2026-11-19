#!/usr/bin/env Rscript
# Análisis diferencial (DESeq2) + enriquecimiento GO/KEGG por contraste.
# Uso: Rscript rnapipe.R <project> <go_mode> <padj_cutoff> <lfc_cutoff> [organism] [gene_id_type]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
  stop("Uso: Rscript rnapipe.R <project> <go_mode> <padj_cutoff> <lfc_cutoff> [organism] [gene_id_type]")

PROJECT    <- args[1]
GO_MODE    <- args[2]                       # BP / MF / CC / ALL
PVALUE_CUT <- as.numeric(args[3])           # umbral sobre padj
LFC_CUT    <- as.numeric(args[4])           # umbral |log2FC|
ORGANISM   <- tolower(if (length(args) >= 5) args[5] else "auto")
GENEID     <- toupper(if (length(args) >= 6) args[6] else "auto")

if (!dir.exists(PROJECT)) stop("No existe el directorio del proyecto: ", PROJECT)
if (!GO_MODE %in% c("BP", "MF", "CC", "ALL")) stop("go_mode debe ser BP, MF, CC o ALL (recibido: ", GO_MODE, ")")
if (is.na(PVALUE_CUT) || is.na(LFC_CUT)) stop("Los umbrales padj y lfc deben ser numéricos")
setwd(PROJECT)
suppressPackageStartupMessages(library(DESeq2))

############################
## ORGANISMO              ##
############################
ORGS <- list(
  ath = list(pkg = "org.At.tair.db", kegg = "ath", key = "TAIR"),
  hsa = list(pkg = "org.Hs.eg.db",   kegg = "hsa", key = "ENSEMBL"),
  mmu = list(pkg = "org.Mm.eg.db",   kegg = "mmu", key = "ENSEMBL")
)

if (ORGANISM == "auto") {   # las cabeceras de un GTF rara vez lo indican: mejor fijar `organism` en params
  hdr <- tolower(readLines("genome/annot.gtf", n = 200))
  ORGANISM <- if (any(grepl("arabidopsis|tair|araport", hdr))) {
    "ath"
  } else if (any(grepl("homo sapiens|hg38|grch38", hdr))) {
    "hsa"
  } else if (any(grepl("mus musculus|mm10|grcm", hdr))) {
    "mmu"
  } else "none"
}

cfg <- ORGS[[ORGANISM]]
if (is.null(cfg)) {
  message("Organismo no soportado/detectado: se omite GO y KEGG.")
} else {
  suppressPackageStartupMessages({
    library(clusterProfiler); library(enrichplot); library(cfg$pkg, character.only = TRUE)
  })
  OrgDb <- get(cfg$pkg)
  KEY   <- if (GENEID == "AUTO") cfg$key else GENEID
}

# Quita sufijos de prepDE ("id|nombre") y versiones Ensembl (".3")
clean_ids <- function(x) {
  x <- sub("\\|.*$", "", x)
  if (KEY == "ENSEMBL") x <- sub("\\.[0-9]+$", "", x)
  unique(x)
}

kegg_ids <- function(ids) {
  if (cfg$kegg == "ath" || KEY == "ENTREZID") return(ids)   # KEGG de Arabidopsis usa locus TAIR
  tryCatch(unique(bitr(ids, KEY, "ENTREZID", OrgDb)$ENTREZID), error = \(e) character())
}

# Abre el dispositivo, imprime el gráfico y lo cierra aunque falle
save_png <- function(file, plot) {
  png(file, width = 1500, height = 1500, res = 200)
  on.exit(invisible(dev.off()))
  print(plot)
}

############################
## COUNTS Y DISEÑO        ##
############################
counts <- read.csv("results/gene_count_matrix.csv", row.names = 1, check.names = FALSE)
design <- read.table("results/design.txt", sep = "\t", col.names = c("sample", "condition"),
                     stringsAsFactors = FALSE)
design$condition <- factor(design$condition, levels = unique(design$condition))  # 1ª = control
rownames(design) <- design$sample
if (!all(design$sample %in% colnames(counts)))
  stop("Muestras del diseño ausentes en gene_count_matrix.csv: ",
       paste(setdiff(design$sample, colnames(counts)), collapse = ", "))
counts <- counts[, design$sample]

############################
## DESEQ2                 ##
############################
dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

write.csv(counts(dds, normalized = TRUE), "results/normalized_counts.csv")

vsd <- tryCatch(vst(dds, blind = TRUE), error = \(e) varianceStabilizingTransformation(dds, blind = TRUE))
save_png("results/PCA.png", plotPCA(vsd, intgroup = "condition"))

control    <- levels(design$condition)[1]
treatments <- levels(design$condition)[-1]

############################
## FUNCIONES AUXILIARES   ##
############################
run_go <- function(genes, universe) {
  if (length(genes) == 0) return(NULL)
  tryCatch(enrichGO(gene = genes, universe = universe, OrgDb = OrgDb, keyType = KEY, ont = GO_MODE,
                    pvalueCutoff = PVALUE_CUT, qvalueCutoff = PVALUE_CUT, readable = KEY != "SYMBOL"),
           error = \(e) { message("enrichGO falló: ", conditionMessage(e)); NULL })
}

run_kegg <- function(genes) {
  genes <- kegg_ids(genes)
  if (length(genes) == 0) return(NULL)
  tryCatch(enrichKEGG(genes, organism = cfg$kegg, pvalueCutoff = PVALUE_CUT, qvalueCutoff = PVALUE_CUT),
           error = \(e) { message("enrichKEGG falló (¿sin internet?): ", conditionMessage(e)); NULL })
}

save_enrich <- function(obj, kind, dir, ct) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) {
    writeLines(paste("No enriched", kind, "terms detected"),
               sprintf("results/%s_%s_%s_EMPTY.txt", kind, dir, ct))
    return(invisible())
  }
  write.csv(as.data.frame(obj), sprintf("results/%s_%s_%s.csv", kind, dir, ct), row.names = FALSE)
  for (p in c("barplot", "dotplot")) {
    save_png(sprintf("results/%s_%s_%s_%s.png", kind, dir, p, ct),
             if (p == "barplot") barplot(obj, showCategory = 20) else dotplot(obj, showCategory = 20))
  }
}

############################
## ANÁLISIS POR CONTRASTE ##
############################
for (treat in treatments) {
  ct  <- paste0(treat, "_vs_", control)
  res <- results(dds, contrast = c("condition", treat, control))
  write.csv(as.data.frame(res), paste0("results/DE_", ct, ".csv"))

  ok  <- !is.na(res$padj)
  up  <- ok & res$padj < PVALUE_CUT & res$log2FoldChange >  LFC_CUT
  dn  <- ok & res$padj < PVALUE_CUT & res$log2FoldChange < -LFC_CUT
  activated <- rownames(res)[up]
  repressed <- rownames(res)[dn]
  write.csv(data.frame(gene = activated), paste0("results/activated_", ct, ".csv"), row.names = FALSE)
  write.csv(data.frame(gene = repressed), paste0("results/repressed_", ct, ".csv"), row.names = FALSE)
  message(ct, ": ", length(activated), " activados, ", length(repressed), " reprimidos")

  # Volcano
  x <- res$log2FoldChange; y <- -log10(res$padj)
  png(paste0("results/volcano_", ct, ".png"), width = 1500, height = 1500, res = 200)
  plot(x[ok], y[ok], pch = 19, col = "grey", xlab = "log2 Fold Change", ylab = "-log10 adjusted p-value")
  points(x[up], y[up], col = "red",  pch = 19)
  points(x[dn], y[dn], col = "blue", pch = 19)
  abline(v = c(-LFC_CUT, LFC_CUT), h = -log10(PVALUE_CUT), lty = 2)
  invisible(dev.off())

  if (is.null(cfg)) next

  # GO (universo = genes testados)
  universe <- clean_ids(rownames(res)[ok])
  save_enrich(run_go(clean_ids(activated), universe), "GO", "up",   ct)
  save_enrich(run_go(clean_ids(repressed), universe), "GO", "down", ct)

  # KEGG
  save_enrich(run_kegg(clean_ids(activated)), "KEGG", "up",   ct)
  save_enrich(run_kegg(clean_ids(repressed)), "KEGG", "down", ct)
}

cat("Análisis en R completado correctamente.\n")
