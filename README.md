# rnapipe

Pipeline de RNA-seq comparativo para clúster SLURM. Parte de FASTQ, un genoma y una anotación, y produce expresión diferencial y enriquecimiento funcional:

**FastQC → STAR → StringTie → gffcompare → prepDE → DESeq2 → GO/KEGG**

Este documento explica cómo se usa, qué hace cada script y dónde tocar para adaptarlo.

---

## 1. Flujo

```
params.txt
   │
   ▼
rnapipe.sh  (nodo de login)
   ├─ valida parámetros, muestras (SE/PE) y diseño
   ├─ crea el espacio de trabajo y copia genoma, anotación y reads
   ├─ GFF → GTF (gffread) si hace falta
   ├─ STAR genomeGenerate (índice)
   ├─ sbatch × N  ──►  sample_proc.sh   (FastQC → STAR → samtools index → StringTie -e)
   └─ sbatch --dependency=afterok:<todos los jobs>
                       │
                       ▼
              transcriptome_merge.sh   (StringTie --merge → gffcompare → prepDE.py)
                       │
                       ▼
                  rnapipe.R            (DESeq2 → volcano/PCA → GO/KEGG por contraste)
```

Si cualquier muestra falla, el merge **no arranca** (`afterok` + `--kill-on-invalid-dep=yes`).

---

## 2. Archivos

| Archivo | Función |
|---|---|
| `rnapipe.sh` | Punto de entrada. Lee `params.txt`, prepara todo y envía los jobs. |
| `sample_proc.sh` | Job por muestra: FastQC, STAR, índice del BAM, StringTie `-e`. |
| `transcriptome_merge.sh` | Job único: merge de GTFs, gffcompare, `prepDE.py`, llama a R. |
| `rnapipe.R` | Análisis diferencial y enriquecimiento. |
| `params.txt` | Ejemplo de fichero de parámetros. |

Los cuatro scripts deben estar en el mismo directorio (`installation_dir`).

---

## 3. Requisitos

En el PATH de los nodos de cómputo (y del login para `gffread` y `STAR`):

- `STAR`, `samtools`, `fastqc`, `stringtie`, `gffcompare`, `gffread`
- `prepDE.py` o `prepDE.py3` (incluido con StringTie)
- `Rscript` con los paquetes: `DESeq2`, `clusterProfiler`, `enrichplot` y el `org.*.db` de tu especie (`org.At.tair.db`, `org.Hs.eg.db`, `org.Mm.eg.db`)
- SLURM (`sbatch`)
- KEGG necesita **conexión a internet** en el nodo donde corre R. Sin ella se omite KEGG con un aviso y el resto sigue.

---

## 4. Uso

```bash
cp params.txt mi_experimento.params     # editar con tu nombre de experimento
bash /ruta/a/rnapipe/rnapipe.sh mi_experimento.params
```

Sigue la ejecución con `squeue -u $USER` y los logs en `<working_dir>/<experiment_name>/logs/`.

`rnapipe.sh` se niega a correr si `<working_dir>/<experiment_name>` ya existe, para no sobrescribir resultados.

---

## 5. Parámetros (`params.txt`)

Formato: `clave valor`, una por línea. Las líneas que empiezan por `#` y los comentarios tras ` #` se ignoran. Se toleran finales de línea de Windows.

### Obligatorios

| Clave | Descripción |
|---|---|
| `working_dir` | Carpeta donde se creará el experimento. |
| `experiment_name` | Nombre de la carpeta del experimento. |
| `installation_dir` | Carpeta con los scripts de rnapipe. |
| `genome_path` | FASTA del genoma. |
| `annot_path` | Anotación `.gtf`, `.gff` o `.gff3` (GFF se convierte a GTF). |
| `experimental_design` | Una condición por muestra, en el orden de las muestras. **La primera condición es el control.** |

### Muestras (una de las dos opciones)

| Clave | Descripción |
|---|---|
| `path_sample_N` | Una línea por muestra (`path_sample_1`, `path_sample_2`, …). Basta indicar el R1. |
| `samples_path` | Carpeta que contiene solo las muestras. Se usa el orden alfabético. |

Si hay `path_sample_*`, tiene prioridad sobre `samples_path`.

### Opcionales

| Clave | Por defecto | Descripción |
|---|---|---|
| `author_name` | `$USER` | Se añade al nombre de los jobs. |
| `organism` | `auto` | `ath`, `hsa`, `mmu`, `none` o `auto`. Con `none` (o si `auto` no detecta nada) se omiten GO y KEGG. Recomendado fijarlo a mano: las cabeceras de un GTF rara vez indican la especie. |
| `gene_id_type` | según organismo | Tipo de ID de gen para GO: `TAIR` (ath), `ENSEMBL` (hsa/mmu), `SYMBOL`, `ENTREZID`. |
| `go_mode` | `BP` | `BP`, `MF`, `CC` o `ALL`. |
| `pvalue_cutoff` | `0.05` | Umbral sobre `padj`. También se usa en GO/KEGG. |
| `lfc_cutoff` | `1` | Umbral de \|log2FC\|. |
| `threads` | `8` | Hilos de STAR/StringTie/FastQC y `--cpus-per-task`. |
| `sa_index_nbases` | `11` | `--genomeSAindexNbases` de STAR. Ajústalo al tamaño del genoma (≈ min(14, log2(tamaño)/2 − 1)). |
| `read_length` | `75` | Longitud de lectura para `prepDE.py -l`. |

---

## 6. Reglas de entrada

- **Reads comprimidos**: deben ser `.gz`. Se copian como `sample_N_R1.fq.gz` (y `_R2`).
- **Detección paired-end**: a partir del R1 se busca el mate sustituyendo `_R1`→`_R2` o `_1`→`_2` (antes de `.fq`/`.fastq`). Si el mate existe → PE; si no → SE.
- En modo `samples_path` y `path_sample_*` se ignoran los ficheros que parecen R2, para no contarlos como muestra.
- Se permite mezclar SE y PE en un experimento, pero ten en cuenta que no es ideal comparar entre ambos tipos.
- **Diseño**: nº de condiciones = nº de muestras; mínimo 2 condiciones y 2 réplicas por condición. Un número desigual de réplicas solo da aviso.

---

## 7. Estructura de salida

```
<working_dir>/<experiment_name>/
├── genome/
│   ├── genome.fa, annot.gtf (y annot.gff si venía en GFF)
│   └── index/                          # índice de STAR
├── samples/sample_N/
│   ├── sample_N_R1.fq.gz [_R2.fq.gz]
│   ├── *_fastqc.html / .zip
│   ├── sample_N_Aligned.sortedByCoord.out.bam (+ .bai)
│   ├── sample_N.gtf, sample_N.tsv      # StringTie
│   └── logs de STAR
├── results/
│   ├── samples.tsv                     # id, SE/PE, condición, R1, R2 originales
│   ├── design.txt, mergelist.txt, prepde_list.txt
│   ├── stringtie_merged.gtf, merged.*  # merge y gffcompare
│   ├── gene_count_matrix.csv, transcript_count_matrix.csv
│   ├── normalized_counts.csv, PCA.png
│   ├── DE_<trat>_vs_<ctrl>.csv
│   ├── activated_*.csv, repressed_*.csv, volcano_*.png
│   ├── GO_up_* / GO_down_*   (.csv, barplot, dotplot)
│   └── KEGG_up_* / KEGG_down_* (.csv, barplot, dotplot)
└── logs/
    ├── out_sample_N, err_sample_N
    └── out_tm, err_tm                  # merge + R
```

Cuando un enriquecimiento no da términos se escribe un `*_EMPTY.txt` en lugar del CSV y los gráficos.

---

## 8. Análisis estadístico

- Diseño `~ condition`; la primera condición del `experimental_design` es el control y el resto se compara contra ella (un contraste por condición).
- Se filtran genes con menos de 10 lecturas  totales. Aunque tenga condiciones con 0 lecturas.
- Activados: `padj < pvalue_cutoff` y `log2FC > lfc_cutoff`. Reprimidos: `padj < pvalue_cutoff` y `log2FC < -lfc_cutoff`.
- GO: universo = genes testados. IDs limpios de sufijos de `prepDE` (`id|nombre`) y de versiones Ensembl.
- KEGG: para `hsa`/`mmu` los IDs se convierten a ENTREZ con `bitr`; para `ath` se usan los locus TAIR.
- La cuantificación es contra la anotación de referencia (`stringtie -e -G annot.gtf`). El merge solo sirve para comparar con gffcompare; los transcritos novedosos no entran en los conteos.

---

## 9. Modificar el pipeline y resolver problemas
 
### Contrato entre scripts
 
Si cambias argumentos en uno, ajusta el que lo llama.
 
| Llamada | Argumentos (en orden) |
|---|---|
| `rnapipe.sh` → `sample_proc.sh` | `sampledir`, `sample`, `PE\|SE`, `R1`, `R2\|NA`, `threads` |
| `rnapipe.sh` → `transcriptome_merge.sh` | `project`, `insdir`, `go_mode`, `padj`, `lfc`, `organism`, `gene_id_type`, `read_length` |
| `transcriptome_merge.sh` → `rnapipe.R` | `project`, `go_mode`, `padj`, `lfc`, `organism`, `gene_id_type` |
 
### Receta de cambios habituales
 
| Quiero… | Dónde tocar |
|---|---|
| Cambiar parámetros de STAR (introns, RAM, etc.) | Bloque `STAR` de `sample_proc.sh`. |
| Cambiar memoria, tiempo o partición de los jobs | Líneas `sbatch` de `rnapipe.sh` (añade `--mem`, `--time`, `--partition`). |
| Añadir un paso por muestra (p. ej. recorte de adaptadores) | `sample_proc.sh`, antes de STAR. Haz que STAR lea los ficheros recortados. |
| Cambiar cómo se cuantifica | Llamada a `stringtie` en `sample_proc.sh`; si cambia el formato de salida, revisa `transcriptome_merge.sh`. |
| Añadir otra especie para GO/KEGG | Tabla `ORGS` al principio de `rnapipe.R` (paquete OrgDb, código KEGG, tipo de ID) y, si quieres detección automática, el bloque `auto`. |
| Añadir covariables al modelo (lote, sexo…) | `design.txt` y la fórmula de `DESeqDataSetFromMatrix` en `rnapipe.R`. |
| Aplicar el umbral de log2FC sobre valores encogidos | Sustituir `results()` por `lfcShrink()` en el bucle de contrastes de `rnapipe.R`. Cambia los genes que pasan el corte. |
| Cambiar el filtro de genes poco expresados | Línea `rowSums(counts(dds)) >= 10` de `rnapipe.R`. |
| Aceptar lecturas sin comprimir | `rnapipe.sh` (comprobación de `.gz` y nombres `.fq.gz`) y `--readFilesCommand` en `sample_proc.sh`. |
| Relanzar solo el análisis en R | `Rscript rnapipe.R <project> <go> <padj> <lfc> <organismo> <tipo_id>` sobre un experimento ya procesado. |
 
Tras tocar un script, haz una prueba con un dataset pequeño (2 condiciones × 2 réplicas, un cromosoma) antes de lanzar el experimento completo.
 
### Problemas frecuentes
 
| Síntoma | Causa probable |
|---|---|
| `Ya existe <ruta>` | El experimento ya fue creado. Cambia `experiment_name` o borra la carpeta. |
| `Hay X condiciones … y Y muestras` | `experimental_design` no coincide con el nº de muestras detectado. Revisa el log de detección. |
| Una muestra se detecta SE y debería ser PE | El nombre no sigue `_R1/_R2` o `_1/_2`, o el R2 no está en la misma carpeta. |
| El merge no arranca | Algún job de muestra falló. Revisa `logs/err_sample_N`. |
| `prepDE.py no está en el PATH` | Carga el módulo de StringTie o añade el script al PATH. |
| `Muestras del diseño ausentes en gene_count_matrix.csv` | Los nombres de `design.txt` no coinciden con los de la matriz. |
| `go_mode debe ser BP, MF, CC o ALL` | Valor inválido en `go_mode`. |
| GO vacío con muchos genes | `gene_id_type` no corresponde a los IDs del GTF (p. ej. IDs Ensembl con `organism ath`). |
| `enrichKEGG falló` | El nodo no tiene internet. |

 ---
---
Siguiente paso previsto: migrar a Snakemake (ver `ROADMAP_SNAKEMAKE.md`).