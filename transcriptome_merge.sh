#! /bin/bash
#SBATCH --export=ALL
# Uso (lo lanza rnapipe.sh): transcriptome_merge.sh <project> <insdir> <go> <padj> <lfc> <organism> <geneid> <readlen>
set -euo pipefail

PROJECT=$1
INSDIR=$2
GO_MODE=$3
PVALUE_CUT=$4
LFC_CUT=$5
ORGANISM=$6
GENEID=$7
READLEN=${8:-75}

cd "$PROJECT"

stringtie --merge -G genome/annot.gtf -o results/stringtie_merged.gtf results/mergelist.txt

cd results
gffcompare -r ../genome/annot.gtf -o merged stringtie_merged.gtf

PREPDE=$(command -v prepDE.py || command -v prepDE.py3 || true)
[ -n "$PREPDE" ] || { echo "ERROR: prepDE.py no está en el PATH" >&2; exit 1; }
"$PREPDE" -i prepde_list.txt -l "$READLEN"   # genera gene_count_matrix.csv y transcript_count_matrix.csv

echo "==============================="
echo "Parte de bash terminada, toca R"
echo "==============================="

Rscript "$INSDIR/rnapipe.R" "$PROJECT" "$GO_MODE" "$PVALUE_CUT" "$LFC_CUT" "$ORGANISM" "$GENEID"
