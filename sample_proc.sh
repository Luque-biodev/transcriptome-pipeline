#! /bin/bash
#SBATCH --export=ALL
# Uso (lo lanza rnapipe.sh): sample_proc.sh <sampledir> <sample> <PE|SE> <R1> <R2|NA> <threads>
set -euo pipefail

SAMPLEDIR=$1
SAMPLE=$2
TYPE=$3
R1=$4
R2=$5
THREADS=${6:-8}

cd "$SAMPLEDIR"
PROJ=$(cd ../.. && pwd)

if [[ "$TYPE" == "PE" ]]; then READS=("$R1" "$R2"); else READS=("$R1"); fi

fastqc -t "$THREADS" "${READS[@]}"
echo "QC hecho"

echo "Ejecutando STAR en modo $TYPE para $SAMPLE..."
STAR --genomeDir "$PROJ/genome/index" \
     --readFilesIn "${READS[@]}" \
     --readFilesCommand "gunzip -c" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMstrandField intronMotif \
     --outFilterIntronMotifs RemoveNoncanonical \
     --alignIntronMax 100000 \
     --outFileNamePrefix "${SAMPLE}_" \
     --runThreadN "$THREADS" \
     --limitBAMsortRAM 1584388303
echo "Mapeado completo"

BAM="${SAMPLE}_Aligned.sortedByCoord.out.bam"
samtools index "$BAM"
echo "Índice del BAM hecho"

stringtie -e -p "$THREADS" -G "$PROJ/genome/annot.gtf" \
    -o "${SAMPLE}.gtf" -A "${SAMPLE}.tsv" "$BAM"
echo "Cuantificación hecha"
