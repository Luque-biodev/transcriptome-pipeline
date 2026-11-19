#! /bin/bash
# rnapipe.sh - Pipeline de transcriptómica comparativa: STAR + StringTie + DESeq2 (+ GO/KEGG)
# Uso: rnapipe.sh <params.file>      (ejemplo en params.txt)
set -eo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

if [ $# -ne 1 ] || [ ! -f "$1" ]; then
    echo "Se usa: rnapipe.sh <params.file>"
    echo "Un ejemplo de params.file se facilita en params.txt"
    exit 1
fi
PARAMS=$(realpath "$1")

# Lee "clave valor..." (ignora comentarios y CRLF). Devuelve todo lo que hay tras la clave.
getp() {
    tr -d '\r' < "$PARAMS" | grep -m1 -E "^$1[[:space:]]" \
        | sed -E 's/^[^[:space:]]+[[:space:]]+//; s/[[:space:]]+#.*$//; s/[[:space:]]+$//' || true
}

WD=$(getp working_dir)
EXP=$(getp experiment_name)
GENOME=$(getp genome_path)
ANNOT=$(getp annot_path)
INSDIR=$(getp installation_dir)
DESIGN=$(getp experimental_design)
SAMSDIR=$(getp samples_path)
AUTHOR=$(getp author_name);     AUTHOR=${AUTHOR:-${USER:-$(whoami)}}
ORGANISM=$(getp organism);      ORGANISM=${ORGANISM:-auto}
GENEID=$(getp gene_id_type);    GENEID=${GENEID:-auto}
GO_MODE=$(getp go_mode);        GO_MODE=${GO_MODE:-BP}
PVALUE_CUT=$(getp pvalue_cutoff); PVALUE_CUT=${PVALUE_CUT:-0.05}
LFC_CUT=$(getp lfc_cutoff);     LFC_CUT=${LFC_CUT:-1}
THREADS=$(getp threads);        THREADS=${THREADS:-8}
SAIDX=$(getp sa_index_nbases);  SAIDX=${SAIDX:-11}
READLEN=$(getp read_length);    READLEN=${READLEN:-75}

for v in WD EXP GENOME ANNOT INSDIR DESIGN; do
    [ -n "${!v}" ] || die "Falta el parámetro obligatorio para $v en $PARAMS"
done
[ -f "$GENOME" ] || die "No existe el genoma: $GENOME"
[ -f "$ANNOT" ]  || die "No existe la anotación: $ANNOT"
INSDIR=$(realpath "$INSDIR")
for f in sample_proc.sh transcriptome_merge.sh rnapipe.R; do
    [ -f "$INSDIR/$f" ] || die "No se encuentra $INSDIR/$f"
done

# ---------------------------------------------------------------- MUESTRAS
# Dos modos: lista explícita (path_sample_1, path_sample_2, ...) o carpeta (samples_path).
# Paired-end: se detecta el mate R2 a partir de _R1/_R2 o _1/_2 en el nombre.
R1S=(); R2S=()
RE_R2a='_R2[._]'
RE_R2b='_2\.(fq|fastq)'

add_sample() {
    local r1=$1 d b b2 r2=""
    d=$(dirname "$r1"); b=$(basename "$r1")
    b2=$(echo "$b" | sed -E 's/_R1([._])/_R2\1/; t; s/_1\.(fq|fastq)/_2.\1/')
    if [ "$b2" != "$b" ] && [ -f "$d/$b2" ]; then r2="$d/$b2"; fi
    R1S+=("$r1"); R2S+=("$r2")
}

mapfile -t PSAMS < <(tr -d '\r' < "$PARAMS" | grep -E '^path_sample_' | sort -V -k1,1 | awk '{print $2}')
if [ ${#PSAMS[@]} -gt 0 ]; then
    FILES=("${PSAMS[@]}")
elif [ -n "$SAMSDIR" ] && [ -d "$SAMSDIR" ]; then
    FILES=("$SAMSDIR"/*)
else
    die "Define path_sample_N en el params o una carpeta válida en samples_path"
fi

for f in "${FILES[@]}"; do
    [ -f "$f" ] || die "No existe la muestra: $f"
    [[ "$f" == *.gz ]] || die "Las muestras deben estar comprimidas (.gz): $f"
    b=$(basename "$f")
    if [[ "$b" =~ $RE_R2a || "$b" =~ $RE_R2b ]]; then continue; fi   # los R2 se emparejan desde R1
    add_sample "$f"
done
NSAM=${#R1S[@]}
[ "$NSAM" -gt 0 ] || die "No se ha detectado ninguna muestra"
echo "Detectadas $NSAM muestras."

# ---------------------------------------------------------------- DISEÑO
read -r -a CONDITIONS <<< "$DESIGN"
[ ${#CONDITIONS[@]} -eq "$NSAM" ] || die "Hay ${#CONDITIONS[@]} condiciones en experimental_design y $NSAM muestras."

declare -A COUNT
for c in "${CONDITIONS[@]}"; do COUNT[$c]=$(( ${COUNT[$c]:-0} + 1 )); done
[ ${#COUNT[@]} -ge 2 ] || die "Se necesitan al menos 2 condiciones."
MINR=999999; MAXR=0
for c in "${!COUNT[@]}"; do
    [ "${COUNT[$c]}" -ge 2 ] || die "La condición '$c' tiene una sola réplica; DESeq2 necesita al menos 2."
    [ "${COUNT[$c]}" -lt $MINR ] && MINR=${COUNT[$c]}
    [ "${COUNT[$c]}" -gt $MAXR ] && MAXR=${COUNT[$c]}
done
[ "$MINR" -eq "$MAXR" ] || echo "AVISO: número desigual de réplicas por condición (min $MINR, max $MAXR)."
echo "La primera condición (${CONDITIONS[0]}) se usará como control."

# ---------------------------------------------------------------- ESPACIO DE TRABAJO
echo "=========================="
echo "CREANDO ESPACIO DE TRABAJO"
echo "=========================="

WD=$(realpath -m "$WD")
PROJECT="$WD/$EXP"
[ ! -e "$PROJECT" ] || die "Ya existe $PROJECT"
mkdir -p "$PROJECT"/{samples,genome,results,logs}
: > "$PROJECT/results/mergelist.txt"
: > "$PROJECT/results/prepde_list.txt"
: > "$PROJECT/results/design.txt"
: > "$PROJECT/results/samples.tsv"

TYPES=()
for (( i=1; i<=NSAM; i++ )); do
    S="sample_$i"
    SD="$PROJECT/samples/$S"
    mkdir -p "$SD"
    r1="${R1S[$((i-1))]}"; r2="${R2S[$((i-1))]}"
    cond="${CONDITIONS[$((i-1))]}"

    cp "$r1" "$SD/${S}_R1.fq.gz"
    if [ -n "$r2" ]; then
        cp "$r2" "$SD/${S}_R2.fq.gz"
        TYPES+=("PE"); echo "Muestra $i: paired-end ($cond)"
    else
        TYPES+=("SE"); echo "Muestra $i: single-end ($cond)"
    fi

    printf "%s\t%s\t%s\t%s\t%s\n" "$S" "${TYPES[$((i-1))]}" "$cond" "$r1" "${r2:-NA}" >> "$PROJECT/results/samples.tsv"
    printf "%s\t%s\n" "$S" "$cond" >> "$PROJECT/results/design.txt"
    echo "$SD/$S.gtf" >> "$PROJECT/results/mergelist.txt"
    printf "%s\t%s\n" "$S" "$SD/$S.gtf" >> "$PROJECT/results/prepde_list.txt"
done

cd "$PROJECT/genome"
mkdir index
cp "$GENOME" genome.fa
case "${ANNOT##*.}" in
    gff|gff3)   # STAR/StringTie necesitan GTF
        cp "$ANNOT" annot.gff
        gffread annot.gff -T -o annot.gtf ;;
    *)  cp "$ANNOT" annot.gtf ;;
esac

echo "========================="
echo "ESPACIO DE TRABAJO CREADO"
echo "========================="

echo "=============="
echo "CREANDO ÍNDICE"
echo "=============="
STAR --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir index \
     --genomeFastaFiles genome.fa --sjdbGTFfile annot.gtf --genomeSAindexNbases "$SAIDX"
echo "============="
echo "ÍNDICE CREADO"
echo "============="

# ---------------------------------------------------------------- JOBS
echo "======================="
echo "PROCESANDO LAS MUESTRAS"
echo "======================="

JOBIDS=()
for (( i=1; i<=NSAM; i++ )); do
    S="sample_$i"
    if [ "${TYPES[$((i-1))]}" == "PE" ]; then R2F="${S}_R2.fq.gz"; else R2F="NA"; fi
    jid=$(sbatch --parsable --job-name="${S}_$AUTHOR" --cpus-per-task="$THREADS" \
          --output="$PROJECT/logs/out_$S" --error="$PROJECT/logs/err_$S" \
          "$INSDIR/sample_proc.sh" "$PROJECT/samples/$S" "$S" "${TYPES[$((i-1))]}" \
          "${S}_R1.fq.gz" "$R2F" "$THREADS")
    JOBIDS+=("${jid%%;*}")
    echo "Enviado: $S (job ${jid%%;*})"
done

# Merge + prepDE + R arrancan solo cuando TODAS las muestras terminan bien
DEPS=$(IFS=:; echo "${JOBIDS[*]}")
jid=$(sbatch --parsable --job-name="tm_$AUTHOR" --dependency="afterok:$DEPS" --kill-on-invalid-dep=yes \
      --output="$PROJECT/logs/out_tm" --error="$PROJECT/logs/err_tm" \
      "$INSDIR/transcriptome_merge.sh" "$PROJECT" "$INSDIR" "$GO_MODE" "$PVALUE_CUT" "$LFC_CUT" \
      "$ORGANISM" "$GENEID" "$READLEN")
echo "Enviado: merge + DESeq2 (job ${jid%%;*}), depende de: $DEPS"
echo "Logs en $PROJECT/logs  |  Resultados en $PROJECT/results"
