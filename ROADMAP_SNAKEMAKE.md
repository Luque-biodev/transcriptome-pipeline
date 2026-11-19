# Road to Snakemake

Orden de implementación para migrar `rnapipe` (bash + SLURM) a Snakemake. Cada etapa deja el pipeline funcionando y se da por terminada con un criterio comprobable, no con una fecha. No empieces una etapa sin haber cerrado la anterior.

## Qué se gana

| Hoy | Con Snakemake |
|---|---|
| El orden de ejecución está escrito a mano (`sbatch` + `--dependency`) | Se deduce de entradas y salidas |
| Si falla una muestra hay que relanzar todo | Se reanuda solo lo pendiente |
| Parámetros en un `params.txt` leído con `grep` | `config.yaml` + `samples.tsv` validables |
| Herramientas del PATH del clúster | Entornos por regla con versiones fijas |
| Añadir un paso implica tocar varios scripts | Añadir una regla |

---

## Etapa 0. Línea base

**Objetivo:** poder demostrar que la migración no cambia los resultados.

- Etiquetar la versión actual en git (`v1-bash`).
- Preparar un dataset pequeño (2 condiciones × 2 réplicas; genoma reducido o un cromosoma) y guardar los resultados de la versión bash como referencia.
- Instalar Snakemake (≥ 8) y el plugin de SLURM en un entorno propio.
- Decidir el mecanismo de herramientas: conda o contenedores.

**Hecho cuando:** existe un dataset de prueba con resultados de referencia y `snakemake --version` funciona.

## Etapa 1. Entradas declarativas

**Objetivo:** sustituir `params.txt` por ficheros estructurados.

- `config/samples.tsv` con columnas `sample`, `condition`, `r1`, `r2` (el `results/samples.tsv` actual ya tiene casi esa forma; falta la cabecera).
- `config/config.yaml` con genoma, anotación, condición de control (explícita, ya no "la primera"), organismo, `gene_id_type`, `go_mode`, umbrales, hilos, `sa_index_nbases`, `read_length`.
- La detección de PE por nombre pasa a ser una utilidad opcional que genera `samples.tsv`; el workflow no adivina nada.

```yaml
# config/config.yaml
samples: config/samples.tsv
genome: /ruta/genoma.fa
annotation: /ruta/anotacion.gtf
control: control
organism: ath
gene_id_type: auto
go_mode: BP
padj: 0.05
lfc: 1
threads: 8
sa_index_nbases: 11
read_length: 75
```

**Hecho cuando:** un script de validación comprueba que las rutas existen, que hay ≥ 2 réplicas por condición y que la condición de control está presente.

## Etapa 2. Reglas equivalentes a los scripts actuales

**Objetivo:** reproducir el pipeline actual, regla a regla, sin mejoras.

| Hoy | Regla |
|---|---|
| `gffread` (si la anotación es GFF) | `gff_to_gtf` |
| `STAR --runMode genomeGenerate` | `star_index` |
| `fastqc` | `fastqc` |
| `STAR` + `samtools index` | `star_align`, `bam_index` |
| `stringtie -e` | `stringtie_quant` |
| `stringtie --merge` | `stringtie_merge` |
| `gffcompare` | `gffcompare` |
| `prepDE.py` | `prepde` |
| `rnapipe.R` | `deseq2` |

Orden de escritura recomendado: `star_index` → `star_align` → `stringtie_quant` → `stringtie_merge` → `prepde` → `deseq2` → `fastqc` y `gffcompare` al final (no bloquean el resto).

Piezas que sustituyen a lógica del bash:
- Una función `reads(wildcards)` devuelve 1 o 2 ficheros según `samples.tsv`; la regla de STAR no necesita saber si es PE o SE.
- `mergelist.txt`, `prepde_list.txt` y el contador de muestras desaparecen: se generan a partir de `expand()` sobre las muestras.
- `threads:` y `resources:` sustituyen a `--cpus-per-task` y a la memoria.
- Un `log:` por regla sustituye a la carpeta `logs/`.

```python
import pandas as pd
configfile: "config/config.yaml"
samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("sample", drop=False)
SAMPLES = list(samples.index)

def reads(wc):
    r = samples.loc[wc.sample]
    return [r.r1] + ([r.r2] if isinstance(r.r2, str) and r.r2 else [])

rule star_align:
    input: r=reads, idx="genome/index"
    output: "samples/{sample}/{sample}.bam"
    threads: config["threads"]
    resources: mem_mb=24000
    shell: "STAR ... --readFilesIn {input.r} ... --runThreadN {threads}"
```

**Hecho cuando:** `snakemake -n` muestra el DAG esperado y una ejecución local sobre el dataset de prueba produce `gene_count_matrix.csv`.

## Etapa 3. Ejecución en SLURM

**Objetivo:** que el workflow corra en el clúster con las mismas capacidades que hoy.

- Crear un perfil (`profiles/slurm/config.yaml`) con executor, número máximo de jobs, partición y recursos por defecto. Los nombres exactos de las claves dependen de la versión del plugin: compruébalos en su documentación.
- Dar a cada regla pesada sus `resources:` (memoria y tiempo).

**Hecho cuando:** el dataset de prueba corre en SLURM, y una muestra que falla a propósito no bloquea a las demás ni obliga a repetir lo ya hecho al relanzar.

## Etapa 4. Validación de equivalencia

**Objetivo:** confirmar que Snakemake y la versión bash dan lo mismo.

- Comparar `gene_count_matrix.csv`, `DE_*.csv` y las listas `activated_`/`repressed_` con la línea base de la etapa 0.
- Documentar cualquier diferencia y su causa (versiones de herramientas, orden de muestras).

**Hecho cuando:** los resultados coinciden o cada diferencia está explicada.

A partir de aquí la versión bash puede archivarse.

## Etapa 5. Entornos reproducibles

**Objetivo:** que el resultado no dependa de lo que haya instalado en el clúster.

- Un fichero de entorno por herramienta con versiones fijadas (`star`, `stringtie`, `samtools`, `fastqc`, `gffread`/`gffcompare`, R con DESeq2 y clusterProfiler).
- Para R/Bioconductor, valorar un contenedor: suele resolver mejor las dependencias que conda.
- Precrear los entornos antes de lanzar jobs.

**Hecho cuando:** el workflow corre en una cuenta o nodo sin las herramientas en el PATH.

## Etapa 6. R integrado en el workflow

**Objetivo:** que las salidas de R sean declaradas y el fallo de una parte no invalide el resto.

- Separar en dos reglas: `deseq2` (DE, PCA, volcano) y `enrichment` (GO/KEGG por contraste). Un fallo de red en KEGG no debe tumbar el análisis diferencial.
- Declarar como salidas los ficheros clave. Los contrastes salen de `samples.tsv` (`expand("results/DE_{t}_vs_{c}.csv", ...)`).
- Pasar parámetros desde `config` (con `script:` y `snakemake@params`, o manteniendo `Rscript` con argumentos desde `shell:`; lo segundo es el cambio más pequeño).

**Hecho cuando:** relanzar solo `enrichment` no repite DESeq2.

## Etapa 7. Calidad y trazabilidad

**Objetivo:** detectar problemas y dejar rastro de cada ejecución.

- Validar `config.yaml` y `samples.tsv` con esquemas (`validate()` de `snakemake.utils`).
- Regla `multiqc` sobre FastQC y los logs de STAR.
- Informe `snakemake --report`.
- Prueba automática con el dataset pequeño que compare con los resultados de referencia.
- Guardar `sessionInfo()` de R en `results/`.

**Hecho cuando:** un cambio en cualquier regla se puede comprobar ejecutando una sola prueba.

## Etapa 8. Mejoras opcionales

Se abordan una por una y solo después de la etapa 4, porque cambian los resultados:

| Mejora | Efecto |
|---|---|
| Recortar adaptadores (fastp o cutadapt) antes de STAR | Nueva regla; cambia los alineamientos. |
| Re-cuantificar contra el GTF fusionado (`stringtie -e -B -G merged.gtf`) | Incluye transcritos nuevos en los conteos. |
| `featureCounts` en lugar de `prepDE.py` | Cuantificación más habitual con DESeq2. |
| `lfcShrink` antes de aplicar el umbral de log2FC | Cambia qué genes pasan el corte. |
| Covariables en el modelo (lote, sexo…) | Nueva columna en `samples.tsv` y cambio en la fórmula. |

## Estructura final

```
rnapipe-smk/
├── workflow/
│   ├── Snakefile
│   ├── rules/            # una regla (o grupo) por fichero
│   ├── envs/             # un entorno por herramienta
│   └── scripts/rnapipe.R
├── config/
│   ├── config.yaml
│   └── samples.tsv
├── profiles/slurm/config.yaml
└── README.md
```

## Resumen del orden

| # | Etapa | Depende de | Resultado |
|---|---|---|---|
| 0 | Línea base | — | Dataset de prueba + referencia |
| 1 | Entradas declarativas | 0 | `config.yaml` + `samples.tsv` validados |
| 2 | Reglas equivalentes | 1 | DAG correcto y matriz de conteos |
| 3 | SLURM | 2 | Ejecución en clúster con reanudación |
| 4 | Equivalencia | 3 | Resultados iguales a la versión bash |
| 5 | Entornos | 4 | Independencia del PATH |
| 6 | R integrado | 4 | DE y enriquecimiento separados |
| 7 | Calidad | 5, 6 | Validación, MultiQC, informe, test |
| 8 | Mejoras | 4 | Cambios de método, uno a uno |
