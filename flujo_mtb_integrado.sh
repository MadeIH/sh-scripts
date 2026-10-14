#!/usr/bin/env bash
# =============================================================================
# Flujo M. tuberculosis integrado (Illumina pareado)
# Base: run_pipeline_masivo.sh  +  mejoras de flujo_mtb.sh
# SIN Kraken2 (sin base de datos) y SIN genoma de referencia (sin Snippy ni
# QUAST con -r).
#
# Por muestra (mismo orden y mismos parametros que run_pipeline_masivo.sh):
#   1. FastQC (crudas)             -> <proyecto>/qc
#   2. fastp                       -> <proyecto>/trimmed
#   3. FastQC (recortadas)         -> <proyecto>/qc_trim      (opcional)
#   4. SPAdes                      -> <proyecto>/spades_out
#   5. Filtro >=200 pb + QUAST     -> <proyecto>/eval
#   6. Prokka                      -> <proyecto>/prokka_out
#   7. Abricate (resfinder, vfdb)  -> <proyecto>/abricate_out
#   8. TB-Profiler                 -> <proyecto>/tbprofiler   (opcional, apagado)
# Al final: MultiQC (lecturas y ensambles) y resumen de Abricate.
#
# Uso:
#   bash flujo_mtb_integrado.sh PRJEB44165              # lista_limpia.txt si existe;
#                                                       # si no, todos los *_1.fastq.gz
#   bash flujo_mtb_integrado.sh PRJEB44165 mi_lista.txt # lista explicita
#   bash flujo_mtb_integrado.sh /ruta/completa/proyecto
#
# Ajustes por variable de entorno (todos opcionales):
#   HILOS=48 MEM=30            # SPAdes/Prokka; run_pipeline_masivo.sh usaba 48 y 30
#   RUN_TBPROFILER=1           # activa TB-Profiler
#   RUN_FASTQC_TRIM=0 RUN_PROKKA=0 RUN_ABRICATE=0
#   ENV_QC=bioasm ENV_QUAST=bioasm ENV_PROKKA=bioasm   # todo en un solo entorno
#   STOP_ON_ERROR=1            # detenerse en la primera muestra que falle
#   BASE_DIR=/ruta/a/proyectos # donde buscar la carpeta del proyecto
#
# Reanudable: cada paso se omite si su resultado ya existe. Para repetir un
# paso, borra su carpeta/archivo de salida de esa muestra.
#
# -----------------------------------------------------------------------------
# BASES DE DATOS
#   - Kraken2 y referencia H37Rv: NO se usan.
#   - Prokka: base propia incluida en su instalacion.
#   - Abricate: bases incluidas; si es una instalacion nueva ejecuta una vez
#       conda activate prokka_env && abricate --setupdb
#   - TB-Profiler: trae su base de mutaciones (y su referencia) incluida;
#       actualizar con:  tb-profiler update_tbdb
# =============================================================================

set -o pipefail

USO="Uso: bash $(basename "$0") <proyecto> [lista.txt]
  <proyecto>   nombre de carpeta (Ej. PRJEB44165) o ruta completa
  [lista.txt]  accesiones, una por linea (por defecto <proyecto>/lista_limpia.txt;
               si no existe se procesan todos los *_1.fastq.gz)"

if [[ -z "${1:-}" ]]; then
  echo "Error: Debes indicar el nombre de la carpeta (Ej. PRJEB44165)." >&2
  echo "$USO" >&2
  exit 1
fi

# ----- CONFIGURACION ---------------------------------------------------------
HILOS="${HILOS:-20}"                 # SPAdes y Prokka (en masivo: 48 y 30)
MEM="${MEM:-16}"                     # GB maximos para SPAdes (en masivo: 30)
HILOS_FASTQC="${HILOS_FASTQC:-4}"
HILOS_FASTP="${HILOS_FASTP:-4}"
HILOS_QUAST="${HILOS_QUAST:-4}"
MIN_CONTIG="${MIN_CONTIG:-200}"

# Opciones de SPAdes de run_pipeline_masivo.sh. OJO: --careful e --isolate son
# incompatibles entre si; si prefieres --isolate quita --careful.
SPADES_OPTS="${SPADES_OPTS:---careful -k 21,33,55,77}"

GENUS="${GENUS:-Mycobacterium}"
SPECIES="${SPECIES:-tuberculosis}"

ENV_QC="${ENV_QC:-hugo}"             # fastqc, fastp, spades, multiqc
ENV_QUAST="${ENV_QUAST:-quast_env}"  # seqkit, quast
ENV_PROKKA="${ENV_PROKKA:-prokka_env}"  # prokka, abricate
ENV_TBPROFILER="${ENV_TBPROFILER:-tbprofiler}"
ENV_FALLBACK="${ENV_FALLBACK:-bioasm}"  # se usa si un entorno de arriba no existe

RUN_FASTQC_TRIM="${RUN_FASTQC_TRIM:-1}"
RUN_PROKKA="${RUN_PROKKA:-1}"
RUN_ABRICATE="${RUN_ABRICATE:-1}"
RUN_TBPROFILER="${RUN_TBPROFILER:-0}"
STOP_ON_ERROR="${STOP_ON_ERROR:-0}"

# ----- LOCALIZAR EL PROYECTO -------------------------------------------------
BASES=("${BASE_DIR:-}" "$HOME/Desktop/hugo/secuencias" "$HOME/datos")
PROYECTO_DIR=""
if [[ -d $1 ]]; then
  PROYECTO_DIR="$(cd "$1" && pwd)"
else
  for b in "${BASES[@]}"; do
    [[ -n $b && -d $b/$1 ]] && { PROYECTO_DIR="$b/$1"; break; }
  done
fi
if [[ -z $PROYECTO_DIR ]]; then
  echo "Error: la carpeta $1 no existe (se busco tambien en: ${BASES[*]})." >&2
  exit 1
fi

# Crudas: en <proyecto>/fastq si existe (estilo flujo_mtb), si no en <proyecto> (estilo masivo)
if compgen -G "$PROYECTO_DIR/fastq/*_1.fastq.gz" > /dev/null; then
  RAW_DIR="$PROYECTO_DIR/fastq"
else
  RAW_DIR="$PROYECTO_DIR"
fi

# Carpetas de salida DENTRO del proyecto (misma estructura que masivo)
mkdir -p "$PROYECTO_DIR/trimmed" "$PROYECTO_DIR/spades_out" "$PROYECTO_DIR/qc" \
         "$PROYECTO_DIR/qc_trim" "$PROYECTO_DIR/eval" "$PROYECTO_DIR/prokka_out" \
         "$PROYECTO_DIR/abricate_out" "$PROYECTO_DIR/tbprofiler" \
         "$PROYECTO_DIR/multiqc" "$PROYECTO_DIR/logs"

exec > >(tee -a "$PROYECTO_DIR/logs/flujo_mtb.log") 2>&1

echo "=========================================================="
echo " PREPARANDO ENTORNO PARA: $1"
echo " Proyecto : $PROYECTO_DIR"
echo " Crudas   : $RAW_DIR"
echo "=========================================================="

# ----- MUESTRAS --------------------------------------------------------------
if [[ -n ${2:-} ]]; then
  LISTA_FILE="$2"; [[ -f $LISTA_FILE ]] || LISTA_FILE="$PROYECTO_DIR/$2"
  [[ -f $LISTA_FILE ]] || { echo "Error: no existe la lista $2" >&2; exit 1; }
else
  LISTA_FILE="${LISTA:-$PROYECTO_DIR/lista_limpia.txt}"
fi

if [[ -f $LISTA_FILE ]]; then
  mapfile -t ACCS < <(tr -d '\r' < "$LISTA_FILE" | sed '/^[[:space:]]*$/d')
  echo "Muestras: ${#ACCS[@]} (lista: $LISTA_FILE)"
else
  ACCS=()
  for f in "$RAW_DIR"/*_1.fastq.gz; do
    [[ -e $f ]] || { echo "No se encontraron archivos _1.fastq.gz en $RAW_DIR." >&2; exit 1; }
    b="$(basename "$f")"; ACCS+=("${b%_1.fastq.gz}")
  done
  echo "Muestras: ${#ACCS[@]} (todos los *_1.fastq.gz de $RAW_DIR)"
fi
(( ${#ACCS[@]} )) || { echo "La lista de muestras esta vacia." >&2; exit 1; }

# ----- CONDA -----------------------------------------------------------------
if [[ -z ${CONDA_SH:-} ]]; then
  if command -v conda > /dev/null 2>&1; then
    CONDA_SH="$(conda info --base)/etc/profile.d/conda.sh"
  else
    for c in /opt/conda "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/mambaforge"; do
      [[ -f $c/etc/profile.d/conda.sh ]] && { CONDA_SH="$c/etc/profile.d/conda.sh"; break; }
    done
  fi
fi
[[ -f ${CONDA_SH:-} ]] || { echo "Error: no se encontro conda.sh (define CONDA_SH=/ruta/etc/profile.d/conda.sh)." >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONDA_SH"

ENVS_DISP="$(conda env list | awk '!/^#/ && NF {print $1}')"

# Activa el entorno pedido; si no existe usa ENV_FALLBACK
activar() {
  local e="$1"
  if grep -qx "$e" <<< "$ENVS_DISP"; then
    conda activate "$e"
  elif grep -qx "$ENV_FALLBACK" <<< "$ENVS_DISP"; then
    echo "AVISO: el entorno '$e' no existe; se usa '$ENV_FALLBACK'." >&2
    conda activate "$ENV_FALLBACK"
  else
    echo "ERROR: no existe el entorno '$e' ni '$ENV_FALLBACK'." >&2
    return 1
  fi
}

paso() { echo "--- [$1] $2 ---"; }

# Recortadas: acepta el nombre de masivo (.trim.fastq.gz) o el de flujo_mtb (.trim.fq.gz)
buscar_trim() {
  local m="$1" ext
  for ext in trim.fastq.gz trim.fq.gz; do
    if [[ -s $PROYECTO_DIR/trimmed/${m}_1.$ext && -s $PROYECTO_DIR/trimmed/${m}_2.$ext ]]; then
      T1="$PROYECTO_DIR/trimmed/${m}_1.$ext"; T2="$PROYECTO_DIR/trimmed/${m}_2.$ext"
      TRIM_EXISTE=1; return 0
    fi
  done
  T1="$PROYECTO_DIR/trimmed/${m}_1.trim.fastq.gz"; T2="$PROYECTO_DIR/trimmed/${m}_2.trim.fastq.gz"
  TRIM_EXISTE=0
}

# ----- PIPELINE POR MUESTRA --------------------------------------------------
procesar_muestra() {
  local M="$1"
  local R1="$RAW_DIR/${M}_1.fastq.gz"
  local R2="$RAW_DIR/${M}_2.fastq.gz"
  local L="$PROYECTO_DIR/logs"
  local SP="$PROYECTO_DIR/spades_out/$M"
  local FIL="$SP/filtrado.fasta"
  local PREF="M.tuberculosis_${M}"

  [[ -f $R1 && -f $R2 ]] || { echo "SIN CRUDAS PAREADAS: $M" >&2; return 3; }

  # ==========================================
  # AMBIENTE: ENV_QC (hugo)
  # ==========================================
  activar "$ENV_QC"

  paso "1/8" "Ejecutando FastQC"
  if [[ -s $PROYECTO_DIR/qc/${M}_1_fastqc.zip && -s $PROYECTO_DIR/qc/${M}_2_fastqc.zip ]]; then
    echo "   ya existe, se omite"
  else
    fastqc "$R1" "$R2" -o "$PROYECTO_DIR/qc/" -t "$HILOS_FASTQC" > "$L/${M}_fastqc.log" 2>&1
  fi

  paso "2/8" "Ejecutando fastp"
  buscar_trim "$M"
  if (( TRIM_EXISTE )); then
    echo "   ya existe, se omite"
  else
    fastp -i "$R1" \
          -I "$R2" \
          -o "$T1" \
          -O "$T2" \
          --detect_adapter_for_pe \
          --qualified_quality_phred 20 \
          --length_required 50 \
          --n_base_limit 5 \
          --cut_right -W 4 -M 20 -w "$HILOS_FASTP" \
          -j "$PROYECTO_DIR/qc/${M}_fastp.json" \
          -h "$PROYECTO_DIR/qc/${M}_fastp.html" > "$L/${M}_fastp.log" 2>&1
  fi

  if (( RUN_FASTQC_TRIM )); then
    paso "3/8" "Ejecutando FastQC sobre lecturas recortadas"
    if [[ -s $PROYECTO_DIR/qc_trim/${M}_1.trim_fastqc.zip && -s $PROYECTO_DIR/qc_trim/${M}_2.trim_fastqc.zip ]]; then
      echo "   ya existe, se omite"
    else
      fastqc "$T1" "$T2" -o "$PROYECTO_DIR/qc_trim/" -t "$HILOS_FASTQC" > "$L/${M}_fastqc_trim.log" 2>&1
    fi
  fi

  paso "4/8" "Ejecutando SPAdes"
  if [[ -s $SP/contigs.fasta ]]; then
    echo "   ya existe, se omite"
  else
    rm -rf "$SP"        # evita restos de una corrida interrumpida
    # shellcheck disable=SC2086
    spades.py -1 "$T1" \
              -2 "$T2" \
              -o "$SP" \
              -t "$HILOS" \
              -m "$MEM" \
              $SPADES_OPTS > "$L/${M}_spades.log" 2>&1
  fi

  conda deactivate

  # ==========================================
  # AMBIENTE: ENV_QUAST (quast_env)
  # ==========================================
  activar "$ENV_QUAST"

  paso "5/8" "Filtrando contigs >= ${MIN_CONTIG} pb con SeqKit y evaluando con QUAST (sin referencia)"
  if [[ -s $FIL ]]; then
    echo "   filtrado ya existe, se omite"
  else
    seqkit seq -m "$MIN_CONTIG" "$SP/contigs.fasta" > "$FIL"
  fi

  if [[ -s $PROYECTO_DIR/eval/quast_${M}/report.tsv ]]; then
    echo "   QUAST ya existe, se omite"
  else
    quast.py "$FIL" \
             -o "$PROYECTO_DIR/eval/quast_${M}" \
             -t "$HILOS_QUAST" \
             --labels "$M" \
             --min-contig "$MIN_CONTIG" > "$L/${M}_quast.log" 2>&1
  fi

  conda deactivate

  # ==========================================
  # AMBIENTE: ENV_PROKKA (prokka_env)
  # ==========================================
  if (( RUN_PROKKA || RUN_ABRICATE )); then
    activar "$ENV_PROKKA"

    if (( RUN_PROKKA )); then
      paso "6/8" "Anotando con Prokka"
      if [[ -s $PROYECTO_DIR/prokka_out/${M}/${PREF}.txt ]]; then
        echo "   ya existe, se omite"
      else
        prokka "$FIL" \
               --outdir "$PROYECTO_DIR/prokka_out/${M}" \
               --prefix "$PREF" \
               --genus "$GENUS" \
               --species "$SPECIES" \
               --cpus "$HILOS" \
               --force > "$L/${M}_prokka.log" 2>&1
      fi
    fi

    if (( RUN_ABRICATE )); then
      paso "7/8" "Analizando resistencia y virulencia con Abricate"
      if [[ ! -s $PROYECTO_DIR/abricate_out/${M}_resistencia.tab ]]; then
        abricate --db resfinder "$FIL" > "$PROYECTO_DIR/abricate_out/${M}_resistencia.tab"
      fi
      if [[ ! -s $PROYECTO_DIR/abricate_out/${M}_virulencia.tab ]]; then
        abricate --db vfdb "$FIL" > "$PROYECTO_DIR/abricate_out/${M}_virulencia.tab"
      fi
    fi

    conda deactivate
  fi

  # ==========================================
  # AMBIENTE: ENV_TBPROFILER (opcional)
  # Usa su propia base y referencia incluidas; no requiere nada externo.
  # ==========================================
  if (( RUN_TBPROFILER )); then
    paso "8/8" "TB-Profiler (linaje y resistencia)"
    if ! grep -qx "$ENV_TBPROFILER" <<< "$ENVS_DISP"; then
      echo "AVISO: no existe el entorno '$ENV_TBPROFILER'; se omite TB-Profiler." >&2
    elif [[ -s $PROYECTO_DIR/tbprofiler/results/${M}.results.json ]]; then
      echo "   ya existe, se omite"
    else
      conda activate "$ENV_TBPROFILER"
      ( cd "$PROYECTO_DIR/tbprofiler" && \
        tb-profiler profile -1 "$T1" -2 "$T2" -p "$M" \
          --threads "$HILOS" --txt > "$L/${M}_tbprofiler.log" 2>&1 )
      conda deactivate
    fi
  fi
}

# ----- BUCLE PRINCIPAL -------------------------------------------------------
OK=(); OMITIDAS=(); FALLIDAS=()
for M in "${ACCS[@]}"; do
  echo "=========================================================="
  echo "   INICIANDO PIPELINE PARA LA CEPA: ${M}"
  echo "=========================================================="

  # Subshell con set -e: si un paso falla, se corta ESA muestra y se sigue con
  # la siguiente (o se detiene todo con STOP_ON_ERROR=1).
  (
    set -eE
    # (el codigo 3 = muestra sin crudas pareadas: no es un fallo, se informa aparte)
    trap 's=$?; (( s == 3 )) || echo "FALLO en la muestra $M (linea $LINENO). Revisa $PROYECTO_DIR/logs/" >&2' ERR
    procesar_muestra "$M"
  )
  rc=$?

  case $rc in
    0) OK+=("$M") ;;
    3) OMITIDAS+=("$M") ;;
    *) FALLIDAS+=("$M"); (( STOP_ON_ERROR )) && { echo "STOP_ON_ERROR=1: se detiene el flujo." >&2; break; } ;;
  esac
done

# ----- RESUMENES FINALES -----------------------------------------------------
echo "=========================================================="
echo " RESUMENES (MultiQC, Abricate, TB-Profiler)"
echo "=========================================================="

# MultiQC solo con las muestras de la lista
activar "$ENV_QC"
command -v multiqc > /dev/null 2>&1 || { conda deactivate; activar "$ENV_QUAST"; }
if command -v multiqc > /dev/null 2>&1; then
  MQ_LEC=(); MQ_ENS=()
  for M in "${ACCS[@]}"; do
    for f in "$PROYECTO_DIR/qc/${M}"_*_fastqc.zip \
             "$PROYECTO_DIR/qc/${M}_fastp.json" \
             "$PROYECTO_DIR/qc_trim/${M}"_*_fastqc.zip; do
      [[ -e $f ]] && MQ_LEC+=("$f")
    done
    [[ -d $PROYECTO_DIR/eval/quast_${M} ]]   && MQ_ENS+=("$PROYECTO_DIR/eval/quast_${M}")
    [[ -d $PROYECTO_DIR/prokka_out/${M} ]]   && MQ_ENS+=("$PROYECTO_DIR/prokka_out/${M}")
  done
  (( ${#MQ_LEC[@]} )) && multiqc "${MQ_LEC[@]}" -n multiqc_lecturas  -o "$PROYECTO_DIR/multiqc" --force > "$PROYECTO_DIR/logs/multiqc_lecturas.log" 2>&1
  (( ${#MQ_ENS[@]} )) && multiqc "${MQ_ENS[@]}" -n multiqc_ensambles -o "$PROYECTO_DIR/multiqc" --force > "$PROYECTO_DIR/logs/multiqc_ensambles.log" 2>&1
else
  echo "AVISO: multiqc no esta instalado en $ENV_QC ni en $ENV_QUAST; se omite." >&2
fi
conda deactivate

# Tabla resumen de Abricate (solo muestras de la lista)
if (( RUN_ABRICATE )); then
  activar "$ENV_PROKKA"
  if command -v abricate > /dev/null 2>&1; then
    AB_RES=(); AB_VIR=()
    for M in "${ACCS[@]}"; do
      [[ -s $PROYECTO_DIR/abricate_out/${M}_resistencia.tab ]] && AB_RES+=("$PROYECTO_DIR/abricate_out/${M}_resistencia.tab")
      [[ -s $PROYECTO_DIR/abricate_out/${M}_virulencia.tab ]]  && AB_VIR+=("$PROYECTO_DIR/abricate_out/${M}_virulencia.tab")
    done
    (( ${#AB_RES[@]} )) && abricate --summary "${AB_RES[@]}" > "$PROYECTO_DIR/abricate_out/resumen_resistencia.tab"
    (( ${#AB_VIR[@]} )) && abricate --summary "${AB_VIR[@]}" > "$PROYECTO_DIR/abricate_out/resumen_virulencia.tab"
  fi
  conda deactivate
fi

# Tabla resumen de TB-Profiler
if (( RUN_TBPROFILER )) && grep -qx "$ENV_TBPROFILER" <<< "$ENVS_DISP" \
   && compgen -G "$PROYECTO_DIR/tbprofiler/results/*.results.json" > /dev/null; then
  conda activate "$ENV_TBPROFILER"
  ( cd "$PROYECTO_DIR/tbprofiler" && tb-profiler collate > "$PROYECTO_DIR/logs/tbprofiler_collate.log" 2>&1 )
  conda deactivate
fi

# ----- CIERRE ----------------------------------------------------------------
if (( ${#FALLIDAS[@]} )); then
  printf '%s\n' "${FALLIDAS[@]}" > "$PROYECTO_DIR/logs/muestras_fallidas.txt"
fi

echo "=========================================================="
echo " TODOS LOS ARCHIVOS DE $1 HAN SIDO PROCESADOS"
echo " OK: ${#OK[@]} | sin crudas pareadas: ${#OMITIDAS[@]} | fallidas: ${#FALLIDAS[@]}"
(( ${#OMITIDAS[@]} )) && echo " Sin crudas : ${OMITIDAS[*]}"
(( ${#FALLIDAS[@]} )) && echo " Fallidas   : ${FALLIDAS[*]}  (lista en logs/muestras_fallidas.txt)"
echo " Los resultados se guardaron dentro de: $PROYECTO_DIR"
echo "=========================================================="

exit $(( ${#FALLIDAS[@]} > 0 ))
