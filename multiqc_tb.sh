set -o pipefail

USO="Uso: bash $(basename "$0") <proyecto> [lista.txt]
  <proyecto>   nombre de carpeta (Ej. PRJEB44165) o ruta completa
  [lista.txt]  accesiones, una por linea (por defecto <proyecto>/lista_limpia.txt;
               si no existe se usan las muestras halladas en <proyecto>/qc y eval)"

if [[ -z "${1:-}" ]]; then
  echo "Error: Debes indicar el nombre de la carpeta (Ej. PRJEB44165)." >&2
  echo "$USO" >&2
  exit 1
fi

# ----- CONFIGURACION (mismos nombres que flujo_mtb_integrado.sh) --------------
ENV_QC="${ENV_QC:-hugo}"
ENV_QUAST="${ENV_QUAST:-quast_env}"
ENV_FALLBACK="${ENV_FALLBACK:-hugo}"
MQ_TODO="${MQ_TODO:-0}"
MQ_OPTS="${MQ_OPTS:-}"

# ----- LOCALIZAR EL PROYECTO --------------------------------------------------
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

mkdir -p "$PROYECTO_DIR/multiqc" "$PROYECTO_DIR/logs"
exec > >(tee -a "$PROYECTO_DIR/logs/multiqc_solo.log") 2>&1

echo "=========================================================="
echo " MULTIQC PARA: $1"
echo " Proyecto : $PROYECTO_DIR"
echo "=========================================================="

# ----- MUESTRAS ---------------------------------------------------------------
if [[ -n ${2:-} ]]; then
  LISTA_FILE="$2"; [[ -f $LISTA_FILE ]] || LISTA_FILE="$PROYECTO_DIR/$2"
  [[ -f $LISTA_FILE ]] || { echo "Error: no existe la lista $2" >&2; exit 1; }
else
  LISTA_FILE="$PROYECTO_DIR/clearlist.txt"
fi

if [[ -f $LISTA_FILE ]]; then
  mapfile -t ACCS < <(tr -d '\r' < "$LISTA_FILE" | sed '/^[[:space:]]*$/d')
  echo "Muestras: ${#ACCS[@]} (lista: $LISTA_FILE)"
else
  # Sin lista: se deducen de los resultados existentes (fastp .json y carpetas QUAST)
  ACCS=()
  for f in "$PROYECTO_DIR"/qc/*_fastp.json; do
    [[ -e $f ]] || continue
    b="$(basename "$f")"; ACCS+=("${b%_fastp.json}")
  done
  for d in "$PROYECTO_DIR"/eval/quast_*/; do
    [[ -d $d ]] || continue
    b="$(basename "$d")"; ACCS+=("${b#quast_}")
  done
  if (( ${#ACCS[@]} )); then
    mapfile -t ACCS < <(printf '%s\n' "${ACCS[@]}" | sort -u)
  fi
  echo "Muestras: ${#ACCS[@]} (deducidas de $PROYECTO_DIR/qc y /eval)"
fi
(( ${#ACCS[@]} )) || { echo "No hay muestras ni resultados para resumir." >&2; exit 1; }

# ----- CONDA ------------------------------------------------------------------
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

# multiqc: primero ENV_QC; si no lo tiene, ENV_QUAST
activar "$ENV_QC"
command -v multiqc > /dev/null 2>&1 || { conda deactivate; activar "$ENV_QUAST"; }
if ! command -v multiqc > /dev/null 2>&1; then
  echo "ERROR: multiqc no esta instalado en $ENV_QC ni en $ENV_QUAST." >&2
  echo "       Instalalo con: conda install -n $ENV_QC -c bioconda multiqc" >&2
  conda deactivate
  exit 1
fi
echo "MultiQC: $(multiqc --version)"

# ----- RECOLECTAR ARCHIVOS (solo de las muestras de la lista) ------------------
MQ_LEC=(); MQ_ENS=(); SIN_LEC=(); SIN_ENS=()
for M in "${ACCS[@]}"; do
  n_lec=${#MQ_LEC[@]}; n_ens=${#MQ_ENS[@]}

  for f in "$PROYECTO_DIR/qc/${M}_1_fastqc.zip" \
           "$PROYECTO_DIR/qc/${M}_2_fastqc.zip" \
           "$PROYECTO_DIR/qc/${M}_fastp.json" \
           "$PROYECTO_DIR/qc_trim/${M}_1.trim_fastqc.zip" \
           "$PROYECTO_DIR/qc_trim/${M}_2.trim_fastqc.zip"; do
    [[ -e $f ]] && MQ_LEC+=("$f")
  done
  [[ -s $PROYECTO_DIR/eval/quast_${M}/report.tsv ]] && MQ_ENS+=("$PROYECTO_DIR/eval/quast_${M}")
  [[ -d $PROYECTO_DIR/prokka_out/${M} ]]            && MQ_ENS+=("$PROYECTO_DIR/prokka_out/${M}")

  (( ${#MQ_LEC[@]} == n_lec )) && SIN_LEC+=("$M")
  (( ${#MQ_ENS[@]} == n_ens )) && SIN_ENS+=("$M")
done

echo "Archivos de lecturas : ${#MQ_LEC[@]}"
echo "Carpetas de ensambles: ${#MQ_ENS[@]}"
(( ${#SIN_LEC[@]} )) && echo "AVISO: sin resultados de lecturas : ${SIN_LEC[*]}" >&2
(( ${#SIN_ENS[@]} )) && echo "AVISO: sin resultados de ensambles: ${SIN_ENS[*]}" >&2

# ----- MULTIQC ----------------------------------------------------------------
OUT="$PROYECTO_DIR/multiqc"
RC=0

# shellcheck disable=SC2086
if (( ${#MQ_LEC[@]} )); then
  echo "--- MultiQC: lecturas ---"
  multiqc "${MQ_LEC[@]}" -n multiqc_lecturas -o "$OUT" --force $MQ_OPTS \
    > "$PROYECTO_DIR/logs/multiqc_lecturas.log" 2>&1 || { echo "FALLO multiqc_lecturas (ver logs/)" >&2; RC=1; }
else
  echo "AVISO: no hay resultados de lecturas; se omite multiqc_lecturas." >&2
fi

# shellcheck disable=SC2086
if (( ${#MQ_ENS[@]} )); then
  echo "--- MultiQC: ensambles ---"
  multiqc "${MQ_ENS[@]}" -n multiqc_ensambles -o "$OUT" --force $MQ_OPTS \
    > "$PROYECTO_DIR/logs/multiqc_ensambles.log" 2>&1 || { echo "FALLO multiqc_ensambles (ver logs/)" >&2; RC=1; }
else
  echo "AVISO: no hay resultados de ensambles; se omite multiqc_ensambles." >&2
fi

# shellcheck disable=SC2086
if (( MQ_TODO )) && (( ${#MQ_LEC[@]} + ${#MQ_ENS[@]} )); then
  echo "--- MultiQC: todo junto ---"
  multiqc "${MQ_LEC[@]}" "${MQ_ENS[@]}" -n multiqc_todo -o "$OUT" --force $MQ_OPTS \
    > "$PROYECTO_DIR/logs/multiqc_todo.log" 2>&1 || { echo "FALLO multiqc_todo (ver logs/)" >&2; RC=1; }
fi

conda deactivate

echo "=========================================================="
echo " LISTO. Reportes en: $OUT"
ls -1 "$OUT"/*.html 2>/dev/null | sed 's/^/   /'
echo "=========================================================="

exit $RC
