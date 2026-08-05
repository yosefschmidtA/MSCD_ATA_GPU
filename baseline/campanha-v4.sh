#!/bin/bash
# Campanha de medicao a frio: V2 x V4, intercalados dentro de cada np.
#
#   ./baseline/campanha-v4.sh            # 2 repeticoes, pausa de 45 s
#   REPS=3 PAUSA=60 ./baseline/campanha-v4.sh
#
# APENDA em baseline/escala.csv com janela=<data de hoje>. Nao apaga o CSV:
# as linhas de 2026-08-04 (V0 e V2 da primeira campanha) ficam preservadas.
#
# POR QUE O V2 E' REMEDIDO, se ja existe no CSV:
#   os pontos de V2 de 2026-08-04 foram tomados com binario compilado com
#   -DMSCDTIMER e em outra janela termica. Comparar V4 medido hoje contra
#   aquele V2 e' o erro que o baseline/README.md documenta ("so compare
#   medicoes feitas na mesma janela"). Os dois binarios desta campanha saem
#   da mesma build de producao, sem -DMSCDTIMER.
#   O V0 NAO e' remedido: ele entra na figura como referencia historica, e a
#   comparacao V2-de-hoje contra V2-de-ontem e' o que diz se isso e' legitimo.
#
# As demais precaucoes sao as mesmas do campanha.sh: aborta com a maquina
# ocupada, descarta uma rodada de aquecimento, pausa entre medicoes.
set -u
REPS=${REPS:-2}
PAUSA=${PAUSA:-45}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BASE=$ROOT/baseline
CSV=$BASE/escala.csv
V2BIN=$BASE/randmscd_parallel.v2
V4BIN=$BASE/randmscd_parallel.v4
JANELA=$(date +%Y-%m-%d)
TMP=${TMPDIR:-/tmp}/campanha4.$$

cd "$ROOT" || exit 1
[ -x "$V2BIN" ] || { echo "falta $V2BIN"; exit 1; }
[ -x "$V4BIN" ] || { echo "falta $V4BIN"; exit 1; }
[ -f "$CSV" ] || { echo "falta $CSV"; exit 1; }

# As ancoras ^ sao obrigatorias: sem elas o padrao casa com a linha de comando
# do proprio grep e o script aborta sempre, mesmo com a maquina livre.
ocupada() { ps -eo args | grep -cE "^mpirun|^randmscd"; }
if [ "$(ocupada)" -gt 0 ]; then
  echo "ABORTADO: ja existe mpirun/randmscd rodando. Medicao seria contaminada."
  ps -eo pid,etimes,args | grep -E "mpirun|randmscd" | grep -v grep
  exit 1
fi

# Carga geral nao aborta (navegador aberto e' o caso comum nesta maquina), mas
# fica registrada: e' o que permite julgar a dispersao depois.
CARGA=$(cut -d' ' -f1 /proc/loadavg)
echo "load ao iniciar: $CARGA   janela: $JANELA"

roda() { # $1=binario $2=np  -> imprime "wall factors"
  /usr/bin/time -f "%e" -o "$TMP.t" \
    mpirun --use-hwthread-cpus -np "$2" "$1" Cov0.txt > "$TMP.out" 2> "$TMP.err"
  # tail -1 e nao cat: quando o comando sai com status nao-zero, o GNU time
  # escreve "Command exited with non-zero status N" ANTES do %e.
  echo "$(tail -1 "$TMP.t") $(grep -h 'factors =' "$TMP.out" | tail -1 |
    sed 's/.*factors = *//; s/  */ /g')"
}

echo "aquecimento (descartado)..."
roda "$V4BIN" 6 > /dev/null
sleep "$PAUSA"

for rep in $(seq 1 "$REPS"); do
  for np in 1 2 4 6 8 12; do
    for ver in V2 V4; do
      [ "$ver" = V2 ] && BIN=$V2BIN || BIN=$V4BIN
      out=$(roda "$BIN" "$np")
      w=${out%% *}; f=${out#* }
      echo "$JANELA,$ver,$np,$rep,$w,$f" >> "$CSV"
      printf '%-3s np=%-3s rep=%s  %8ss  [%s]\n' "$ver" "$np" "$rep" "$w" "$f"
      sleep "$PAUSA"
    done
  done
done

rm -f "$TMP".*
echo "--- CSV em $CSV ---"
python3 "$BASE/escala.py"
