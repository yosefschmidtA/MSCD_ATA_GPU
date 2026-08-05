#!/bin/bash
# Campanha de medicao a frio: V0 x V2 em varios np.
#
#   ./baseline/campanha.sh            # 2 repeticoes, pausa de 45 s
#   REPS=3 PAUSA=60 ./baseline/campanha.sh
#
# Grava baseline/escala.csv e gera baseline/escalabilidade.png no fim.
# Duracao aproximada: ~27 min por repeticao.
#
# POR QUE AS PRECAUCOES (todas custaram tempo real nesta campanha):
#  - maquina ocupada contamina a medicao: o script ABORTA se achar outro
#    mpirun/randmscd rodando. (A rodada de 897 s que motivou esta checagem NAO
#    era contencao -- era o bug de laco infinito do phase.cpp:320, ver
#    OTIMIZACAO.md. A checagem continua valendo; o motivo registrado estava
#    errado. Se um np travar de novo, amostre a pilha antes de matar.)
#  - deriva sob carga continua: depois de ~1,5 h medindo, a mesma versao sem
#    mudanca de codigo passou de 48,4 s para 56,5 s. Dai a pausa entre rodadas.
#  - primeira rodada fria: a primeira execucao apos o boot/build custa ~40% a
#    mais (100,45 s contra 70 s). Dai a rodada de aquecimento descartada.
set -u
REPS=${REPS:-2}
PAUSA=${PAUSA:-45}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BASE=$ROOT/baseline
CSV=$BASE/escala.csv
V0BIN=$BASE/randmscd_parallel.baseline
V2BIN=$ROOT/randmscd_parallel
TMP=${TMPDIR:-/tmp}/campanha.$$

cd "$ROOT" || exit 1
[ -x "$V0BIN" ] || { echo "falta $V0BIN (binario V0 preservado)"; exit 1; }
[ -x "$V2BIN" ] || { echo "falta $V2BIN — rode o make antes"; exit 1; }

# As ancoras ^ sao obrigatorias: sem elas o padrao casa com a linha de comando
# do proprio grep e o script aborta sempre, mesmo com a maquina livre.
ocupada() { ps -eo args | grep -cE "^mpirun|^randmscd"; }
if [ "$(ocupada)" -gt 0 ]; then
  echo "ABORTADO: ja existe mpirun/randmscd rodando. Medicao seria contaminada."
  ps -eo pid,etimes,args | grep -E "mpirun|randmscd" | grep -v grep
  exit 1
fi

roda() { # $1=binario $2=np  -> imprime "wall factors"
  /usr/bin/time -f "%e" -o "$TMP.t" \
    mpirun --use-hwthread-cpus -np "$2" "$1" Cov0.txt > "$TMP.out" 2> "$TMP.err"
  # tail -1 e nao cat: quando o comando sai com status nao-zero, o GNU time
  # escreve "Command exited with non-zero status N" ANTES do %e. Com cat, o
  # campo wall virava a palavra "Command" e o float() do escala.py quebrava,
  # corrompendo o CSV inteiro da campanha.
  echo "$(tail -1 "$TMP.t") $(grep -h 'factors =' "$TMP.out" | tail -1 |
    sed 's/.*factors = *//; s/  */ /g')"
}

echo "aquecimento (descartado)..."
roda "$V2BIN" 6 > /dev/null
sleep "$PAUSA"

echo "versao,np,rep,wall,factors" > "$CSV"
for rep in $(seq 1 "$REPS"); do
  for np in 1 2 4 6 8 12; do
    for ver in V0 V2; do
      [ "$ver" = V0 ] && BIN=$V0BIN || BIN=$V2BIN
      out=$(roda "$BIN" "$np")
      w=${out%% *}; f=${out#* }
      echo "$ver,$np,$rep,$w,$f" >> "$CSV"
      printf '%-3s np=%-3s rep=%s  %8ss  [%s]\n' "$ver" "$np" "$rep" "$w" "$f"
      sleep "$PAUSA"
    done
  done
done

rm -f "$TMP".*
echo "--- CSV em $CSV ---"
python3 "$BASE/escala.py"
