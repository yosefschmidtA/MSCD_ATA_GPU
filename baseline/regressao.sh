#!/bin/bash
# Roda o randmscd_parallel e compara com a linha de base em baseline/.
# Uso: ./baseline/regressao.sh [np]     (np default 6)
#
# Criterio: as linhas de DADOS de saida1Co-alterado-alexandre.txt tem de ser
# identicas byte a byte. Descartadas as linhas de data/duracao, as unicas que
# variam entre rodadas iguais: o rodape ("This calculation took ... starting
# on ...") e o cabecalho ("calculated by <user> on <data>"). O cabecalho entrou
# no filtro em 05/08/2026 -- sem ele a regressao acusa 1 linha de diferenca
# sempre que a base e a rodada caem em dias diferentes.
set -u
NP=${1:-6}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BASE=$ROOT/baseline
OUT=$ROOT/saida1Co-alterado-alexandre.txt
TMP=${TMPDIR:-/tmp}/regr.$$

corpo() { grep -v -E "This calculation took|starting on|and ending on|calculated by" "$1"; }

cd "$ROOT" || exit 1
[ -f "$BASE/saida.txt" ] || { echo "sem linha de base em $BASE"; exit 1; }

/usr/bin/time -f "WALL %e s  MAXRSS %M kB" -o "$TMP.time" \
  mpirun --use-hwthread-cpus -np "$NP" randmscd_parallel Cov0.txt > "$TMP.stdout" 2>&1
rc=$?

echo "--- tempo (np=$NP) ---"; cat "$TMP.time"
echo "--- factors ---"
echo "  agora: $(grep -h 'factors =' "$TMP.stdout" | tail -1)"
echo "  base:  $(grep -h 'factors =' "$BASE/stdout.txt" | tail -1)"
echo "--- curva ---"
if [ $rc -ne 0 ]; then
  echo "FALHOU: exit=$rc"; tail -20 "$TMP.stdout"; exit 1
fi
corpo "$BASE/saida.txt" > "$TMP.a"; corpo "$OUT" > "$TMP.b"
if cmp -s "$TMP.a" "$TMP.b"; then
  echo "IDENTICA byte a byte   OK"
  rm -f "$TMP".*
else
  echo "DIFERE — $(diff "$TMP.a" "$TMP.b" | grep -c '^[<>]') linhas:"
  diff "$TMP.a" "$TMP.b" | head -10
  echo "(saidas completas em $TMP.*)"
  exit 1
fi
