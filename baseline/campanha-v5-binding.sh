#!/bin/bash
# Campanha do V5: separa o ganho do OpenMP no pathcut do custo que o modo de
# binding cobra no laco dos pontos. Os tres modos tem de ser medidos na MESMA
# janela -- a deriva entre janelas ja foi de 6,5% neste projeto, maior que o
# efeito procurado.
#
# Uso: ./baseline/campanha-v5-binding.sh [np] [repeticoes]
#
# Requer binario compilado COM -fopenmp:
#   rm -f *.o && make randmscd_parallel \
#     CPPFLAGS="-O3 -std=c++98 -w -fpermissive -fopenmp"
set -u
NP=${1:-1}
REP=${2:-2}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
CSV=$ROOT/baseline/binding.csv
JANELA=$(date +%Y%m%d-%H%M)

cd "$ROOT" || exit 1

n=$(ps -eo args | grep -cE "^mpirun|^randmscd") || true
if [ "$n" -ne 0 ]; then
  echo "ABORTADO: ha $n processo(s) mpirun/randmscd rodando. Meca com a maquina vazia."
  exit 1
fi
echo "load inicial: $(cat /proc/loadavg)"

corpo() { grep -v -E "This calculation took|starting on|and ending on|calculated by" "$1"; }

[ -f "$CSV" ] || echo "janela,np,modo,rep,wall_s,rss_kb,curva" > "$CSV"

# Rodada de aquecimento, descartada: a primeira sempre paga cache frio de disco.
echo "--- aquecimento (descartado) ---"
mpirun --use-hwthread-cpus -np "$NP" randmscd_parallel Cov0.txt >/dev/null 2>&1

for modo in "padrao:" "bindnone:--bind-to none" "pe:--map-by slot:PE=$(nproc)"; do
  nome=${modo%%:*}
  flags=${modo#*:}
  for r in $(seq 1 "$REP"); do
    sleep 45
    # A saida e' apagada ANTES de cada rodada. Sem isso, uma rodada que falha
    # (ex.: --map-by slot:PE=N nao cabe no np pedido) deixa o arquivo da rodada
    # anterior no lugar e a comparacao devolve "identica" para um run que nunca
    # aconteceu -- falso positivo observado em 05/08/2026, np=12, braco "pe".
    rm -f "$ROOT/saida1Co-alterado-alexandre.txt"
    out=$( { /usr/bin/time -f "%e %M" mpirun --use-hwthread-cpus $flags \
             -np "$NP" randmscd_parallel Cov0.txt >/dev/null; } 2>&1 | tail -1 )
    rc=$?
    wall=${out%% *}
    rss=${out##* }
    if [ $rc -ne 0 ] || [ ! -s "$ROOT/saida1Co-alterado-alexandre.txt" ]; then
      curva=NAO-RODOU
    elif corpo "$ROOT/baseline/saida.txt" > /tmp/v5.a.$$ &&
         corpo "$ROOT/saida1Co-alterado-alexandre.txt" > /tmp/v5.b.$$ &&
         cmp -s /tmp/v5.a.$$ /tmp/v5.b.$$; then curva=identica; else curva=DIFERE; fi
    rm -f /tmp/v5.a.$$ /tmp/v5.b.$$
    printf "%-10s np=%-3s rep%-2s  %8s s  %8s kB  curva %s\n" \
      "$nome" "$NP" "$r" "$wall" "$rss" "$curva"
    echo "$JANELA,$NP,$nome,$r,$wall,$rss,$curva" >> "$CSV"
  done
done

echo
echo "gravado em $CSV (coluna janela: so' compare linhas da mesma janela)"
echo "load final: $(cat /proc/loadavg)"
