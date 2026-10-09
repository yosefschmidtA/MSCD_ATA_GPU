#!/bin/bash
# Campanha a frio da Fase 5. Rodar com a maquina livre (sem protetor de tela).
# Alterna referencia e binario final, 3 rodadas cada, e confere bit a bit.
cd /home/madarame/Projetos/MSCDATA || exit 1
n=$(ps -eo args | grep -cE "^(mpirun|prterun|randmscd|gw/)")
[ "$n" -gt 0 ] && { echo "ha outro MSCD rodando, abortando"; exit 1; }
echo "carga inicial: $(cat /proc/loadavg)"
for r in 1 2 3; do
  gw/r.sh gw/ref_gpu cov cref_$r 1 gpu
  gw/r.sh ./randmscd_gpu cov cfin_$r 1 gpu
  gw/r.sh ./randmscd_gpu iron cfin_$r 1 gpu
done
gw/r.sh gw/ref_gpu iron cref_1 1 gpu      # ~7 min, uma vez basta
for f in cov iron; do for t in cfin_1 cfin_2 cfin_3; do
  cmp -s <(grep -vE "calculated by|took|starting|ending" gw/${f}_refgpu.out) \
         <(grep -vE "calculated by|took|starting|ending" gw/${f}_$t.out) \
    && echo "$f $t bit a bit" || echo "$f $t DIFERE"
done; done
echo "carga final: $(cat /proc/loadavg)"
