#!/bin/bash
# =====================================================================
# Bateria de desempenho do MSCD: executavel ORIGINAL x versao com GPU.
#
#   ./bateria/bateria.sh            roda (ou retoma) a bateria
#   ./bateria/bateria.sh resumo     mostra a tabela do que ja' foi medido
#   touch bateria/PAUSAR            pausa depois da rodada atual
#   Ctrl+C                          para na hora (a rodada atual e' refeita)
#
# Cada rodada terminada vai para bateria/resultados.csv. Ao rodar de novo, o
# que ja' esta' no CSV e' pulado, entao pode desligar o PC a qualquer momento:
# perde-se so' a rodada que estava em andamento.
#
# Grade: Cov0.txt com raio/profundidade crescentes (profundidade = 1,5 x
# raio) e -np de 1 a 12, primeiro o original e depois a GPU.
# Tudo pode ser mudado por variavel de ambiente, por exemplo
#   NPS="1 2 4" TEMPO_MAX=3600 ./bateria/bateria.sh
# =====================================================================

cd "$(dirname "$0")/.." || exit 1          # raiz do projeto (os ps*.txt)
B=bateria
# Em outra maquina, MAQUINA=nome separa o CSV e as saidas, para nao misturar
# com as medidas desta (resultados_nome.csv, saidas_nome/). Nome curto: o
# MSCD corta caminho de arquivo em 100 caracteres (mscdrun.cpp:70).
if [ -n "$MAQUINA" ]; then
  CSV=$B/resultados_$MAQUINA.csv; SAIDAS=$B/saidas_$MAQUINA
else
  CSV=$B/resultados.csv; SAIDAS=$B/saidas
fi

VERSOES=${VERSOES:-"original gpu"}
# raio:profundidade
TAMANHOS=${TAMANHOS:-"8:12 10:15 12:18 14:21 16:24 17:26"}
# TAMANHOS=centenas: perto de 100, 200, ... 1100 atomos (contados em 09/10/2026)
#   96 210 304 394 508 592 686 785 890 1016 1112
[ "$TAMANHOS" = centenas ] && TAMANHOS="7:10.5 9.25:13.875 10.75:16.125 \
11.75:17.625 13:19.5 13.75:20.625 14.25:21.375 15.25:22.875 15.75:23.625 \
16.5:24.75 17:26"
NPS=${NPS:-"1 2 3 4 5 6 7 8 9 10 11 12"}
TEMPO_MAX=${TEMPO_MAX:-36000}              # segundos por rodada (10 h)
FOLGA_MEM=${FOLGA_MEM:-0.85}               # fracao da RAM livre que pode usar

# ---------------------------------------------------------------------
resumo() {
  [ -f "$CSV" ] || { echo "ainda nao ha resultados"; return; }
  column -s, -t "$CSV"
}
[ "$1" = "resumo" ] && { resumo; exit 0; }

# ---------------------------------------------------------------------
[ -f Cov0.txt ] || { echo "falta o Cov0.txt na raiz do projeto"; exit 1; }
for v in $VERSOES; do
  [ -x $B/bin/randmscd_$v ] || { echo "falta $B/bin/randmscd_$v. Os executaveis nao vem pelo git, rode antes ./bateria/instalar_ubuntu.sh ate o fim"; exit 1; }
done
if ps -eo comm= | grep -q "^randmscd"; then
  echo "ja' ha' um MSCD rodando nesta maquina; a medida sairia contaminada."
  echo "pare o outro primeiro (pgrep -af randmscd)."
  exit 1
fi
rm -f $B/PAUSAR
mkdir -p $B/entradas $SAIDAS
[ -f "$CSV" ] || echo "versao,raio,prof,natoms,np,segundos,status,rfactor_a,rfactor_b,pico_mem_MB,carga_inicial,data" > "$CSV"

# Ctrl+C: mata a rodada atual e sai; ela nao entra no CSV e sera' refeita.
FILHO=""
parar() {
  echo; echo ">> interrompido; a rodada atual sera' refeita na proxima vez."
  [ -n "$FILHO" ] && kill -- -"$FILHO" 2>/dev/null
  pkill -f "[b]ateria/bin/randmscd_" 2>/dev/null
  exit 130
}
trap parar INT TERM

ja_feito() {   # versao raio prof np
  grep -q "^$1,$2,$3,[0-9]*,$4," "$CSV"
}
pulado_tamanho() {   # versao natoms: este tamanho ou um menor ja' estourou o tempo?
  # O que trava no original e' o preparo, que roda so' no processo 0: se uma
  # rodada nao termina em TEMPO_MAX, outro -np do mesmo tamanho tambem nao.
  awk -F, -v v="$1" -v n="$2" '$1==v && $7=="TEMPO_ESGOTADO" && $4<=n {f=1} END{exit !f}' "$CSV"
}

entrada() {    # raio prof -> gera a entrada e devolve o caminho
  local r=$1 p=$2 f=$B/entradas/r${1}p${2}.in
  sed -E -e "s|^pe([[:space:]]+)[^[:space:]]+|pe\1$SAIDAS/SAIDA|" \
         -e "s|^[0-9.]+([[:space:]]+)[0-9.]+([[:space:]]+)([0-9.]+[[:space:]]+radius)|$r\1$p\2\3|" \
         Cov0.txt > "$f"
  echo "$f"
}

natoms_de() {  # entrada -> numero de atomos (le o cabecalho e corta)
  timeout 20 $B/bin/randmscd_original "$1" 2>/dev/null \
    | awk '/natoms emiters/{print $1; exit}'
}

mem_estimada_MB() {   # versao natoms np
  # medido: o original guarda 4 tabelas natoms^3 (16 bytes por trio) e cada
  # processo recebe a sua copia; a GPU com np=1 guarda ~9 bytes por trio.
  awk -v v="$1" -v n="$2" -v np="$3" 'BEGIN{
    t=n*n*n;
    if (v=="gpu" && np==1) m=10*t;
    else if (np==1) m=16*t;
    else m=np*1.6*16*t;
    printf "%d", m/1048576+300 }'
}

rodar() {      # versao raio prof natoms np entrada
  local v=$1 r=$2 p=$3 n=$4 np=$5 in=$6
  local tag=${v}_r${r}p${p}_np${np}
  local run=$B/entradas/$tag.in out=$SAIDAS/$tag.out
  sed "s|$SAIDAS/SAIDA|$out|" "$in" > "$run"
  local extra="" envs=""
  [ "$np" -gt 6 ] && extra="--use-hwthread-cpus"
  if [ "$v" = "gpu" ]; then
    envs="MSCD_GPU=1"
    extra="$extra --bind-to none"
    # o pathcut na placa precisa de 8 bytes por trio de VRAM
    local vram=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1)
    local prec=$(awk -v n="$n" 'BEGIN{printf "%d", 8*n*n*n/1048576+600}')
    [ -n "$vram" ] && [ "$prec" -gt "$vram" ] && envs="$envs MSCD_PATHCUTCPU=1"
  fi
  local carga=$(cut -d' ' -f1 /proc/loadavg)
  local t0=$(date +%s.%N)
  setsid env $envs /usr/bin/time -f "%M" -o $SAIDAS/$tag.mem \
    timeout --kill-after=30 $TEMPO_MAX \
    mpirun $extra -np $np $B/bin/randmscd_$v "$run" \
    > $SAIDAS/$tag.log 2> $SAIDAS/$tag.err &
  FILHO=$!
  wait $FILHO; local rc=$?
  FILHO=""
  local t1=$(date +%s.%N)
  local seg=$(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.1f", b-a}')
  local status=OK
  [ $rc -eq 124 ] || [ $rc -eq 137 ] && status=TEMPO_ESGOTADO
  [ $status = OK ] && [ $rc -ne 0 ] && status="ERRO_$rc"
  local rf=$(grep -h "r-factora =" "$out" 2>/dev/null | awk '{print $3","$6}')
  [ -z "$rf" ] && rf=","
  [ $status = OK ] && [ "$rf" = "1.0000,-1.0000" ] && status=CURVA_ZERADA
  local mem=$(tail -1 $SAIDAS/$tag.mem 2>/dev/null | awk '{printf "%d", $1/1024}')
  echo "$v,$r,$p,$n,$np,$seg,$status,$rf,$mem,$carga,$(date +%F_%T)" >> "$CSV"
  printf "   %-8s %4s atomos  np=%-2s  %8s s  %s\n" "$v" "$n" "$np" "$seg" "$status"
  rm -f "$run"
}

# ---------------------------------------------------------------------
echo "bateria: versoes [$VERSOES], tamanhos [$TAMANHOS], np [$NPS], tempo max ${TEMPO_MAX}s"
echo "para pausar depois desta rodada: touch $B/PAUSAR    (Ctrl+C para na hora)"
for v in $VERSOES; do
  for rp in $TAMANHOS; do
    r=${rp%%:*}; p=${rp##*:}
    in=$(entrada $r $p)
    n=$(natoms_de "$in")
    [ -z "$n" ] && { echo "r=$r p=$p: nao consegui ler o numero de atomos (passa de 1250?)"; continue; }
    echo "== $v  raio $r  prof $p  ($n atomos)"
    for np in $NPS; do
      ja_feito $v $r $p $np && continue
      if [ -f $B/PAUSAR ]; then
        rm -f $B/PAUSAR
        echo ">> pausado. Rode o script de novo para continuar daqui."
        exit 0
      fi
      if pulado_tamanho $v $n; then
        echo "$v,$r,$p,$n,$np,,PULADO_JA_ESGOTOU_NESTE_TAMANHO_OU_MENOR,,,,,$(date +%F_%T)" >> "$CSV"
        continue
      fi
      livre=$(awk '/MemAvailable/{printf "%d", $2/1024}' /proc/meminfo)
      prec=$(mem_estimada_MB $v $n $np)
      if [ $(awk -v a=$prec -v b=$livre -v f=$FOLGA_MEM 'BEGIN{print (a>b*f)}') = 1 ]; then
        echo "$v,$r,$p,$n,$np,,PULADO_MEMORIA_${prec}MB_de_${livre}MB,,,,,$(date +%F_%T)" >> "$CSV"
        printf "   %-8s %4s atomos  np=%-2s  pulado: precisaria ~%s MB, livres %s MB\n" "$v" "$n" "$np" "$prec" "$livre"
        continue
      fi
      rodar $v $r $p $n $np "$in"
    done
  done
done
echo ">> bateria completa."
resumo
