#!/bin/bash
# =====================================================================
# Prepara um Ubuntu limpo (22.04, 24.04 ou 26.04) com placa NVIDIA para rodar a
# bateria da GPU.
#
#   ./bateria/instalar_ubuntu.sh
#
# 1. instala compilador, Open MPI e o pacote time
# 2. instala o driver da NVIDIA, se faltar. Ai' pede para reiniciar e rodar
#    este script de novo, que continua de onde parou
# 3. instala o CUDA da NVIDIA, na versao mais nova que o driver aceita
# 4. compila os dois executaveis em bateria/build e copia para bateria/bin
#      randmscd_gpu       o codigo do ultimo commit (HEAD)
#      randmscd_original  o commit d4c8408 com o limite de 300 para 1250
#                         atomos. A bateria so' usa para contar os atomos
# 5. roda 135 atomos na GPU e compara com a saida desta bateria medida no
#    PC de origem (bateria/saidas/gpu_r8p12_np1.out)
#
# Os dois executaveis foram conferidos no PC de origem em 09/10/2026. O
# original reconstruido assim sai com o mesmo md5 do da bateria, e o de GPU
# da' a mesma saida byte a byte (so' a data e a hora mudam).
# =====================================================================
set -eo pipefail
trap 'echo; echo "ERRO: o script parou na linha $LINENO, no comando: $BASH_COMMAND"' ERR
cd "$(dirname "$0")/.." || exit 1
B=bateria
MAQ=${MAQUINA:-$(hostname -s | tr -cd 'A-Za-z0-9_-' | cut -c1-12)}
msg() { echo; echo ">> $*"; }
erro() { echo; echo "ERRO: $*"; exit 1; }

# ---------------------------------------------------------------------
. /etc/os-release
[ "$ID" = ubuntu ] || erro "feito para Ubuntu, aqui e' $ID"
case "$VERSION_ID" in 22.04|24.04|26.04) ;;
  *) echo "aviso: pensado para 22.04, 24.04 e 26.04, aqui e' $VERSION_ID" ;; esac
git rev-parse -q --verify d4c8408^{commit} >/dev/null \
  || erro "o clone nao tem o commit d4c8408. Clone sem --depth."
[ -f $B/saidas/gpu_r8p12_np1.out ] \
  || erro "falta $B/saidas/gpu_r8p12_np1.out, a referencia do teste final"

# ---------------------------------------------------------------------
msg "1/5 pacotes do sistema"
sudo apt-get update
sudo apt-get install -y build-essential git openmpi-bin libopenmpi-dev \
  time pciutils wget ca-certificates ubuntu-drivers-common

# ---------------------------------------------------------------------
msg "2/5 driver da NVIDIA"
lspci | grep -qi nvidia || erro "nenhuma placa NVIDIA no lspci"
if ! nvidia-smi >/dev/null 2>&1; then
  sudo ubuntu-drivers install || sudo ubuntu-drivers autoinstall
  echo
  echo "Driver instalado. Reinicie (sudo reboot) e rode este script de novo."
  echo "Se o Secure Boot estiver ligado, na reinicializacao aparece uma tela"
  echo "azul (MOK). Escolha 'Enroll MOK' e digite a senha que o apt pediu."
  exit 0
fi
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

# ---------------------------------------------------------------------
msg "3/5 CUDA"
CUDA_MAX=$(nvidia-smi | grep -oP 'CUDA (UMD )?Version: \K[0-9]+\.[0-9]+' || true)
[ -n "$CUDA_MAX" ] || erro "nao consegui ler a versao do CUDA no nvidia-smi"
echo "o driver aceita ate' CUDA $CUDA_MAX"
achar_cuda() {   # o do repositorio da NVIDIA primeiro, o do Ubuntu (/usr) por ultimo
  local d
  for d in /usr/local/cuda $(ls -d /usr/local/cuda-* 2>/dev/null | sort -V -r) /usr; do
    [ -x "$d/bin/nvcc" ] && { echo "$d"; return; }
  done
  return 0   # sem isto, nao achar o nvcc derruba o script pelo set -e
}
toolkit_nvidia() {   # repositorio da NVIDIA de uma versao do Ubuntu (ex. 2604)
  local url=https://developer.download.nvidia.com/compute/cuda/repos/ubuntu$1/x86_64
  wget -q -O /tmp/cuda-keyring.deb $url/cuda-keyring_1.1-1_all.deb || return 1
  sudo dpkg -i /tmp/cuda-keyring.deb && sudo apt-get update || return 1
  # o toolkit mais novo que o driver aceita
  local pkg=$(apt-cache pkgnames cuda-toolkit- | grep -E '^cuda-toolkit-[0-9]+-[0-9]+$' \
    | awk -F- -v max="$CUDA_MAX" '{ split(max,m,".");
        if ($3<m[1] || ($3==m[1] && $4<=m[2])) print }' | sort -V | tail -1)
  [ -n "$pkg" ] || return 1
  echo "instalando $pkg do repositorio ubuntu$1 da NVIDIA"
  sudo apt-get install -y "$pkg"
}
CUDA_PATH=$(achar_cuda)
if [ -z "$CUDA_PATH" ]; then
  # 1o o repositorio desta versao. Se a NVIDIA ainda nao tiver um para ela,
  # o do 24.04. Por ultimo o pacote do proprio Ubuntu.
  toolkit_nvidia ${VERSION_ID/./} || toolkit_nvidia 2404 \
    || sudo apt-get install -y nvidia-cuda-toolkit \
    || erro "nao consegui instalar o CUDA"
  CUDA_PATH=$(achar_cuda)
  [ -n "$CUDA_PATH" ] || erro "instalei o CUDA mas nao achei o nvcc"
fi
NVCC=$CUDA_PATH/bin/nvcc
CUDA_LIB=$CUDA_PATH/lib64
[ -d "$CUDA_LIB" ] || CUDA_LIB=/usr/lib/x86_64-linux-gnu
$NVCC --version | tail -2

# arquitetura da placa (8.9 -> sm_89, a RTX 4050 e a 4060 sao 8.9)
CC=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '. ')
ARCH=${CC:+sm_$CC}; ARCH=${ARCH:-native}

# O nvcc recusa g++ mais novo do que conhece. Se recusar o padrao, tenta os
# anteriores, do mais novo para o mais velho.
mkdir -p $B/build
printf '__global__ void k(){}\nint main(){k<<<1,1>>>();return 0;}\n' > $B/build/t.cu
nvcc_ok() { $NVCC -arch=$ARCH $B/build/t.cu -o $B/build/t 2>$B/build/t.err; }
if ! nvcc_ok; then
  for v in 15 14 13 12; do
    apt-cache show g++-$v >/dev/null 2>&1 || continue
    echo "o nvcc recusou o g++ padrao, tentando com g++-$v"
    sudo apt-get install -y g++-$v
    export NVCC_CCBIN=g++-$v
    nvcc_ok && break
    unset NVCC_CCBIN
  done
  [ -n "$NVCC_CCBIN" ] || { cat $B/build/t.err; erro "o nvcc nao aceitou nenhum g++"; }
fi
rm -f $B/build/t $B/build/t.cu $B/build/t.err

# ---------------------------------------------------------------------
msg "4/5 compilando em $B/build"
rm -rf $B/build/original $B/build/gpu
mkdir -p $B/build/original $B/build/gpu $B/bin

git archive d4c8408 | tar -x -C $B/build/original
( cd $B/build/original
  sed -i -e 's/300\*12/1250*12/g' -e 's/n>=300)/n>=1250)/' \
         -e 's/natoms>300)/natoms>1250)/' mscdrun.cpp mscdruna.cpp mscdrun.h
  n=$(cat mscdrun.cpp mscdruna.cpp mscdrun.h | grep -c 1250)
  [ "$n" = 12 ] || { echo "troca do limite mudou $n linhas, esperava 12"; exit 1; }
  make randmscd_parallel CPPFLAGS="-O3 -std=c++98 -w -fpermissive" \
    > build.log 2>&1 || { tail -20 build.log; exit 1; }
) || erro "falhou a compilacao do original ($B/build/original/build.log)"
cp $B/build/original/randmscd_parallel $B/bin/randmscd_original

# -fmad=false nao e' opcional (ver o makefile e o PLANO_CUDA.md)
git archive HEAD | tar -x -C $B/build/gpu
( cd $B/build/gpu
  make randmscd_gpu CPPFLAGS="-O3 -std=c++98 -w -fpermissive -fopenmp -DMSCDGPU" \
    NVCC="$NVCC" NVFLAGS="-O3 -arch=$ARCH -lineinfo -fmad=false" \
    CUDALIBS="-L$CUDA_LIB -lcudart -Wl,-rpath,$CUDA_LIB" \
    > build.log 2>&1 || { tail -20 build.log; exit 1; }
) || erro "falhou a compilacao da GPU ($B/build/gpu/build.log)"
cp $B/build/gpu/randmscd_gpu $B/bin/randmscd_gpu
( cd $B/bin && md5sum randmscd_gpu randmscd_original > md5_$MAQ.txt )
echo "compilado para $ARCH com $($NVCC --version | grep -o 'release [0-9.]*')"

# ---------------------------------------------------------------------
msg "5/5 teste com 135 atomos"
# mesma troca que o bateria.sh faz. Caminhos curtos, o MSCD corta em 100
# caracteres.
sed -E -e "s|^pe([[:space:]]+)[^[:space:]]+|pe\1$B/build/teste.out|" \
       -e "s|^[0-9.]+([[:space:]]+)[0-9.]+([[:space:]]+)([0-9.]+[[:space:]]+radius)|8\112\2\3|" \
       Cov0.txt > $B/build/teste.in
n=$(timeout 20 $B/bin/randmscd_original $B/build/teste.in 2>/dev/null \
    | awk '/natoms emiters/{print $1; exit}' || true)
[ "$n" = 135 ] || erro "o original contou '$n' atomos, esperava 135"
echo "original conta 135 atomos, ok"

rm -f $B/build/teste.out
MSCD_GPU=1 mpirun --bind-to none -np 1 $B/bin/randmscd_gpu $B/build/teste.in \
  > $B/build/teste.log 2> $B/build/teste.err \
  || { tail $B/build/teste.err; erro "a GPU terminou com erro"; }
rf=$(grep -h "r-factora =" $B/build/teste.out | awk '{print $3, $6}' || true)
echo "r-factor $rf (no PC de origem 0.6836 0.8648)"
[ "$rf" != "1.0000 -1.0000" ] \
  || { cat $B/build/teste.err; erro "curva zerada, a GPU falhou no setup"; }

# chical e' a coluna 4 das linhas de dados, criterio do baseline/regressao-gpu.sh
chi() { awk 'NF==5 && $1+0==$1 && $2 ~ /[eE]/ {print $4}' "$1"; }
paste <(chi $B/saidas/gpu_r8p12_np1.out) <(chi $B/build/teste.out) | awk '
  { d=$2-$1; if (d<0) d=-d; if (d>m) m=d; n++ }
  END { printf "%d pontos, max|dchi| = %.2e contra o PC de origem (criterio 1e-4)\n", n, m;
        exit !(n>100 && m<=1e-4) }' \
  || erro "a curva difere da do PC de origem acima de 1e-4"

# ---------------------------------------------------------------------
nt=$(nproc)
msg "pronto. Para rodar a bateria so' da GPU"
echo
echo "    MAQUINA=$MAQ VERSOES=gpu TAMANHOS=centenas ./bateria/bateria.sh"
echo
echo "Os resultados vao para $B/resultados_$MAQ.csv e $B/saidas_$MAQ/,"
echo "separados dos do PC de origem."
[ "$nt" -lt 12 ] && echo "aviso: esta maquina tem $nt threads, e -np acima disso falha. Use NPS=\"\$(seq -s' ' 1 $nt)\"."
echo "Notebook: deixe na tomada e no modo de desempenho maximo durante a bateria."
