# Plano do port de CUDA

Escrito em 05/08/2026 sobre o V5; **Fases 0 e 1 executadas no mesmo dia.** Tudo
que está aqui vem de medição, não de leitura do código — e a leitura já errou
**quatro** vezes neste projeto: a fase C do `symtrivert`, o V3, o `allrotation`
que eu suspeitei ser gargalo e é 1,4%, e o desenho `bloco ↦ par` que **este
próprio arquivo** exigia e a Fase 1 derrubou. A quarta é a mais instrutiva: o
plano estava certo sobre *onde* o tempo estava e errado sobre *o que o kernel
precisava*, e só a inspeção das funções chamadas mostrou isso.

**Objetivo final:** a figura **V0 × V5 × CUDA** do paper. Ver a última seção.

## A meta, posta pelo usuário em 05/08/2026

**Chegar a ~10 s de tempo total.** A trajetória até aqui, em `np` de produção:

| | tempo | como |
|---|---:|---|
| V0 | 70–80 s | original |
| V4 / V5 | ~40 s | `symtrivert`, corte do `tevencut`, OpenMP no `pathcut` |
| GPU Fase 1 | 60–64 s (`np=1`) | **ainda perde para o V5 em `np=12`** |
| **meta** | **~10 s** | Fases 2 e 3, talvez a 4 |

**É plausível, e o piso é conhecido.** Com as Fases 2 e 3 o laço inteiro sai do
host e sobra o **preparo serial: 5,4 s**. A Fase 4 (`pathcut` na placa) atacaria
1,71 s desses. Então o chão realista é ~4–6 s de preparo mais o tempo dos
kernels — 10 s é apertado mas não é fantasia.

**O que decide se dá ou não é a Fase 3** (`summation`, 34% do laço). Ela tem uma
dependência **serial em `m`** (7 passos, `bsum(m)` alimenta `asum(m+1)`) que não
some; o paralelismo tem de vir dos pares dentro de cada `m`. Se a Fase 3 não
render, a meta vira ~15–20 s. **Meça a Fase 3 antes de prometer os 10 s a
alguém.**

---

## COMO CONTINUAR

**Fases 0 e 1 estão FEITAS e validadas (05/08/2026). `np=1` foi de 130,58 s para
64,15 s — 2,04×, com `max|Δχ| = 1,0×10⁻⁵`.**

### Leia isto antes de comemorar o 2,04×

**O V5 com todos os núcleos ainda é mais rápido que o build de GPU.** Na janela
de 05/08 tarde: V5 `np=12` ≈ **42 s** contra GPU `np=1` ≈ **60–64 s**.

Não há contradição com o 2,04×: aquele número é GPU × V5 **ambos em `np=1`**, que
é o controle certo para isolar o efeito do kernel — e é **errado** citá-lo como
"a GPU ganhou do V5". Enquanto 43% do laço estiver serial no host, um núcleo
sozinho não alcança doze ranks dividindo o laço inteiro.

**E não confunda com o "Analyzing/Reanalyzing".** Aquilo é o `symtrivert`, do
preparo serial, alvo do V1/V2. A GPU **não encosta nele**. O que subiu foi o
`alldblevent`, de dentro do laço dos 779 pontos.

O dado que importa para decidir o que fazer: dos 66,4 s removidos em `np=1`, o
kernel **mais** o tráfego de PCIe custaram **~5 s** para os 779 pontos. **A placa
está praticamente ociosa.** Com as Fases 2 e 3 sobra o preparo serial (~5,4 s) e
o alvo plausível cai para a casa de **10–15 s** — aí sim passa o `np=12` com
folga. É por isso que parar na Fase 1 não é uma opção.

**Ideia barata e não testada:** rodar o build de GPU com `np=4` ou mais. Cada
rank chamaria o kernel para os seus pontos; as tabelas são ~330 MB por rank e
caberiam 4–8 ranks nos 8 GB. Não se sabe se a contenção na placa mata o ganho.
É **uma rodada de medição**, não um port — e pode valer mais que a Fase 2 se der
certo. Meça antes de assumir qualquer coisa.

### A próxima ação, concreta

**Fase 2: `allevendetec` (`mscdrunc.cpp:893`) para a placa.** Nesta ordem:

1. Ler a função inteira. Ela é `natoms²×radim` = 915 015 elementos, cópia com
   escala complexa, zero dependência entre elementos.
2. **A razão de fazer não são os 8,8%.** É que hoje `mscdgpu_alldblevent`
   termina com um `cudaMemcpy` de volta de 4,9 MB por ponto — **3,8 GB de PCIe
   na corrida** — só porque o consumidor seguinte está no host. Com a Fase 2 na
   placa, o `devenelem` nasce e morre lá.
3. Mudar a interface: `mscdgpu_alldblevent` para de devolver `devenelem`.
   **Manter um `mscdgpu_get_devenelem()`** só para o `MSCD_GPU=validate`, senão
   se perde a única ferramenta que achou os bugs da Fase 1.
4. Validar com `MSCD_GPU=validate` **antes** de encostar em `summation`, e
   fechar com `./baseline/regressao-gpu.sh 1`.

**Antes de escrever qualquer linha, leia a seção "Armadilhas específicas do
port".** Os quatro primeiros itens custaram a depuração inteira da Fase 1, e o
sintoma deles é sempre "quase certo, com uma fração minúscula muito errada".

**Estado do código:** V5 + Fase 1 de CUDA. Dois executáveis:

```bash
# CPU, produção (inalterado)
make randmscd_parallel CPPFLAGS="-O3 -std=c++98 -w -fpermissive -fopenmp"

# GPU (exige rm -f *.o: os .o precisam de -DMSCDGPU)
rm -f *.o && make randmscd_gpu \
  CPPFLAGS="-O3 -std=c++98 -w -fpermissive -fopenmp -DMSCDGPU"

MSCD_GPU=1        mpirun --use-hwthread-cpus --bind-to none -np 1 randmscd_gpu Cov0.txt
MSCD_GPU=validate mpirun ... randmscd_gpu Cov0.txt   # roda CPU e GPU, compara
./baseline/regressao-gpu.sh 1                        # o teste, com o critério
```

Sem `MSCD_GPU` o `randmscd_gpu` roda 100% no host — o mesmo caminho do
`randmscd_parallel`.

**O que já foi decidido e não deve ser refeito:**

- **Fase 0** — critério `max |Δχ| ≤ 1e-4`, medido, não estimado. Ver a seção.
- **O atalho do `beta` inteiro é REAL** — 30 726 810 de 30 726 810. Ver a seção.
- **A otimização de CPU que saía de graça daí (V6) foi medida e REPROVADA.**
  Byte a byte idêntica, +1,5% em `np=1` e −3,5% em `np=12`. Está revertida e o
  motivo está em `OTIMIZACAO.md`, seção V6. **Não retente sem ler.**
- **`bloco ↦ par` foi descartado a favor de `thread ↦ par`** — o motivo que
  justificava blocos evaporou. Ver "O desenho mudou".
- **`-fmad=false` é obrigatório** no `nvcc`. Ver "Armadilhas do port".

---

## Onde estamos (05/08/2026, medido)

| | valor |
|---|---:|
| **melhor tempo `np=1` com a Fase 1 na placa** | **59,84 s** |
| melhor tempo `np=1` só CPU (V5) | 122,02 s |
| preparo serial | 5,4 s (4,4%) |
| laço dos 779 pontos | ~115 s (95,6%) |
| teto de Amdahl se só o laço for para a GPU | 22,3× |
| teto se só a Fase 1 for para a GPU | 2,20× |
| **medido com a Fase 1** | **2,04×** (91% do teto dela) |

**Perfil do laço** (`np=1`, acumuladores fechando em 99,7%):

| bloco | % do laço | onde |
|---|---:|---|
| **`alldblevent`** | **57,0%** | `mscdrunc.cpp:753` |
| `summation` — laço de `m` | 22,1% | `mscdrund.cpp:108` |
| `summation` — bloco final | 9,0% | `mscdrund.cpp:183` |
| `allevendetec` | 8,8% | `mscdrunc.cpp:879` |
| `summation` — init do `asum` | 3,1% | `mscdrund.cpp:97` |

**Hardware:** RTX 4060 Laptop, 8188 MB, driver 596.49 (CUDA 13.2), `nvcc` 12.0,
WSL2 com `/dev/dxg`. Compute 8.9, ~3072 cores FP32, ~15 TFLOPS FP32.
**FP64 é 1/64 nesta placa** — o código é `float` inteiro e tem de continuar
sendo. Não promova nada para `double` "para melhorar a precisão".

---

## Fase 0 — o critério de validação. **FEITO.**

**A validação byte a byte não sobrevive ao port**, porque a redução soma em outra
ordem e em `float` isso mexe no último bit. O critério novo tinha de ser fixado
*antes* de existir kernel, senão não há como distinguir "reordenou a soma" de
"errou a física".

**O piso de ruído foi medido**, não estimado. Inverter a ordem da soma sobre `al`
em `evenelem` (`mscdrunc.cpp:34`) é fisicamente idêntico e muda só o
arredondamento — exatamente o que uma redução de GPU faz:

| | valor |
|---|---:|
| max \|Δχ\| | **1,0×10⁻⁵** |
| rms \|Δχ\| | 3,9×10⁻⁷ |
| max \|ΔI/I\| | 1,7×10⁻⁵ |

E 1,0×10⁻⁵ é **um dígito na última casa impressa** para os χ maiores (o arquivo
tem 5 algarismos significativos): a reordenação mexe, no máximo, no último bit
que a saída consegue mostrar.

Para calibrar a escala: os 779 pontos têm χ de **−0,49043 a +0,42295**
(amplitude 0,91338), rms 0,14195, e o menor |χ| não nulo é 1,6×10⁻⁴.

### O critério, decidido

    max |Δχ| <= 1e-4   sobre os 779 pontos      <- criterio primario
    rms |Δχ| <= 1e-5

1×10⁻⁴ é **10× o piso de ruído medido** e **0,011% da amplitude** da curva.
Qualquer coisa acima disso não é reordenação de soma: é bug.

**Três coisas que este critério NÃO é:**

- **Não é tolerância relativa por ponto.** O menor |χ| é 1,6×10⁻⁴; exigir 1% dele
  seria exigir 1,6×10⁻⁶ absoluto, abaixo do piso de ruído, e o teste reprovaria
  código correto. Tolerância relativa por ponto está **descartada por medição**.
- **Não é o R-factor.** Ele mal se moveu quando 783 das 787 linhas mudaram na
  correção do k. Serve de sanidade, nunca de critério.
- **Não são os `factors`.** Mesma razão — são fatores de escala do ajuste.

### O que fazer com o `regressao.sh`

Ele continua valendo **sem alteração** para qualquer mudança de CPU (V6 em
diante): byte a byte é o critério mais forte e não custa nada manter. Para o
caminho de GPU, escrever um `baseline/regressao-gpu.sh` separado com o critério
acima. **Não relaxe o `regressao.sh` existente** — perder o teste byte a byte no
CPU seria pagar duas vezes.

---

## Fase 1 — `alldblevent` → `evenelem`. **FEITA e validada.**

`mscdgpu.cu` + `mscdgpu.h`, ligados por `-DMSCDGPU`. `mscdrunc.cpp:768`
(`Mscdrun::gpudblevent`) é a cola; `mscdrund.cpp:298` despacha.

**Resultado**, `./baseline/regressao-gpu.sh 1`:

```
max|dchi| = 1.000e-05   (criterio 1e-4)   <- e' o piso de ruido medido na Fase 0
rms|dchi| = 6.432e-07   (criterio 1e-5)
pontos acima de 1e-4: 0
```

O erro ficou **exatamente no piso de ruído de ponto flutuante do programa**, 10×
abaixo do critério. `factors` inalterados.

**Tempo**, campanha pareada e alternada V5 × GPU, `np=1`, mesma janela
(05/08/2026 tarde, `$CLAUDE_JOB_DIR/tmp/gpuab.csv` — copiar para `baseline/` se
for usar no paper):

| rep | V5 | GPU |
|---|---:|---:|
| 1 | 130,58 s | 64,15 s |
| 2 | 133,45 s | 69,13 s |
| 3 | 145,08 s | 69,58 s |
| **mínimo** | **130,58 s** | **64,15 s** |

**2,04× pelos mínimos**, 2,02× pelas médias. Melhor rodada isolada depois:
**59,84 s**.

**Isso é ~91% do teto teórico da Fase 1**: se só o `alldblevent` (57,0% do laço,
que é 95,6% do total) for para a placa, o máximo é
`1/(1−0,570×0,956) = 2,20×`. Para passar disso é preciso a Fase 2 e a Fase 3 —
não há mais nada a extrair deste kernel.

**Aviso de leitura:** a terceira repetição do V5 (145,08 s) destoa; a máquina
derivou durante a campanha. Por isso a tabela é pelo **mínimo**, como o resto do
projeto. Para o paper isto tem de ser remedido junto com V0 e V5 numa campanha
única — ver a última seção.

### O atalho do `beta` inteiro — confirmado por contagem

Medido em 05/08/2026 com `-DBETATEST` (o contador continua em
`rotamat.cpp:336`, desligado):

| origem | reconstruções de `rotmatb` | com `(xa-i) != 0` |
|---|---:|---:|
| **dentro de `alldblevent`** | 30 726 810 | **0** |
| fora | 11 286 832 | 11 180 536 |

As 346 exceções são `beta=180` exato (fator 1,0, ainda cópia de linha). E
**`beta` nunca é negativo**: o ramo de troca de sinal do `makerotation` é código
morto aqui. O kernel indexa `rotmata` direto e **nunca monta `rotmatb`**.

Detalhe: em `beta=180` a CPU satura `i` em 359 e calcula `a+1,0*(b−a)`, que não é
bit a bit `b`; o kernel lê a linha 360. São 346 pares em 30,7 milhões e o valor
do kernel é o certo.

**A otimização de CPU que parecia sair de graça daí (V6) foi medida e reprovada
— ver `OTIMIZACAO.md`.** O achado vale para a GPU, não para o host.

### O desenho mudou: **thread ↦ par**, não bloco ↦ par

A versão anterior deste plano exigia bloco ↦ par porque cada thread teria de
reconstruir ~2,5 KB de tabela. **Esse motivo evaporou**, por dois achados:

- **`rotmatb` (500 floats) não precisa existir** — é o atalho do `beta`.
- **`hankarg` (100 complexos) também não.** `fhankelfaca` é uma interpolação de
  **dois pontos** de `hankmat`; o cache existe na CPU só para amortizar as 15
  chamadas por par. No device sai inline: 2 leituras e um lerp por `(al,am)`,
  contra 100 entradas construídas. `hankmat` tem 205 KB e vive em L2.
- **`fsinexpa` é keyed em `akin`, constante** ⇒ tabela fixa, sobe uma vez.

Sobra um kernel só de registrador, sem *shared memory* nenhuma. 128 threads por
bloco.

### Três coisas que a leitura do código não mostrava

Custaram medição; estão no `mscdgpu.cu` como comentário colado na linha.

1. **`vkb` é sempre `0.0f`** (`mscdrunc.cpp:814`) e **`mb=nb=0` nas 15
   chamadas** ⇒ `kb=0`. A tabela do `hankb` é **constante na corrida inteira**.
   Ela é *fotografada* do objeto da CPU em vez de reproduzida — o cache de
   `fhankelfaca` é keyed no argumento e outras funções do laço mexem nele, e
   copiar 100 complexos é mais barato que reproduzir o autômato.
2. **As 15 saídas vêm de 9 `evenelem` distintos**, com `(ma,na)` em
   `{(0,0),(1,0),(0,1),(2,0),(1,1),(3,0),(0,2),(2,1),(4,0)}`. Compartilham o
   laço de `al`: 9 acumuladores complexos numa passada só.
3. **`getkelem(ma,0) == getkharm(ma,0) == ma·(ma+1)`**, sempre ≥ 0 ⇒ o ramo de
   sinal de `rotharma` (`kelem<0`) nunca dispara.

### A geometria: o `onerotation` refazia 779× o que é fixo

Em `alldblevent` o terceiro ponto é `patom[ib]+xdetec`, então o vetor
`(patomc−patomb)` **é o `xdetec`** — igual para todos os 40 562 pares num dado
ponto. E a perna `(patomb−patoma)` só depende do par, **não do detector**.

Então `cosb,sinb,phib` são 3 escalares por ponto (host) e `cosa,sina,phia` são um
precomputo de uma vez só (`k_pairgeo`). O `alpha` não é consumido por
`alldblevent`, e seus dois `atan2` saem junto. Sobram 4 transcendentais por par
por ponto.

**Isto também é uma otimização de CPU não colhida** — ver "O que sobrou".

### Como foi validado

`MSCD_GPU=validate` roda o `alldblevent` da CPU e o kernel lado a lado e compara
`devenelem` **elemento a elemento**, sem tocar no resultado. Foi assim que os
dois bugs de promoção de tipo apareceram; sem isso teriam virado "a curva está
meio diferente" no fim do laço.

---

## Fase 2 — `allevendetec` (8,8%). **PRÓXIMA.**

Cópia com escala complexa, `natoms²×radim` = 915 015 elementos, zero dependência.
Trivial. **A razão de fazer não são os 8,8%** — é que hoje o `mscdgpu_alldblevent`
termina com um `cudaMemcpy` de volta de 4,9 MB por ponto, **3,8 GB de PCIe na
corrida**, só porque o consumidor seguinte está no host. Subindo a Fase 2, o
`devenelem` nasce e morre na placa.

Isso muda a interface: `mscdgpu_alldblevent` deixa de devolver `devenelem` e
passa a deixá-lo no device. O modo `MSCD_GPU=validate` precisa continuar
funcionando — mantenha um `mscdgpu_get_devenelem()` só para ele.

---

## Fase 3 — `summation` (34% somando os três pedaços)

O laço de `m` (`mscdrund.cpp:108`) é **serial em `m`** — 7 passos, `msorder=8`,
`bsum(m)` alimenta `asum(m+1)`. Dentro de cada `m`, os pares `(ia,ib)` são
independentes.

Duas coisas saem de graça:

- **A cópia `bsum ← asum` some.** Ping-pong de ponteiros. É o item que ficou
  pendente do V4 (~199 GB de tráfego na corrida inteira); na GPU resolve-se por
  construção, não por otimização.
- **Compactação do `tevencut`.** Em vez de lançar `natoms²` threads e descartar
  99,87%, lançar sobre a lista dos pares sobreviventes. Como o `tevencut` não
  muda na corrida (energia constante), **essa lista é construída uma vez**.

**Não invista em aritmética aqui.** Das 90 231 570 visitas `(m,ia,ib,ic)`,
**99,87% são podadas** e sobram 1 455 trios com 14 724 MACs no total. O custo
deste bloco é tráfego de memória, não conta. **A divergência de `evedim` não é
problema — não há trabalho para divergir.**

---

## Fase 4 — `pathcut` na GPU (opcional, 1,71 s)

Só depois das fases 1-3, e talvez nunca. O V5 já o levou de 7,13 s para 1,71 s
com OpenMP, e o preparo inteiro é 5,4 s de 122 s. **A decomposição da GPU é a
mesma que o V5 já usa e validou** (thread ↦ `(ib,ic)`, `ia` sequencial por
dentro), então o risco é baixo — mas o retorno também.

---

## O que sobe para a GPU, e quando

**A energia nunca muda** (`kmin=kmax=13,63`), então `alltrievent` roda uma vez e
as tabelas de geometria são constantes na corrida inteira. **Sobem uma vez.**

| array | tamanho | nota |
|---|---:|---|
| `tevenadd`, `tevendim` (int, `natoms³`) | 60,3 MB cada | |
| `talpha`, `tgamma` (float, `natoms³`) | 60,3 MB cada | |
| `tevenelem`, `tevenpar`, `tevencut` | ~50 MB | |
| `rotmata` (361×20×25) | 722 KB | cabe em L2 |
| `hankmat` (256×5×20 cplx) | 205 KB | cabe em L2 |
| phase shift, `expix` | ~70 KB | memória constante |
| **total constante** | **~330 MB** | **4% da placa** |

Estado **por ponto**: `devenelem` 4,9 MB + `devendetec` 7,3 MB + `asum`/`bsum`
7,3 MB cada ≈ **20 MB**. Com 7,6 GB livres cabem **~350 dos 779 pontos em voo**.
Na prática nem precisa: um único ponto já expõe 40 562 pares independentes, o que
satura 24 SMs. Lotes são só para esconder latência.

---

## Armadilhas específicas do port

> As quatro primeiras foram **pagas em depuração** na Fase 1, não previstas.
> O sintoma delas é sempre o mesmo e é traiçoeiro: a física fica *quase* certa,
> com uma fração minúscula dos valores muito errada. Isso porque este programa
> tira **índices inteiros de contas em float** em três lugares — a linha de
> `rotmata` (`beta`), o `k` de `fexpix` e o `i` de `fhankelfaca`. Perto da
> fronteira, 1 ulp troca a entrada da tabela inteira. **Erro de 1 ulp aqui não
> é erro de arredondamento, é erro de índice.**

- **`-fmad=false` no `nvcc` é obrigatório.** O `nvcc` contrai `a*b+c` em FMA por
  padrão; o g++ compila para x86-64 base, que não tem FMA, e não contrai. Está
  no `makefile` com o porquê.
- **Cuidado com a ordem das promoções float→double.** O C de 1998 aqui mistura
  as duas o tempo todo, e copiar a fórmula "com o mesmo significado" não basta.
  Dois casos reais:
  - `sqrt(1.0-cosa*cosa)` (`mscdrunc.cpp:694`) — `cosa*cosa` é **float**, só a
    subtração promove. Fazendo o quadrado em double, `zc` muda no 8º dígito; como
    `beta=acos(zc)` e a derivada `1/√(1−zc²)` explode perto de `zc=±1`, virou
    erro de **1 grau** na linha de `rotmata`.
  - `akin*akin*xd*(1.0-cosbeta)` (`mscdrunc.cpp:795`) — a associação à esquerda
    faz `akin*akin*xd` ser **float**; só o último fator promove.
- **O que sobra depois de acertar tudo isso** são ~60 elementos de 608 430 por
  ponto com erro relativo de **exatamente π/1800**. Essa assinatura é o `k` do
  `fexpix` errando por **um**, numa tabela de passo 0,1°, empurrado por 1 ulp de
  diferença no `gamma` entre o `atan2` da glibc e o da CUDA. **Não se conserta**
  e não precisa: somado sobre 40 562 pares, some — o `max|Δχ|` final ficou no
  piso de ruído.
- **`Fcomplex` é `{float re, im}`** (`fcomplex.h:31`), mesmo layout de `float2` —
  mas `fcomplex.h:46` declara `friend Fcomplex polar(float,float=0)`, argumento
  padrão em `friend` que não é definição. É o que obriga `-fpermissive` no g++, e
  **barra o front-end de device**. Confirmado. Por isso `mscdgpu.h` define um
  `Gcplx` POD e **nenhum cabeçalho do programa entra no `.cu`**.
- **`lnum` difere entre as espécies** (`psAg111.txt` e `psl9.txt`), mas o
  `phasec` é sempre alocado com **61 entradas fixas** (`phase.cpp:28`). O passo
  do array no device tem de ser o da alocação, não o do arquivo.
- **`tevenpar` é AoS de passo 10 floats** e o laço lê 2 campos. Separar em arrays
  próprias é obrigatório.
- **As quatro arrays `natoms³` já são coalescidas** se thread ↦ `ic`: são
  indexadas pelo mesmo `id` com `ic` variando mais rápido. Não reorganize.
- **Com uma GPU, `np>1` perde o sentido** para o laço — e isso é bônus: some a
  espera ociosa dos trabalhadores e some o `sendjobs`, que é o que faz o RSS ir de
  313 MB para 581 MB.
- **O bug do `phase.cpp` continua sem correção** (`phase.cpp:312-313`, `i=0` passa
  e a linha 316 lê `phasea[-60]`). Hoje está sem gatilho. Se o port mexer em
  `kmin`, ele volta.

## Build

`nvcc` para os `.cu`, `mpic++` para o resto, link com `-lcudart`. **Não tente
compilar o programa inteiro com `nvcc`** — `-fpermissive` e C++98 de 1998 não vão
sobreviver ao front-end. Um `.cu` isolado com interface `extern "C"` e tipos
POD é o caminho.

---

## A figura final do paper: V0 × V5 × CUDA

O script `baseline/figura-v5.py` já produz V0 × V5 e é a base. Três coisas a
resolver quando o CUDA existir:

1. **O eixo x muda de sentido.** V0 e V5 são tempo × `np`; a versão CUDA roda com
   `np=1` e uma GPU, então não tem curva em `np`. **Proposta: linha horizontal
   anotada** atravessando o painel, com o texto dizendo `np=1 + 1 GPU`. A
   alternativa (barras do melhor tempo de cada versão) esconde a escalabilidade,
   que é metade da história.
2. **Remedir V0 e V5 na mesma janela do CUDA.** A figura atual traz no rodapé o
   aviso de que o ganho é limite superior por causa de ~5% de deriva entre
   janelas. **Para o paper isso não serve**: as três versões têm de sair da mesma
   campanha, na mesma janela, máquina verificadamente vazia. É ~40 min de
   máquina e resolve o problema de vez.
3. **A legenda tem de dizer o critério de validação.** V0 e V5 são byte a byte;
   CUDA é `max |Δχ| ≤ 1e-4`. Omitir isso num paper seria esconder a única
   diferença metodológica entre as barras.

Paleta: os slots categóricos já validados são azul `#2a78d6`, laranja `#eb6834`,
carmim `#c2255c` (o terceiro entra para o CUDA). Revalidar com
`validate_palette.js --pairs all` ao passar de 2 para 3 séries — com três linhas
qualquer par pode ser confundido, não só os adjacentes.

---

## O que sobrou da Fase 1

Achados que valem para o **host** e não foram colhidos, em ordem de tamanho:

1. **`onerotation` refaz 779× o que é fixo.** Em `alldblevent`, `cosa,sina,phia`
   dependem só do par e `cosb,sinb,phib` só do ponto — ver "A geometria". Na CPU
   isso são 31,6 milhões de `sqrt`+`atan2` redundantes por corrida. É o candidato
   a V7 mais promissório que este trabalho produziu, e **é independente de GPU**.
   Validação: `regressao.sh` byte a byte — deve passar, porque os valores são
   literalmente os mesmos, só calculados uma vez.
2. **`fhankelfaca` reconstrói 100 complexos por par** para amortizar 15 chamadas
   que leem 3 colunas (`am ∈ {0,1,2}`). Construir só o que se lê é ~60 entradas
   em vez de 100.
3. **A cópia de linha do `makerotation`** — medida e reprovada (V6), mas o
   meio-termo "copiar em vez de apontar" não foi testado. Ver `OTIMIZACAO.md`.

## Ordem, em uma linha

Fase 0 **feita** → teste do `beta` **feito** → Fase 1 **feita e validada** →
**Fase 2** → Fase 3 → campanha final das três versões na mesma janela → figura.
