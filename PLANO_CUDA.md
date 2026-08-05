# Plano do port de CUDA

Escrito em 05/08/2026, sobre o V5. **Nada de CUDA foi escrito ainda.** Tudo que
está aqui vem de medição, não de leitura do código — e a leitura já errou três
vezes neste projeto (a fase C do `symtrivert`, o V3, e o `allrotation` que eu
suspeitei ser gargalo e é 1,4%).

**Objetivo final:** a figura **V0 × V5 × CUDA** do paper. Ver a última seção.

---

## COMO CONTINUAR

**Próxima ação, decidida em 05/08/2026:** rodar o **teste do `beta` inteiro** da
Fase 1 (~15 min) e, confirmado o resultado, emendar direto na Fase 1.

O teste está descrito na seção "O atalho do `beta` inteiro". Em resumo: contador
em `makerotation` (`rotamat.cpp:299`) contando quantas chamadas vindas de
`alldblevent` têm `xa != floor(xa)`. Se der zero, `rotmatb` é a linha `2·|beta|`
de `rotmata` e as 500 interpolações por par somem — na GPU **e** no CPU (seria o
V6). **Não implemente por leitura**: três palpites de leitura já foram derrubados
neste projeto.

**Estado do código:** V5 aplicado e validado, nada de CUDA escrito. Build:

```bash
make randmscd_parallel CPPFLAGS="-O3 -std=c++98 -w -fpermissive -fopenmp"
mpirun --use-hwthread-cpus --bind-to none -np 1 randmscd_parallel Cov0.txt
```

**Fase 0 está FEITA** — o critério de validação já foi medido e decidido
(`max |Δχ| ≤ 1e-4`). Não refaça essa discussão; a seção "Fase 0" traz o número e
por que as duas alternativas foram descartadas.

---

## Onde estamos (05/08/2026, medido)

| | valor |
|---|---:|
| melhor tempo hoje, `np=1` | **122,02 s** |
| preparo serial | 5,4 s (4,4%) |
| laço dos 779 pontos | ~115 s (95,6%) |
| **teto de Amdahl se só o laço for para a GPU** | **22,3×** |

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

## Fase 1 — `alldblevent` → `evenelem` (57% do laço)

O alvo, com folga. `mscdrunc.cpp:753` varre **40 562 pares distintos**
(`ndbleven`), cada um chamando `evenelem` (`mscdrunc.cpp:22`) até 15 vezes, e
`evenelem` é a soma sobre `l` com `lnum=20`.

### O desenho é **bloco ↦ par**, não thread ↦ par

Isto não é preferência, é consequência de um fato que só aparece lendo as
funções chamadas: **`evenelem` parece chamar quatro funções puras de tabela e
nenhuma das quatro é pura.** Todas têm cache mutável keyed no argumento:

| chamada | onde | o que reconstrói quando o argumento muda |
|---|---|---|
| `rotharma` | `rotamat.cpp:445` | `rotmatb` inteiro: `lnum×lamdum` = 20×25 = **500** |
| `fhankelfaca` | `msfuncs.cpp:173` | `hankarg`: `lnum×cmnum` = 20×5 = **100 complexos** |
| `fsinexpa` | `phase.cpp:431` | keyed em `akin`, **que é constante** ⇒ nunca reconstrói |

Por par são ~700 operações de reconstrução contra ~5400 flops no laço de `al`.
Num kernel ingênuo (thread ↦ par) **cada thread reconstruiria 2,5 KB de tabela** —
desastre. Com bloco ↦ par, o bloco monta `rotmatb` e `hankarg` em *shared memory*
cooperativamente (~3 KB por bloco, folgado) e as 15×20 iterações leem de lá.

### O atalho do `beta` inteiro — **verificar primeiro, é barato**

`alldblevent:769` faz `beta=(float)floor((beta*10.0+0.5)/10.0)`, que é
`floor(beta+0.05)`: **`beta` é sempre inteiro** ali. Como `betanum=361`
(`rotamat.cpp:27`, porque `sizeof(int)>=4`), o `makerotation` calcula
`xa = |beta|·(361−1)/180 = |beta|·2`, um inteiro par, e o fator de interpolação
`(xa−i)` dá **exatamente zero**.

Se isso se confirmar, aquelas 500 interpolações por par são **uma cópia de linha
disfarçada**: `rotmatb` é a linha `2·|beta|` de `rotmata`, com uns sinais
trocados. Na GPU indexa-se `rotmata` direto e o custo some.

**E vale no CPU também** — é candidato a V6, independente de GPU.

**Teste (15 min):** contador em `makerotation` contando quantas chamadas vindas
de `alldblevent` têm `xa != floor(xa)`. Se for zero, o atalho é real.
**Não implemente por leitura** — a regra deste projeto existe porque três
palpites já foram derrubados.

### Como validar a Fase 1 isoladamente

Rodar o kernel e o `evenelem` da CPU lado a lado, comparando `devenelem` par a
par, **antes** de mexer em qualquer outro bloco. É o único momento em que o
kernel é testável sem depender de nada — não pule.

---

## Fase 2 — `allevendetec` (8,8%)

Cópia com escala complexa, `natoms²×radim` = 915 015 elementos, zero dependência.
Trivial. **A razão de fazer é manter os dados na placa** — descer para o host só
para isto anularia o ganho da Fase 1.

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

- **`Fcomplex` é `{float re, im}`** (`fcomplex.h:31`), mesmo layout de `float2` —
  mas `fcomplex.h:46` declara `friend Fcomplex polar(float,float=0)`, argumento
  padrão em `friend` que não é definição. É o que obriga `-fpermissive` no g++, e
  **vai barrar o front-end de device.** Escreva um tipo próprio para o device; não
  tente incluir `fcomplex.h` em `.cu`.
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

## Ordem, em uma linha

Fase 0 **feita** → teste do `beta` inteiro (15 min) → Fase 1 com validação
isolada → Fase 2 → Fase 3 → campanha final das três versões na mesma janela →
figura.
