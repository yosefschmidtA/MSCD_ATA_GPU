# MSCDATA — mapa do projeto

**MSCD 1.37** (Van Hove / LBNL, 1997-98): difração e dicroísmo de fotoelétrons.
C++ antigo, 69 arquivos, ~19 mil linhas, paralelo por MPI. Não é o MSCD de
fábrica: tem uma extensão **ATA** da UNICAMP (`Mscdrun::ATAevenelem`,
`mscdrunc.cpp:46`), desligada no `Cov0.txt` (`ATA=0`).

Repositório: <https://github.com/yosefschmidtA/MSCD_ATA_GPU>. **O git é do
usuário** — não rode `add`, `commit`, `tag` nem `push`; entregue o comando pronto.

**Antes de abrir qualquer `.cpp`, leia estes** — eles existem justamente
para economizar leitura de código:

- **`PLANO_CUDA.md`** — o plano do port. Traz o **critério de validação medido**
  (`max |Δχ| ≤ 1e-4`, piso de ruído do programa em 1,0×10⁻⁵), as quatro fases na
  ordem, e o desenho da figura final do paper. **Fases 0 e 1 estão FEITAS e
  validadas (05/08/2026); a próxima é a Fase 2.** A seção "Armadilhas
  específicas do port" tem quatro itens que **custaram depuração real** — leia
  antes de escrever qualquer `.cu`, especialmente a regra de que **1 ulp aqui é
  erro de índice, não de arredondamento**.

- **`OTIMIZACAO.md`** — **o mais atual, comece por aqui**. Abre com a seção
  **COMO CONTINUAR**, que diz em que estado o código está, qual é a próxima ação
  pendente e o que já foi decidido e não deve ser refeito. Depois: diagnóstico,
  o que mudou em cada versão, medições por fase, escalabilidade, e os resultados
  negativos com o motivo.
- **`baseline/README.md`** — a linha de base congelada e o critério de validação
  (`baseline/regressao.sh`: as 787 linhas de intensidade têm de sair idênticas
  byte a byte).
- **`README.md`** — caderno de laboratório, **atualizado em 05/08/2026**. A
  fatoração caro/barato da eq. (46) que organiza o programa inteiro, a tabela das
  versões, o perfil do laço por dentro, o grafo das equações, e onde está (e onde
  não está) o paralelismo. Tem blocos "Correção de 05/08/2026" marcando duas
  afirmações que o perfil desmentiu — leia-os, o erro é a lição.
- **`EQUACOES.md`** — as ~40 equações do manual (`MANUAL_MSCD_CEA.pdf`) mapeadas
  para arquivo:linha, nos dois sentidos.

## Estado atual (05/08/2026, noite)

**A Fase 3 do port de CUDA está feita e validada.** Contudo, antes na Fase 2 o tempo aumentou para **66,19 s** em `np=1`. Isso ocorreu porque a transferência PCIe aumentou levemente (copiávamos os 7,3 MB do `devendetec` de volta para a CPU). Na Fase 3:

- **Fase 3** (Summation para GPU): **CONCLUÍDO.** (05/08/2026). Tráfego pesado contornado! A GPU agora resolve todo o loop `m` com matriz esparsa e devolve apenas as linhas dos átomos emissores (31 KB/ponto). Tempo despencou de 66.19s (Fase 2) para fenomenais **37.69s**. Gargalo principal do PCIe aniquilado.
- **Fase 4** (pathcut na GPU): (Opcional, 1.7s).

## Atualização de 06/08/2026 — capacidade e o teto de latência

**O limite de átomos subiu de 300 para 1250** (`mscdruna.cpp:317` e `:368`,
`patom[1250*12]` em `mscdrun.cpp` e `mscdrun.h`). Não é otimização, é capacidade:
sem isso o `1x2iron.in` de 316 átomos dava erro 602. **O 1250 não é número
redondo à toa** — `tevenadd` é indexado por `ia*natoms²+ib*natoms+ic` em `int`, e
acima de **natoms ≈ 1290** isso estoura o int de 32 bits e corrompe em silêncio
(`OTIMIZACAO.md:1088`). 1250³ = 1,95e9 contra o teto de 2,15e9: sobra 9% de
folga. **Não suba mais esse limite sem passar o índice para 64 bits.**

**Cache do `tevenelem` na GPU** (`mscdgpu.cu:625`, `last_akin`): só reenvia o
array quando `akin` muda. Correto **porque a energia é fixa** neste modo
(`scanmode=223`, `kmin=kmax`). **Atenção:** se algum dia ligar ajuste de
geometria com energia fixa, `tevenelem` muda sem `akin` mudar e o cache serve
lixo. Hoje não há gatilho (`trymax=0`).

**O ganho do cache foi ~6 s, não o que se esperava** — a hipótese dos 240 GB de
PCIe estava errada, porque o `pathcut` já poda 99,8% dos caminhos e o que trafega
é minúsculo. **O gargalo real com 316 átomos é latência de lançamento de kernel**:
a GPU aparece com <1% de uso porque passa o tempo escalonando pedidos liliputianos
de vários ranks. É por isso que `np>4` piora. **A próxima etapa é *batching***:
resolver centenas de ângulos por kernel, em vez de um. Derivação na seção "Fase 4"
do `OTIMIZACAO.md`.

## Histórico de Armadilhas

**05/08/2026 - O "Falso Positivo" ZERADO na Fase 3**
- **O erro:** Durante a implementação da Fase 3 (`summation`), a versão CPU do loop foi completamente removida do bloco compilado com `-DMSCDGPU`, em vez de ter seu desvio em tempo de execução via `getenv("MSCD_GPU")`. 
- **O sintoma:** O script `./baseline/regressao-gpu.sh 1` passava com precisão perfeita, mas isso porque ele forçava a execução na GPU (`MSCD_GPU=1`). Quando o usuário rodou o binário manualmente *sem a variável de ambiente* (como exigido no fallback de CPU detalhado no `PLANO_CUDA.md`), a execução falhou silenciosamente e cuspiu zeros em todo o arquivo `.chi`.
- **A correção:** O usuário percebeu que o arquivo de saída gerado estava preenchido com zeros. A solução foi restaurar a rota original da CPU e inserir um bloco `if (getenv("MSCD_GPU"))` dentro do macro, garantindo que o binário suporte as duas vias em tempo de execução. Nunca remova o caminho CPU da função, ele deve coexistir!

**A GPU passou o V5 na Fase 3, e por larga margem.** *(Este parágrafo dizia o
contrário até 06/08/2026 — era verdade na Fase 1, quando só 57% do laço estava na
placa. A Fase 3 inverteu, e a versão antiga mandaria você otimizar o lado errado.)*

| `-np` | V5 CPU | Fase 3 GPU |
|---:|---:|---:|
| 1 | 122,02 s | 35,65 s |
| **4** | 40,50 s | **24,43 s** ← melhor marca do projeto |
| 12 | **38,77 s** | 37,79 s |

**`MSCD_GPU=1 mpirun -np 4` é o comando de produção.** O ótimo saiu de `np=12`
(CPU) para **`np=4`** (GPU), e a razão é arquitetural: com 4 ranks as CPUs se
revezam orquestrando a placa em *overlap* — um copia memória, outro dispara
kernel, outro faz a física serial. Com 12 ranks brigando pela mesma janela de
PCIe o escalonador serializa e o tempo volta para 37,79 s, empatando com o V5.
**Mais ranks pioram, não melhoram.** Medição e derivação em `OTIMIZACAO.md:1093`,
figura em `baseline/gpu-v0-v5.png`.

**A GPU continua não tocando no "Analyzing/Reanalyzing"** (`symtrivert`, do
preparo) — esse segue sendo alvo do V1/V2.

**V6 (`beta` inteiro no `makerotation`) foi medido e REVERTIDO** — byte a byte
idêntico, +1,5% em `np=1` e −3,5% em `np=12`. O achado que ele produziu é o que
sustenta o kernel; a otimização de CPU, não. **O nome V6 está queimado; a próxima
é V7.** Ver `OTIMIZACAO.md`, seção V6.

## Estado do CPU (05/08/2026, manhã)

Cluster de trabalho: **raio 9 / profundidade 20, 247 átomos** (`Cov0.txt`).
Três otimizações aplicadas, todas com curva idêntica byte a byte:

- **V1/V2** (`mscdrunb_not_reanalize.cpp`): `symtrivert` reescrito — fase C
  redundante eliminada, hash no lugar da busca linear, `std::sort` no lugar dos
  *selection sorts*.
- **V3**: tentativa de inverter os laços do `pathcut`. **Piorou 13%, revertida.**
  O nome está queimado — a próxima é V5.
- **V4** (`mscdrund.cpp:128`): o teste do `tevencut` movido para **antes** de
  encher o `csum`. 99,87% dos pares saem no `continue`, então encher antes era
  ler 15 complexos para jogar fora — ~199 GB de tráfego na corrida inteira.
- **V5** (`mscdrunc.cpp:394` e `:427`): dois `#pragma omp parallel for` no
  `pathcut`. **Laço de 7,27 s → 1,71 s (4,3×)**; total em `np=1` de 125,69 s
  para **120,44 s (4,2%)**, curva idêntica byte a byte nas seis rodadas.
  **Só vale com `--bind-to none`** — ver Armadilhas. Sem `-fopenmp` os pragmas
  somem e o binário é byte a byte igual ao V4 (verificado com `cmp`).
  **Não ajuda com `np>1`**: o preparo é só do rank 0 e os outros ranks giram em
  espera ativa, então não há núcleo livre.

Campanha de 05/08: **V0 69,70 s → V4 44,05 s (1,58×)** em `np=12`.

**Preparo serial cronometrado inteiro em 05/08** (`np=1`, máquina leve, ~10,8 s
antes do V5): `pathcut` 7,13 s (66%), `symtrivert` 1,28 s, `alltrievent` 1,18 s,
`symdblvert` 0,80 s, `allrotation` 0,15 s, resto 0,29 s. **`allrotation` não é
gargalo** — são 15 milhões de iterações com `atan2`/`acos`, mas o guarda
`if (eledim<2)` (`mscdrunc.cpp:742`) poda quase tudo. Hipótese testada e
descartada. **Com o V5 o preparo é ~5,4 s**, e o teto de Amdahl para o port de
GPU subiu de **12,2× para 22,3×** — sem uma linha de CUDA. O preparo deixou de
ser o assunto; o próximo item serial é o `symtrivert` com 1,28 s, que não vale o
risco.

**O laço dos 779 pontos foi perfilado por dentro em 05/08** e o resultado
contraria o que se supunha por leitura do código: **`alldblevent` é 57% do
laço**; o laço de `m` do `summation` é 22% e **quase não calcula física** —
99,87% das visitas são podadas pelo `tevencut`, sobram 1 455 trios com 14 724
MACs, e o custo é tráfego de memória. Detalhes em `OTIMIZACAO.md`, seção "Perfil
do laço paralelo".

**A configuração de física mudou em 04/08/2026.** O `ps01` era
`psAg111-slab.txt`, que cobre k 16,00–18,00 — e o `kconfine` (`phase.cpp:419`)
levantava em silêncio o k de 13,63 para 16,00, rodando o cálculo a 975 eV contra
dado experimental medido a 708 eV. Agora é **`psAg111.txt`** (mesma prata, k
5,00–17,75, `lnum=20`), e o log confirma `13.63 13.63`. Config antiga em
`Cov0.txt.k16-slab.bak`. **Nenhum número medido antes disso é comparável com os de
depois**, e a referência do `baseline/regressao.sh` está obsoleta. Derivação
completa na seção "O k errado" do `OTIMIZACAO.md`.

Compilar com `-DMSCDTIMER` liga os cronômetros do `mscdtimer.h`, que imprimem em
stderr; sem a macro não sobra instrução nenhuma. São dois níveis: cronômetros de
**fase** (`MSCDT`, um `fprintf` por marca) e **acumuladores** (`MSCDT_A`, que
somam e imprimem uma vez — obrigatório dentro do laço, onde `summation` roda 3895
vezes). Os acumuladores são `static` por unidade de tradução: marcar e relatar têm
de ficar no mesmo `.cpp`. **Build instrumentado não serve para medir tempo total**:
o `summation` faz um passe extra de histograma na primeira chamada (~1 s).

O único executável usado é o **`randmscd_parallel`**. O `makefile` constrói outros
doze; estão compilados e funcionando, e não interessam.

## Compilar e rodar

```bash
# CPU (produção)
make randmscd_parallel CPPFLAGS="-O3 -std=c++98 -w -fpermissive -fopenmp"
mpirun --use-hwthread-cpus --bind-to none -np 1 randmscd_parallel Cov0.txt

# GPU (Fase 3). O rm é obrigatório: os .o precisam de -DMSCDGPU, e o makefile
# não tem alvo clean nem sabe distinguir as duas variantes.
rm -f *.o && make randmscd_gpu \
  CPPFLAGS="-O3 -std=c++98 -w -fpermissive -fopenmp -DMSCDGPU"
./baseline/regressao-gpu.sh 1

# GPU em produção — np=4 é o ótimo medido, e np>4 PIORA (latência de kernel)
MSCD_GPU=1 mpirun --use-hwthread-cpus -np 4 randmscd_gpu Cov0.txt
```

**Os dois builds compartilham os `.o` e se atropelam.** Depois de mexer no
`randmscd_gpu`, `rm -f *.o` de novo antes de reconstruir o `randmscd_parallel`,
senão ele sai com `-DMSCDGPU` dentro.

Arquivos do port: **`mscdgpu.cu`** (kernel), **`mscdgpu.h`** (interface POD),
cola em **`mscdrunc.cpp:768`** (`gpudblevent`), despacho em
**`mscdrund.cpp:298`**. Acessores `gpu_*` em `rotamat.h`, `msfuncs.h`,
`phase.h`, `vibrate.h` — só leitura, existem para não incluir cabeçalho do
programa no `.cu`.

As duas primeiras flags são obrigatórias, não preferência:

- **`-fpermissive`** — `fcomplex.h:46` declara `friend Fcomplex polar(float,float=0)`.
  Argumento padrão em `friend` que não é definição: ilegal no padrão, aceito pelos
  compiladores de 1998, recusado pelo g++ 13. É o **único** erro em 69 arquivos.
- **`--use-hwthread-cpus`** — o Open MPI conta núcleos **físicos** (6 num
  i5-13420H), não threads (12), e desde a série 3.x recusa superalocar em vez de
  aceitar calado como a 1.10 fazia. Sem ela, `-np 10` morre em "not enough slots".

As outras duas são do V5: **`-fopenmp`** liga os pragmas do `pathcut` (sem ela o
binário é byte a byte idêntico ao V4), e **`--bind-to none`** é o que faz o
OpenMP existir de verdade — ver Armadilhas.

**`--bind-to none` só com `np=1`.** Medido em 05/08/2026: em `np=1` ela vale
5,25 s; em `np=12` os dois modos empatam (`padrao` 39,00 s, `bindnone` 42,08 s) —
o preparo é do rank 0 e os outros ranks ocupam os núcleos, então não há o que o
OpenMP pegar. **`--map-by slot:PE=$(nproc)` não serve para `np` alto**: com
`-np 12` ela pede 144 PEs e o `mpirun` morre em 0,03 s.

Referência de saída correta: `factors = 0.6710 0.8642`, curva em
`saida1Co-alterado-alexandre.txt`. **Os `factors` não servem de teste de
regressão** — mudaram só de 0.6719 0.8649 para 0.6710 0.8642 quando 783 das 787
linhas de intensidade mudaram, porque são fatores de escala do ajuste.

**Use `-np ≥ 4`, e evite `-np 6`.** Mínimo das repetições. **Quatro janelas
diferentes**: V0 de 04/08; V2 e V4 de 05/08 madrugada (`baseline/campanha-v4.sh`);
V5 de 05/08 tarde (`baseline/campanha-v5-escala.sh`). Build de produção sem
`-DMSCDTIMER`. **Compare colunas com cuidado** — a deriva entre janelas medida
entre a madrugada e a tarde de 05/08 foi de ~5% com o mesmo binário:

| `-np` | V0 | V2 | V4 | V5 |
|------:|---------:|---------:|---------:|---------:|
| 1 | 146,49 s | 144,52 s | 132,45 s | **122,02 s** |
| 2 | 87,53 s | 77,49 s | 71,45 s | 69,82 s |
| 4 | 64,26 s | 51,07 s | 45,93 s | 40,50 s |
| 6 | 70,97 s | 53,29 s | 46,08 s | 39,50 s |
| 8 | 69,06 s | 49,30 s | 44,67 s | 38,84 s |
| 12 | 69,70 s | 47,92 s | 44,05 s | **38,77 s** |

Figura V0 × V5 em `baseline/v0-v5.png` (gerada por `baseline/figura-v5.py`), com
o aviso da deriva no rodapé. **A coluna V5 em `np≥2` não é ganho do V5** — o
OpenMP não age lá; é a mesma física do V4 medida numa janela mais rápida.

**Não afirme um `np` ótimo.** Em 04/08 o `np=4` ganhava do `np=6` nos quatro
pareamentos e foi declarado ótimo; em 05/08 a ordem é `12 < 8 < 4 < 6` e o
`np=12` ganha do `np=4` nos quatro pareamentos, com o **mesmo argumento apontando
para o outro lado**. A diferença (~6%) é do tamanho da deriva entre janelas
(+0,8% em média, +6,5% no pior `np`). O que sobrevive às duas: `np=1` e `np=2`
são piores, e **`np=6` é consistentemente o pior de {4,6,8,12}** (V2 deu 53,28 e
53,29 s nas duas janelas). Pareamento detecta ordenação *dentro* de uma janela,
não a estabilidade dela entre janelas.

**O gargalo está dentro do laço, e tem nome.** O preparo serial é ~11,4 s de
132,45 s em `np=1` (8,6%) ⇒ teto de Amdahl 11,6×, medido 3,0×. Dentro do laço:
`alldblevent` 57%, laço de `m` 22%, bloco final 9%, `allevendetec` 9%.

*(Versões anteriores deste arquivo erraram aqui três vezes: `-np 10` = 41,4 s;
depois `-np 6` = 41,4 s; depois `-np 4` como ótimo estabelecido. As duas
primeiras mediram mal; a terceira mediu bem, mas generalizou uma janela só.)*

## Port de GPU — Fase 1 feita e validada

O objetivo do repositório (`MSCD_ATA_GPU`). **A Fase 1 existe e passa**
(`mscdgpu.cu`, `mscdgpu.h`, cola em `mscdrunc.cpp:768`). Detalhes e a próxima
ação em `PLANO_CUDA.md`. Em uma linha:

```bash
rm -f *.o && make randmscd_gpu \
  CPPFLAGS="-O3 -std=c++98 -w -fpermissive -fopenmp -DMSCDGPU"
./baseline/regressao-gpu.sh 1      # max|dchi| <= 1e-4
```

`MSCD_GPU=1` substitui a CPU, `MSCD_GPU=validate` roda as duas e compara par a
par, sem a variável roda 100% no host. Resultado medido: **`max|Δχ| = 1,0×10⁻⁵`,
que é o piso de ruído do próprio programa** — 10× abaixo do critério.

**Três coisas que custaram depuração e estão em `PLANO_CUDA.md`:** `-fmad=false`
é obrigatório no `nvcc`; a ordem das promoções float→double do código de 1998
tem de ser copiada literalmente (`sqrt(1.0-cosa*cosa)` faz o quadrado em
**float**); e **1 ulp aqui é erro de índice, não de arredondamento**, porque
`beta`, o `k` de `fexpix` e o `i` de `fhankelfaca` são inteiros tirados de
floats.

O resto desta seção é o que a medição de 05/08/2026 estabeleceu **antes** do
kernel existir, e continua valendo.

**Perfil do laço dos 779 pontos** (`np=1`, V4, acumuladores fechando em 99,7%):

| bloco | % do laço | onde |
|---|---:|---|
| **`alldblevent`** | **57,0%** | `mscdrunc.cpp:753` |
| `summation` — laço de `m` | 22,1% | `mscdrund.cpp:108` |
| `summation` — bloco final | 9,0% | `mscdrund.cpp:183` |
| `allevendetec` | 8,8% | `mscdrunc.cpp:879` |
| `summation` — init do `asum` | 3,1% | `mscdrund.cpp:97` |

**O alvo é `alldblevent` → `evenelem`** (`mscdrunc.cpp:22`): 40 562 pares
distintos, cada um chamando `evenelem` até 15 vezes, e `evenelem` é a soma sobre
`l` com `lnum=20`. Tarefas independentes, aritmética densa, escrita coalescida em
`devenelem[j*radim+…]`, tabelas pequenas (deslocamento de fase, Hankel,
harmônicos) que cabem em memória constante.

**Não perca tempo com o laço de `m`.** Ele parece o núcleo e não é: **99,87% das
90,2 milhões de visitas `(m,ia,ib,ic)` são podadas pelo `tevencut`**, sobram
1 455 trios com 14 724 MACs no total. O custo dele é tráfego de memória (as
cópias `bsum ← asum` e do `csum`), não conta. Metade já foi eliminada pelo V4 no
CPU; a outra metade também é problema de CPU, não de GPU. **A divergência de
`evedim` não é problema** — não há trabalho para divergir.

Fatos que o port pode usar:

- **A energia nunca muda** (`kmin=kmax=13,63`), então `alltrievent` roda uma vez
  e as tabelas de geometria são constantes na corrida inteira: sobem para a GPU
  **uma vez**. ~313 MB, ~25 ms em PCIe 4 — irrelevante contra ~125 s.
- **`Fcomplex` é `{float re, im}`** (`fcomplex.h:31`), mesmo layout de `float2`.
  Mas o `friend` com argumento padrão de `fcomplex.h:46` — o que exige
  `-fpermissive` — vai barrar o front-end de device: escreva um tipo próprio.
- **`tevenpar` é AoS de passo 10 floats** e o laço lê 2 campos: separar em arrays
  próprias é obrigatório. Já as quatro arrays `natoms³` são indexadas pelo mesmo
  `id` com `ic` variando mais rápido — **isso já é coalescido** se thread↦`ic`.
- **Amdahl**: preparo serial ~11,4 s de 132,45 s em `np=1` ⇒ **teto 11,6×** se só
  o laço for para a GPU. Para passar disso é preciso levar o `pathcut` junto
  (7,9 s dos 11,4 s). O `symtrivert` (1,4 s) não vale o risco.
- **Com uma GPU, `np>1` perde o sentido** para o laço — e isso é bônus: some a
  espera ociosa dos trabalhadores e some o `sendjobs` (que também é o que faz o
  RSS ir de 313 MB para 581 MB).

**A validação byte a byte não sobrevive ao port** — feito: o critério está em
`baseline/regressao-gpu.sh` (`max|Δχ| ≤ 1e-4`, `rms ≤ 1e-5`, sobre a coluna
`chical`). O `baseline/regressao.sh` byte a byte **continua valendo sem
alteração** para qualquer mudança de CPU e não deve ser relaxado. O R-factor e os
`factors` servem de sanidade, nunca de critério.

## Como medir neste projeto (custou tempo real aprender)

- **A carga da máquina quase dobra as fases seriais** (`symtrivert` 1,37 s com a
  máquina leve, 2,72 s com o navegador aberto). Confira `/proc/loadavg` e
  `ps -eo args | grep -cE "^mpirun|^randmscd"` antes de medir.
- **Normalize por um bloco intocado.** Foi assim que o V4 foi medido sem máquina
  dedicada: o `alldblevent`, que a mudança não tocou, variou 3,3% enquanto o laço
  de `m` caía 43%. Sem isso o tempo de parede não conclui nada.
- **Só compare medições da mesma janela.** O `escala.csv` tem coluna `janela` e o
  `escala.py` imprime a deriva do V2 entre elas justamente para isso.
- **Meça sempre o braço de controle na mesma janela, mesmo quando "já se tem o
  número".** Custa duas rodadas e é a única coisa que separa o efeito da deriva.
  Caso concreto do V5 (05/08/2026): o V4 tinha 132,45 s de `np=1` da janela
  anterior; contra ele o V5 pareceria ganhar **12,01 s (9,1%)**. Medido o
  controle na mesma janela (125,69 s), o ganho real é **5,25 s (4,2%)** — a
  janela inteira estava 5% mais rápida. Sem o controle o anúncio teria sido mais
  que o dobro do verdadeiro.
- **Pareamento não prova estabilidade.** Duas campanhas concluíram `np` ótimos
  diferentes com o mesmo teste. Ver "Use `-np ≥ 4`" acima.
- **Build com `-DMSCDTIMER` não serve para medir tempo total**: o `summation` faz
  um passe extra de histograma na primeira chamada (~1 s).

## Arquitetura do paralelismo

**Mestre–trabalhador com difusão do job**, não divisão de domínio:

```
mscdmain.cpp:13   mpiinit(); mype=mpigetmype(); numpe=mpigetnumpe()
mscdjob.cpp:108   Mscdjob::sendjobs()     — SÓ o rank 0 (mype==0)
                  paexport() serializa o job e envia a todos:
                  tamanho (tag 1024+i), depois o corpo (tag 2048+i)
mscdjob.cpp:132   Mscdjob::receivejobs()  — ranks > 0
                  paimport() reconstrói, e então dispmode=1, displog=0,
                  flogout=NULL — por isso só o rank 0 escreve arquivo
```

O código de física nunca chama `MPI_` direto: passa por uma camada fina
(`mpisend`, `mpireceive`, `mpigetmype`…) implementada em **`userCluster.cpp`**.
Foi ela que permitiu trocar de versão de MPI sem tocar em cálculo nenhum.

## Armadilhas

- **O Open MPI amarra o processo a um núcleo, e isso mata o OpenMP em silêncio.**
  Com `np` baixo o binding padrão é *bind-to core*: as 12 threads do
  `#pragma omp` do V5 ficam empilhadas num núcleo só. Medido em 05/08/2026: o
  laço do `pathcut` deu **7,269 s com binding padrão e 1,706 s com
  `--bind-to none`** — mesmo binário, mesma máquina, ganho zero contra 4,3×.
  Não há aviso nenhum. Se uma paralelização OpenMP "não fez nada", **confira o
  binding antes de culpar o código**: `mpirun --report-bindings`.
- **Se o programa girar a 100% de CPU sem escrever nada, é `phase.cpp`, não
  contenção.** `phase.cpp:312-313` deixa `i=0` passar (o guarda `if (i<0) i=0` é
  código morto; devia ser `if (i<1) i=1`), e a linha 316 então lê `phasea[-60]`,
  fora do array — `phasea` é `float*`, então é lixo do heap. Lixo grande faz o
  `while (xc-xb>90.0) xc-=180.0f` da linha 320 nunca terminar, porque `float`
  grande absorve a subtração de 180. Provado no *disassembly* em 04/08/2026
  (`rip = makephase+0x110`) depois de 52 min travado com a máquina vazia. É
  intermitente: depende do que estiver na memória antes do bloco. **Hoje está sem
  gatilho** (a tabela cobre o k de trabalho), **não corrigido.** Diagnóstico
  completo no `OTIMIZACAO.md`, seção "O `np=1` que não terminava".
  Para amostrar a pilha é preciso `sudo` — `ptrace_scope=1` só permite descendente:
  `sudo gdb -p PID -batch -ex "bt 25" -ex "info registers rip"`.
- **O `makefile` linka `mscdrunb_not_reanalize.o`, não `mscdrunb.o`.** O binário
  usa a variante modificada pelo Abner (tolerância de deduplicação geométrica de
  0,001 → 0,005 Å, limiar `natoms>100` → `>333`). **Editar `mscdrunb.cpp` não tem
  efeito nenhum no executável.** `diff` entre os dois mostra tudo.
- **`usercomp.cpp` e `usert3e.cpp` são cópias byte a byte** do `userCluster.cpp`.
  Os nomes sugerem variantes de plataforma; não são. Não compare os três.
- **O `core` na raiz é pista falsa**: ELF 32-bit i386, e o binário é x86-64. Veio
  na cópia da máquina antiga. `file core` antes de investigar.
- **`pgrep -x randmscd_parallel` dá 0 com o programa rodando** — o kernel trunca
  `/proc/PID/comm` em 15 chars e o nome tem 18. Use `pgrep -f`.
- **Sem alvo `clean` no `makefile`**, e os timestamps vieram todos iguais da cópia:
  para recompilar de verdade, tire os `.o` do caminho à mão.
- **Se voltar a faltar `libmpi.so.12`**: é ABI de Open MPI 1.10 contra o 4.1.6 da
  máquina, não pacote faltando. Recompilar resolve; instalar não.

## O grafo

907 nós, 1744 arestas em `graphify-out/`, construído em 2,6 s sem gastar token.
**388 arestas `calls` cruzam arquivos** — aqui o grafo serve para navegar.

Cobertura de `calls` é **parcial**: `mpisend()` aparece com grau 1 e é chamado 6×
de `mscdjob.cpp`. Grafo para orientar, `grep` para confirmar, `Read` estreito.

## Backups

`backup-original-20260804.tar.gz` e `.obj-antigos/` **não estão mais neste
diretório** (verificado em 04/08/2026). Se precisar dos binários originais,
procure em outro lugar antes de assumir que existem aqui. O código-fonte está
versionado no GitHub desde a tag/commit inicial, então o risco real é baixo.
