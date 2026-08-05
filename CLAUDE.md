# MSCDATA — mapa do projeto

**MSCD 1.37** (Van Hove / LBNL, 1997-98): difração e dicroísmo de fotoelétrons.
C++ antigo, 69 arquivos, ~19 mil linhas, paralelo por MPI. Não é o MSCD de
fábrica: tem uma extensão **ATA** da UNICAMP (`Mscdrun::ATAevenelem`,
`mscdrunc.cpp:46`), desligada no `Cov0.txt` (`ATA=0`).

Repositório: <https://github.com/yosefschmidtA/MSCD_ATA_GPU>. **O git é do
usuário** — não rode `add`, `commit`, `tag` nem `push`; entregue o comando pronto.

**Antes de abrir qualquer `.cpp`, leia estes** — eles existem justamente
para economizar leitura de código:

- **`OTIMIZACAO.md`** — **o mais atual, comece por aqui**. Abre com a seção
  **COMO CONTINUAR**, que diz em que estado o código está, qual é a próxima ação
  pendente e o que já foi decidido e não deve ser refeito. Depois: diagnóstico,
  o que mudou em cada versão, medições por fase, escalabilidade, e os resultados
  negativos com o motivo.
- **`baseline/README.md`** — a linha de base congelada e o critério de validação
  (`baseline/regressao.sh`: as 787 linhas de intensidade têm de sair idênticas
  byte a byte).
- **`README.md`** — caderno de laboratório. A fatoração caro/barato da eq. (46)
  que organiza o programa inteiro, o grafo das equações, onde está (e onde não
  está) o paralelismo. **Atenção: os números dele são da configuração antiga**
  (profundidade 15, 205 átomos) e a descrição das fases do `symtrivert` é do
  código antes da reescrita.
- **`EQUACOES.md`** — as ~40 equações do manual (`MANUAL_MSCD_CEA.pdf`) mapeadas
  para arquivo:linha, nos dois sentidos.

## Estado atual (05/08/2026)

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

Campanha de 05/08: **V0 69,70 s → V4 44,05 s (1,58×)** em `np=12`.

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
make randmscd_parallel CPPFLAGS="-O3 -std=c++98 -w -fpermissive"
mpirun --use-hwthread-cpus -np 4 randmscd_parallel Cov0.txt
```

As duas flags são obrigatórias, não preferência:

- **`-fpermissive`** — `fcomplex.h:46` declara `friend Fcomplex polar(float,float=0)`.
  Argumento padrão em `friend` que não é definição: ilegal no padrão, aceito pelos
  compiladores de 1998, recusado pelo g++ 13. É o **único** erro em 69 arquivos.
- **`--use-hwthread-cpus`** — o Open MPI conta núcleos **físicos** (6 num
  i5-13420H), não threads (12), e desde a série 3.x recusa superalocar em vez de
  aceitar calado como a 1.10 fazia. Sem ela, `-np 10` morre em "not enough slots".

Referência de saída correta: `factors = 0.6710 0.8642`, curva em
`saida1Co-alterado-alexandre.txt`. **Os `factors` não servem de teste de
regressão** — mudaram só de 0.6719 0.8649 para 0.6710 0.8642 quando 783 das 787
linhas de intensidade mudaram, porque são fatores de escala do ajuste.

**Use `-np ≥ 4`, e evite `-np 6`.** Mínimo das repetições; V0 é de 04/08, V2 e V4
de 05/08 (`baseline/campanha-v4.sh`, build de produção sem `-DMSCDTIMER`):

| `-np` | V0 | V2 | V4 |
|------:|---------:|---------:|---------:|
| 1 | 146,49 s | 144,52 s | 132,45 s |
| 2 | 87,53 s | 77,49 s | 71,45 s |
| 4 | 64,26 s | 51,07 s | 45,93 s |
| 6 | 70,97 s | 53,29 s | 46,08 s |
| 8 | 69,06 s | 49,30 s | 44,67 s |
| 12 | 69,70 s | 47,92 s | **44,05 s** |

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
