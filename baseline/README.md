# Linha de base — Ag(111), raio 9 / profundidade 20

> Campanha de otimização e resultados por versão: **`../OTIMIZACAO.md`**.
> `v1/` e `v2/` guardam curva, log e fonte de cada versão. Figuras:
> `antes-depois.png` (`figura.py`) e `escalabilidade.png` (`escala.py`).

**Base regerada em 05/08/2026** na configuração corrigida de k = 13,63
(`ps01 = psAg111.txt`, `lnum=20`). A base anterior, de k = 16,00, está
preservada em **`k16-slab/`** — é registro histórico, não compare com ela.

## O critério

**As 787 linhas de intensidade de `saida1Co-alterado-alexandre.txt` têm de sair
idênticas byte a byte.** Qualquer modificação que mude uma casa decimal mudou o
resultado físico e está errada até prova em contrário.

```bash
./baseline/regressao.sh        # np=6
./baseline/regressao.sh 1      # np=1
```

O script descarta as linhas de data e duração — o rodapé (`This calculation took
… starting on …`) e o cabeçalho (`calculated by … on <data>`). São as únicas que
variam entre rodadas idênticas. *(O cabeçalho entrou no filtro em 05/08/2026:
sem ele a regressão acusa 1 linha de diferença sempre que a base e a rodada caem
em dias diferentes.)*

Os `factors` do terminal **não** servem sozinhos como critério: saem com 4 casas
(`mscdrund.cpp:530`), grosso demais para pegar deriva numérica. Prova disso: a
troca do `ps01` mudou **783 das 787 linhas** e os fatores se moveram só de
0,6724 0,8642 para 0,6710 0,8642, porque são fatores de escala do ajuste.
Servem de conferência rápida; **quem manda é a curva**.

> **Isto vai quebrar no port de GPU.** A redução sobre `ic` soma em outra ordem,
> e em `float` isso mexe no último bit. O critério byte a byte terá de virar
> tolerância relativa — decidir o limiar **antes** de escrever kernel, com a
> curva do V4 na mão.

## Os números

Entrada: `Cov0.txt` com `9 20 4.086 radius,depth,lattice` e `ps01 psAg111.txt`
(a cópia que gerou esta base está aqui ao lado). Curva de referência:
`saida.txt`, o V2 do commit 5cac3b1, reproduzido byte a byte em 05/08.

| | |
|---|---|
| factors (R-factor) | **0.6710 0.8642** |
| pico de RSS | **313 MB** (`/usr/bin/time -f %M`, np=1) |
| `natoms` | 247 |
| `nsymm` final | 1525 (inicial 247²/20 = 3050, **1 Reanalyzing**) |
| `ntrieven` | **795 605** trios distintos, de 247³ = 15 069 223 |
| `ndbleven` | 40 562 |
| `npoint` | 779 (19 θ × 41 φ), × 5 sub-direções = 3895 solves |

O `ntrieven` **depende do k**: `mscdrunb_not_reanalize.cpp:411` faz
`vlenc = kmin/100.0f`, e esse `vlenc` entra na largura do bin de deduplicação.
Era 823 318 em k = 16,00. A dedup do `symtrivert` não é puramente geométrica.

Dispersão típica entre rodadas na campanha fria: 0,4% em `np` baixo, até 7% em
`np` alto. **Ganho abaixo de ~5% em `np` alto é ruído**, não conclua nada com
uma rodada só.

## As campanhas

| script | o que compara | quando |
|---|---|---|
| `campanha.sh` | V0 × V2, todos os `np` | 04/08/2026 20:47 |
| `campanha-v4.sh` | V2 × V4, todos os `np` | 05/08/2026 01:0x |

Ambas **apendam** em `escala.csv`, que tem uma coluna `janela` (a data). O
`escala.py` usa a janela mais recente de cada versão e **imprime a deriva do V2
entre as duas** — é isso que diz se comparar versões medidas em dias diferentes
é legítimo. Na prática deu +0,8% em média, +6,5% no pior `np`.

Binários preservados, um por versão: `randmscd_parallel.baseline` (V0 original),
`.v2` e `.v4` (build de produção de 05/08). **Não recompile por cima.**

## Como medir sem se enganar

Três erros cometidos nesta campanha, todos com custo real de tempo:

1. **Máquina ocupada.** Confira com
   `ps -eo args | grep -cE "^mpirun|^randmscd"` **antes** de medir; tem de dar 0.
   E não dispare dois trabalhos de medição sobrepostos (eu disparei).
   *(Este item dizia que a rodada de `np=1` de 897 s foi contenção de CPU — 6,0×,
   "exatamente a razão esperada". **Estava errado**: era o laço infinito do
   `phase.cpp:320`. A coincidência aritmética é que deu credibilidade à
   explicação errada e fechou a investigação cedo.)*
2. **Deriva sob carga contínua.** Depois de ~1,5 h medindo, a mesma versão sem
   mudança de código passou de 48,4 s para 56,5 s (~15%). **Só compare medições
   feitas na mesma janela**, ou deixe a máquina descansar antes de uma campanha.
   Resolvido na prática: a campanha fria de 04/08/2026 começou com `load 0.04` e
   deu 0,4% de dispersão entre repetições em `np` baixo.
3. **Buffer de saída não é progresso — mas parada prolongada é sintoma.** `stdout`
   redirecionado é bufferizado em bloco, então o arquivo fica atrás do programa.
   Só que **esta regra, aplicada cega, esconde travamento de verdade**: em
   04/08/2026 as "24 linhas paradas em 0,00%" foram descartadas como buffer e eram
   sintoma real de um laço infinito. O que decide não é o arquivo, é o processo:
   `awk '{print $14}' /proc/PID/stat` subindo com o arquivo parado por dezenas de
   minutos = girando em laço. Aí **amostre a pilha antes de matar**, que a
   evidência morre com o processo:
   `sudo gdb -p PID -batch -ex "bt 25" -ex "info registers rip"`.

## Onde vão os 70 s — **do V0, na configuração de k = 16**

> Seção histórica: descreve o V0 antes do V1/V2/V4 e antes da correção do k.
> Guardada porque explica **por que** a fase C do `symtrivert` foi atacada
> primeiro. **Não use estes números.** Os atuais estão no `../OTIMIZACAO.md`:
> preparo serial ~11 s com a máquina leve, e o laço — que aqui era uma caixa
> preta de 30,7 s — está perfilado por dentro desde 05/08/2026, com o
> `alldblevent` respondendo por 57% dele no V4 (49,6% quando medido no V2 —
> a fatia dele cresce à medida que o resto do laço encolhe).

Build com `-DMSCDTIMER` (ver `mscdtimer.h`), `-np 6`. A rodada instrumentada deu
70,42 s e curva idêntica — os cronômetros não custam nada.

| etapa (rank 0) | tempo | |
|---|---:|---|
| leitura de parâmetros | 0,01 s | |
| **`symtrivert`** | **20,87 s** | |
| ├ A — dedup `natoms³` | 6,75 s | inclui 1 Reanalyzing, ou seja ~metade é retrabalho |
| ├ B — ordena/compacta baldes | 1,02 s | |
| ├ C — indexa `natoms³` | **13,09 s** | recalcula a assinatura e busca linear no balde |
| └ D — estatísticas | 0,00 s | |
| `symdblvert` | 0,97 s | |
| `precutable` | 11,20 s | |
| `sendjobs` | 0,88 s | |
| **preparo serial** | **33,9 s** | **48% do total** |
| laço dos 779 pontos | 30,7 s | paralelo |

Teto de Amdahl com o preparo serial deste jeito: **2,1×**, por mais núcleo que se
jogue no problema.

A surpresa é a fase **C**, que sozinha custa o dobro da A. As duas percorrem os
mesmos 15 milhões de trios e calculam a mesma assinatura; a C ainda precisa
achar o índice com busca linear num balde de ~940 entradas de passo 40 bytes
(estoura a L1, enquanto o array da A tem passo 12 bytes). Ela é redundante: se a
fase A guardar o id da assinatura em `tevenadd` — array que já existe e já tem
exatamente `natoms³` ints — a C vira uma passada de remapeamento.

## Atenção

`saida1Co-alterado-alexandre.txt` é ao mesmo tempo o arquivo de saída (campo `pe`
do `Cov0.txt`) e um arquivo versionado. **Toda rodada sobrescreve ele** — por
isso ele aparece modificado no `git status` depois de qualquer execução. As
linhas que mudam entre rodadas idênticas são só as de data e duração; se mudar
mais que isso, a física mudou.
