# Linha de base — Ag(111), raio 9 / profundidade 20

> Campanha de otimização e resultados por versão: **`../OTIMIZACAO.md`**.
> `v1/` e `v2/` guardam curva, log e fonte de cada versão (55,6 s e 48,4 s,
> curvas idênticas). Figura: `antes-depois.png`, gerada por `figura.py`.

Referência congelada em 04/08/2026 para validar qualquer otimização do
`symtrivert` ("Analyzing/Reanalyzing") e do resto do preparo serial.

## O critério

> **OBSOLETO desde 04/08/2026 20:00 — a base tem de ser regerada antes de o
> critério voltar a valer.** A referência guardada aqui foi calculada em
> k = 16,00 Å⁻¹, com o `ps01 = psAg111-slab.txt`. Descobriu-se que aquela tabela
> de phase shift não cobria o k do experimento (Co 2p, 13,63 Å⁻¹) e que o
> `kconfine` levantava o k em silêncio; o `Cov0.txt` agora usa `psAg111.txt` e
> **783 das 787 linhas mudaram**. Nada sai igual ao que está aqui. Ver "O k
> errado" no `../OTIMIZACAO.md`. Os arquivos deste diretório ficam preservados
> como registro da configuração antiga — `Cov0.txt` daqui é a config de k = 16.

**As 787 linhas de intensidade de `saida1Co-alterado-alexandre.txt` têm de sair
idênticas byte a byte.** Qualquer modificação que mude uma casa decimal mudou o
resultado físico e está errada até prova em contrário.

```bash
./baseline/regressao.sh        # np=6
./baseline/regressao.sh 1      # np=1
```

O script descarta só o rodapé (`This calculation took … starting on …`), que traz
data e duração — é a única coisa que varia entre rodadas idênticas. Verificado
em 5 rodadas.

Os `factors = 0.6719 0.8649` do terminal **não** servem sozinhos como critério:
saem com 4 casas (`mscdrund.cpp:470`), grosso demais para pegar deriva numérica.
Servem de conferência rápida; quem manda é a curva.

## Os números

Entrada: `Cov0.txt` com `9 20 4.086 radius,depth,lattice` (a cópia que gerou esta
base está aqui ao lado). Binário `randmscd_parallel` de 04/08 11:49, sem
modificação. `mpirun --use-hwthread-cpus -np 6`.

| | |
|---|---|
| factors (R-factor) | **0.6719 0.8649** |
| tempo, regime | **70 ± 2 s** (68,0 / 69,9 / 70,3 / 71,9) |
| tempo, 1ª rodada fria | 100,5 s — **descartar** |
| pico de RSS | 582 MB |
| `natoms` | 247 (era 205 com profundidade 15) |
| `nsymm` final | 1525 (inicial 3050, **1 Reanalyzing**) |
| `ntrieven` | 823 318 trios distintos, de 247³ = 15 069 223 |
| `ndbleven` | 40 562 |

Dispersão de 3% entre rodadas: **ganho abaixo de ~5% é ruído**, não conclua nada
com uma rodada só.

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

## Onde vão os 70 s

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
do `Cov0.txt`) e um arquivo versionado. **Toda rodada sobrescreve ele.** A versão
no HEAD é da configuração antiga (profundidade 15, 205 átomos, factors
0.6724 0.8647) — por isso ele aparece modificado no `git status`.

As medições do `README.md` da raiz (tabela de `-np`, `symtrivert` = 9,4 s) são da
configuração de profundidade 15 e **não** comparam com esta.
