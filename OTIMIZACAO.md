# Otimização do preparo serial do MSCD

Caderno da campanha de otimização do `symtrivert` ("Analyzing/Reanalyzing") e do
resto do preparo serial. Iniciado em 04/08/2026, cluster Ag(111) raio 9 /
profundidade 20, 247 átomos.

---

## COMO CONTINUAR (estado em 04/08/2026, fim da sessão)

### Onde o código está

- **V1 e V2 aplicados e validados** (`mscdrunb_not_reanalize.cpp`). **V3 foi
  tentado e revertido** — ver a seção V3; ficou um comentário no código avisando.
- O binário `randmscd_parallel` está compilado **com `-DMSCDTIMER`**, então
  imprime cronômetros de fase em stderr. Para reconstruir:
  ```bash
  make randmscd_parallel CPPFLAGS="-O3 -std=c++98 -w -fpermissive -DMSCDTIMER"
  ```
  Sem a macro o binário fica idêntico ao de produção (a instrumentação some).
- **`baseline/randmscd_parallel.baseline` é o binário V0 original preservado.**
  Não recompile por cima; é a referência de comparação.
- Validação de qualquer mudança: `./baseline/regressao.sh [np]` — as 787 linhas
  de intensidade têm de sair idênticas byte a byte. **A referência atual está
  obsoleta** (foi gerada em k = 16,00): ver "O k errado".
- **A configuração mudou**: `ps01` agora é `psAg111.txt` (k 13,63, `lnum=20`). A
  anterior está em `Cov0.txt.k16-slab.bak`. **Nenhum número medido antes de
  04/08/2026 20:47 é comparável com os de depois.**

### Feito: campanha de medição a frio

Rodada em 04/08/2026 20:47, com a configuração corrigida. Resultados na seção
"Escalabilidade". Em uma linha: **`-np 4`, V2 em 49,84 s contra 64,26 s do V0
(1,29×)**, reprodutibilidade de 0,4% nos pontos de `np` baixo.

```bash
./baseline/campanha.sh          # 2 repetições, pausa de 45 s, ~72 min
```

O script embute as precauções que custaram tempo real: aborta se achar outro
`mpirun`/`randmscd` rodando, descarta uma rodada de aquecimento, e pausa entre
medições. Grava `baseline/escala.csv` e chama `baseline/escala.py`, que regenera
`baseline/escalabilidade.png`. **Rodar com a máquina fria** — a campanha limpa
começou com `load 0.04`, e é a isso que se deve a dispersão de 0,4%.

### Próximas ações, nesta ordem

1. **Regerar a base de regressão.** Sem isso o `regressao.sh` acusa diferença em
   toda mudança, porque a referência é de k = 16,00. É o primeiro passo, antes de
   mexer em qualquer código.
2. **Fechar o bug do `phase.cpp`.** `if (i<1) i=1` na linha 313 (o guarda atual é
   código morto) e um limite de iteração nos quatro `while` de 320-325. Hoje o bug
   está sem gatilho, não corrigido — ver "O `np=1` que não terminava".
3. **Descobrir por que o laço paralelo perde 44%.** Mudou de prioridade: com o
   serial em 11,4% o teto de Amdahl é 8,8×, mas o medido é 2,72×, porque o laço
   escala 3,4× em 6 ranks. **Continuar cortando serial rende pouco agora.**
   Instrumentar o laço por `np` antes de propor qualquer coisa.
4. **`precutable` com profiling de verdade.** Segue sendo o maior item serial
   (9,7 s de 13,1 s na config antiga; 13,2 s de 15,5 s na nova), 82% num laço só.
   **Não otimizar por leitura**: dois palpites já foram derrubados pela medição (a
   fase C, que eu estimei igual à A e era o dobro; e o V3, que eu previ mais
   rápido e ficou 13% mais lento). Usar `perf stat -e
   cache-misses,LLC-load-misses` antes de propor transformação. *(`perf` não está
   instalado nesta máquina — verificado em 04/08/2026.)*
5. **`sendjobs` sem `MPI_Bcast`** (`mscdrun.cpp:502`): envia ~190 MB num laço
   ponto a ponto, custo ∝ número de ranks. 0,97 s com 6, ~10 s com 64. É o
   defeito que mais atrapalha em máquina grande.

### Regras de medição (custaram tempo real para aprender)

Estão em `baseline/README.md`, seção "Como medir sem se enganar". Resumo:
máquina verificadamente vazia antes de medir (`ps -eo args | grep -cE
"^mpirun|^randmscd"` tem de dar 0); nunca dois trabalhos de medição
sobrepostos; só comparar medições da mesma janela; e buffer de saída não é
progresso.

### O que foi decidido e não deve ser refeito sem motivo

- **Eliminar o Reanalyzing**: descartado *após medir*. Com o hash, economiza
  ~0,7 s (1,4%, abaixo do ruído) e exige a parte de maior risco do plano. O
  roteiro está registrado na seção V2 caso volte a valer.
- **`--mca mpi_yield_when_idle 1`**: não ajuda (69,87 s contra 70,42 s).
- **Inverter os laços do pathcut (V3)**: piorou 13%. Não repetir sem medir.

---

Linha de base, critério de validação e medições por fase: **`baseline/README.md`**.
Mapa das equações: `EQUACOES.md`. Arquitetura geral: `README.md`.

## O k errado — achado de física, 04/08/2026

Encontrado ao investigar o travamento do `np=1`, e não tem nada a ver com
otimização: **o programa vinha rodando na energia errada.**

`Cov0.txt` pede k = 13,63 Å⁻¹. O log dizia `16.00 16.00 0.00  kmin kmax kstep`.

A causa é o `kconfine` (`phase.cpp:419`), que eleva o `kmin` até o primeiro ponto
tabelado do phase shift: `if (*kmin<xa) *kmin=xa` com `xa=phasea[0]`. O `ps01`
era `psAg111-slab.txt`, que cobre **k 16,00–18,00** — faixa do Ag 3d, não do
Co 2p. Então o k era levantado em silêncio, sem aviso no log nem código de erro.

Que 13,63 é o valor correto se verifica por dois caminhos independentes:

- **Cinemática.** `linitial=1` ⇒ `"2p"` (`mscdruna.cpp:405`). Co 2p₃∕₂ (BE 778,1
  eV) com Al Kα (1486,6 eV) dá KE = 708,5 eV ⇒ k = 0,512331·√708,5 = **13,637**.
- **O dado experimental.** `Cobalto108.mscd` carrega o k na linha de cabeçalho da
  curva: `1  41  13.6300  18.0000  1.00000  0.00000`. O `18.0` é o theta, que
  bate com `dthetamin=18` do `Cov0.txt`, e as linhas de dados começam em
  `111.000`, que é o `dphimin=111` — a leitura dos campos está confirmada.

Ou seja: cálculo a **975 eV** ajustado contra dado medido a **708 eV**.

E a faixa 16–18 tem dono identificável: Ag 3d₅∕₂ (BE 368,3) com Al Kα dá k =
17,13, no meio do intervalo. O `Cov0.txt` ainda traz `sn  Ag(111)-3d` como nome
do sistema. O arquivo é sobra do experimento de Ag 3d, reaproveitado num
experimento de Co 2p.

**Correção aplicada:** `ps01` passou a `psAg111.txt` — a **mesma prata**, cobrindo
k 5,00–17,75 com 256 pontos e `lnum=20`. Confirmado no log: `13.63 13.63`. A
configuração antiga está em `Cov0.txt.k16-slab.bak`.

Consequências:

- **783 das 787 linhas de intensidade mudaram.** Os `factors` mexeram pouco
  (0.6719 0.8649 → 0.6710 0.8642) porque são fatores de escala do ajuste, não
  observável — **não use `factors` como teste de que a física mudou.**
- `lnum` dobrou (10 → 20): laço paralelo +42%, preparo serial +18%, total +8 a
  17%. Barato pelo que se ganha.
- Com 13,63 no meio da tabela, o `i=0` do `phase.cpp:312` deixa de ser alcançado
  — foi o que parou os travamentos.

**Pendente, decisão de física:** o `ps02` é `psl9.txt`, arquivo genérico
(`1  1 phase shift data`, `lnum=9`), enquanto existe `psCoHCP.txt` de cobalto na
mesma árvore (`~/Ag_Co/theory/`), cobrindo 5,00–17,75. Não foi trocado.

## A regra do jogo

**Toda versão tem de reproduzir as 787 linhas de intensidade byte a byte.**
`./baseline/regressao.sh` verifica. Nenhuma otimização aqui troca precisão por
velocidade — todas eliminam trabalho redundante ou trocam a estrutura de dados,
mantendo a mesma aritmética na mesma ordem.

**Atenção:** a referência do `regressao.sh` foi gerada em k = 16,00, antes da
correção acima. Depois dela nada sai igual àquela referência — **a base tem de
ser regerada** antes de o critério voltar a valer.

Isso não é preciosismo: os laços a jusante percorrem `tevenpar` na ordem do
índice (`mscdrunc.cpp:166,358,462,486`, `mscdrund.cpp:605`) e acumulam em ponto
flutuante. **Reordenar `tevenpar` muda o resultado na última casa.** Por isso
toda versão preserva a ordem exata do original, mesmo quando a ordem em si é
arbitrária.

## O diagnóstico

`symtrivert` não é busca de simetria no sentido de teoria de grupos. É
**deduplicação**: reduzir os `natoms³` trios ordenados (emissor, espalhador₁,
espalhador₂) ao conjunto das assinaturas geométricas distintas
`(r₁, 1/r₂, cos β, espécie)` — 15 069 223 → 823 318 no cluster atual.

A tabela hash é feita à mão, com baldes de capacidade fixa, e tem três fases:

| fase | linhas | o que faz |
|---|---|---|
| **A** | 385-454 | varre `natoms³`, calcula a assinatura, procura no balde por **varredura linear**, insere se não achou |
| **B** | 457-511 | ordena cada balde por (cos β, 1/r₂, r₁) com *selection sort* O(n²) e compacta tudo num array contíguo |
| **C** | 534-575 | varre `natoms³` **de novo**, recalcula a mesma assinatura, e acha o índice final por busca linear no balde |

O "Reanalyzing" é o estouro de um balde: `nsymm/=2`, `ia=ib=ic=natoms+10`, e
**tudo recomeça do zero** (linha 449). No cluster atual acontece 1 vez, ou seja,
metade do trabalho da fase A é retrabalho.

Custo total ≈ `natoms³ × ocupação_do_balde`. Como o número de distintos cresce
mais rápido que `natoms²`, o total escala perto de `natoms³·⁸` — é por isso que
interfaces complicadas explodem.

## Medições

`-np 6`, build com `-DMSCDTIMER`. Tempo em regime, descartando a primeira rodada
(fria). Dispersão entre rodadas: 3% — **ganho abaixo de ~5% é ruído**.

| etapa (rank 0) | V0 base | V1 | V2 |
|---|---:|---:|---:|
| `symtrivert` A — dedup `natoms³` | 6,75 s | 7,43 s | **1,38 s** |
| `symtrivert` B — ordena baldes | 1,02 s | 1,07 s | **0,05 s** |
| `symtrivert` C — indexa `natoms³` | 13,09 s | **0,19 s** | 0,16 s |
| **`symtrivert` total** | **20,87 s** | **8,70 s** | **1,59 s** |
| `symdblvert` | 0,97 s | 0,90 s | 0,92 s |
| `precutable` | 11,20 s | 10,89 s | 9,67 s |
| `sendjobs` | 0,88 s | 0,98 s | 0,94 s |
| **preparo serial** | **33,9 s** | **21,5 s** | **13,1 s** |
| laço dos 779 pontos (paralelo) | 30,7 s | 30,7 s | 30,7 s |
| **total, em regime** | **70 ± 2 s** | **55,6 ± 1,2 s** | **48,4 ± 1,8 s** |
| pico de RSS | 582 MB | 582 MB | 582 MB |

Rodadas: V1 = 55,76 / 54,28 / 56,67 s. V2 = 50,38 / 47,22 / 47,46 s.

**`symtrivert` 20,87 → 1,59 s (13×). Total 70 → 48,4 s (−31%).** Em todas as
versões: curva idêntica byte a byte, `nsymm=1525` e `ntrieven=823318` idênticos.

O preparo serial caiu de 48% para 27% do tempo total. O teto de Amdahl subiu de
2,1× para 3,7×.

![antes e depois](baseline/antes-depois.png)

Figura gerada por `baseline/figura.py` a partir das duas curvas gravadas.
**Resíduo máximo entre V0 e V1: 0 exato**, nos 779 pontos — não é "concordância
dentro da tolerância", é o mesmo número em todos eles.

---

## V1 — fundir a fase C na fase A

**A observação.** As fases A e C percorrem os mesmos 15 milhões de trios e
calculam exatamente a mesma assinatura. A fase A **já descobre** em que slot do
balde cada trio cai — e joga a informação fora:

```c
if (...casou...) j=n+tempadd[m]+10;   // linha 408 do original: o slot é perdido
```

A fase C existe só para redescobrir esse índice, e custa o dobro da A (13,1 s
contra 6,8 s) porque busca em `tevenpar`, cujo passo é de 40 bytes (`[j*10+k]`)
e estoura a L1, enquanto o array da fase A tem passo de 12 bytes.

**A mudança.**

1. `tevenadd` (que já existe e já tem exatamente `natoms³` ints) passa a ser
   alocado **antes** da fase A, e a fase A grava nele o slot de cada trio.
2. Dois arrays novos de `nsymm*tscatter` ints: `origem[]` acompanha de onde veio
   cada entrada enquanto a fase B embaralha os baldes; `destino[]` traduz
   slot-da-fase-A → índice compactado final, preenchido durante a compactação.
3. A fase C vira `tevenadd[n] = destino[tevenadd[n]]` — uma passada de
   remapeamento, sem geometria e sem busca.

**Por que o resultado não muda.** A fase B só permuta entradas; `origem[]`
registra a permutação exata. O conjunto, a ordenação e a compactação continuam
idênticos, byte a byte. O erro 621 (assinatura não encontrada) passa a ser
disparado por `destino[slot] < 0`, que é a mesma condição.

**Custo em memória: nenhum, medido.** Eram esperados +100 MB (os 40 MB de
`origem`+`destino`, mais 60 MB porque `tevenadd` agora convive com `tempar` em
vez de ser alocado depois que ele morre). O pico de RSS não mudou: 582 MB nas
duas versões. O pico do processo é fixado em outro ponto do programa, não aqui.

**O que custou.** A fase A subiu de 6,75 s para 7,43 s (+10%): são 15 milhões de
escritas espalhadas num array de 60 MB, que a versão antiga não fazia. Paga-se
0,7 s na A para economizar 12,9 s na C.

**Resultado medido.** `symtrivert` 20,87 → 8,70 s. Total 70 → 55,6 s, **−20,6%**.
Curva idêntica byte a byte.

**Bug corrigido de brinde.** O original refazia `new tempadd`/`new tempar` a cada
Reanalyzing **sem liberar os anteriores** (linha 381): ~130 MB vazados por volta.
Agora libera antes de realocar.

---

## V2 — hash de verdade na fase A, `std::sort` na fase B

Feito. **Fase A 7,43 → 1,38 s. Fase B 1,07 → 0,05 s.** Curva idêntica.

### O que foi feito

**V2a — hash no lugar da varredura linear (fase A).** A busca `for (j=n;
j<n+tempadd[m]; ++j)` percorria o balde inteiro — em média ~270 comparações por
trio, 15 milhões de vezes. Trocada por sondagem linear numa tabela de
endereçamento aberto.

Duas coisas tornam o hash **exato**, não aproximado:

1. As três chaves já chegam quantizadas (`round(...,0.005f)` e
   `round(...,vlenc)`), então padrão de bits idêntico ⟺ valores iguais, que é
   exatamente a comparação `==` que o original faz.
2. Única armadilha: `-0.0f` e `+0.0f` são `==` mas têm bits diferentes. Sem
   normalizar, o mesmo valor entraria duas vezes. `mscd_chave()` normaliza.

**A correção de rota que o código me obrigou a fazer.** Meu plano dizia que o
conjunto distinto não dependia de `nsymm`. Está errado, e a leitura atenta da
linha 424 mostra por quê: o balde é escolhido com `search`, que **inclui
`akind`** (a espécie do átomo do meio), mas a chave guardada em `tempar` é só
`(lena, vlenb, cosbeta)`. Duas espécies com a mesma geometria caem em baldes
diferentes e são **duas entradas distintas** no original. Por isso a chave do
hash inclui o balde `m` — sem isso elas seriam fundidas e a física mudaria.

O balde de uma entrada é recuperado de graça do próprio slot (`j/tscatter`), sem
array extra.

**V2b — `std::sort` na fase B.** Eram três *selection sorts* encadeados, O(n²)
cada, para chegar na ordem lexicográfica por (cos β, 1/r₂, r₁). `std::sort` dá o
mesmo resultado porque dentro do balde as três chaves juntas são únicas (a dedup
garante): a ordem total é a mesma e não há empate para desempatar. 20× mais
rápido.

### O que **não** foi feito, e por quê

**O Reanalyzing continua lá** — `Analyzed symmetries for 2 times`. A eliminação
dele estava no plano, e foi descartada depois de medir: com o hash, a passada
inteira da fase A custa ~0,7 s, então matar o retrabalho economiza ~0,7 s. Isso
é 1,4% do total, abaixo do ruído de medição, e exigiria prever o `nsymm` final
por simulação da regra de halving — a parte de maior risco de todo o plano, com
três modos de falha listados abaixo. **Não vale a troca.** Ficou registrado
porque volta a valer se o custo por passada crescer muito com o tamanho do
cluster.

Se um dia for feito, o roteiro é: (1) uma passada com hash sem capacidade fixa;
(2) achar o `nsymm` final simulando a regra original sobre as ~823 mil
assinaturas distintas em vez dos 15 milhões de trios — `nsymm = natoms²/20`,
`tscatter = mscatter/nsymm/3`, halving enquanto algum balde tiver contagem
`> tscatter`, mais a regra `if (nsymm<10 && mscatter<=5*natoms³)
mscatter=mscatter*3/2`; (3) montar o layout e remapear. Riscos, em ordem: a
ordem final tem de bater com a da compactação (linhas 493-511, onde o teste
`tempar[p*3]==0.0` faz dupla função de "consumida" e "comprimento zero"); o
`nsymm` simulado tem de bater com o do original; e a contagem por balde tem de
ser a de chaves `(l,v,c)` distintas, não a de entradas.

### Custo

+34 MB de tabela hash (dimensionada estritamente maior que o número de slots,
para que uma sondagem sem sucesso nunca entre em laço infinito). Pico de RSS
medido: **inalterado**, 582 MB — o pico do processo é fixado em outro ponto.

---

## Rascunho do plano original do V2 (mantido como registro)

Ataca A (6,75 s) e B (1,02 s) de uma vez. Depende do V1 estar fechado.

> **Este trecho contém um erro**, mantido de propósito porque a correção é o
> achado mais importante do V2: o conjunto distinto **depende** de `nsymm`,
> porque o balde é escolhido usando `akind`, que não faz parte da chave
> guardada. Ver a seção do V2 acima.

**A chave do problema.** O índice do balde `m` é **função pura da assinatura**
(linhas 401-404): assinaturas iguais caem sempre no mesmo balde, para qualquer
`nsymm`. Logo **o conjunto distinto não depende de `nsymm`**. Isso é o que
permite separar duas coisas que hoje estão amarradas: *descobrir quais são os
distintos* (caro, `natoms³`) e *escolher o `nsymm` e arrumar os baldes* (barato,
823 mil itens).

**O segundo fato que destrava.** Todas as três chaves são quantizadas antes de
comparar — `round(...,0.005f)` nas linhas 390/400/402 e `vlenc = kmin/100`. Chave
quantizada é inteiro exato, então dá para hashear com igualdade **idêntica** à
comparação de floats que o código faz hoje. O hash não é aproximação.

**O algoritmo.**

1. **Uma** passada sobre `natoms³` com hash de endereçamento aberto sobre
   `(lena, vlenb, cosbeta)` convertidos a inteiro. Sem capacidade fixa, logo sem
   estouro e **sem Reanalyzing**. Guarda o id da assinatura em `tevenadd` (como
   no V1). Custo por trio: ~1 sonda em vez de ~270 comparações.
2. **Descobrir o `nsymm` final** simulando a regra original sobre as 823 mil
   assinaturas distintas — não sobre os 15 milhões de trios:
   - começa em `nsymm = natoms²/20`;
   - `tscatter = mscatter/nsymm/3`;
   - conta quantos distintos caem em cada balde (recalcular `m` para 823 mil
     itens é barato);
   - se algum balde tem contagem `> tscatter`, então `nsymm/=2` (mais a regra de
     crescimento `if (nsymm<10 && mscatter<=5*natoms³) mscatter=mscatter*3/2`) e
     repete.

   Isso reproduz **exatamente** o `nsymm` a que o original converge. A
   equivalência vale porque o original dispara o estouro quando a contagem
   ultrapassa `tscatter` durante a inserção, e contagem final `> tscatter` ⟺
   ultrapassou em algum momento.
3. **Montar o layout** no `nsymm` final: ordenar cada balde por
   (cos β, 1/r₂, r₁) e concatenar na ordem `m = 1..nsymm-1`. Isso é exatamente o
   que as fases B produzem hoje — verificado relendo a compactação das linhas
   493-511, que emite as entradas na ordem ordenada. Trocar o *selection sort*
   O(n²) por `std::sort` com o mesmo comparador de 3 chaves.
4. **Remapear `tevenadd`**, como no V1.

**Onde isso pode dar errado, em ordem de risco.**

- A ordem final tem de bater com a do original. O ponto delicado é a
  compactação (linhas 493-511), que parece emitir na ordem ordenada mas tem o
  teste `tempar[p*3]==0.0` fazendo dupla função de "entrada consumida" e
  "comprimento zero". Conferir com um dump comparativo de `tevenpar` antes de
  confiar.
- Empates: dentro de um balde as três chaves juntas são únicas (a dedup
  garante), então `std::sort` não precisa ser estável.
- Se `nsymm` simulado não bater com o do original, o layout muda e a curva muda.
  Instrumentar e comparar `nsymm` e `ntrieven` antes de comparar a curva.

---

## Escalabilidade — V0 contra V2

Campanha fria de 04/08/2026, 20:47 (`baseline/campanha.sh`, `REPS=2 PAUSA=45`,
load 0,04 ao iniciar): `np` = 1, 2, 4, 6, 8, 12, duas versões, duas repetições,
**intercaladas dentro de cada `np`** para que deriva térmica atinja as duas
igualmente e a razão entre elas se preserve. Tabela pelo **mínimo** das
repetições (estimador padrão em benchmark: ruído só atrasa, nunca acelera).
Dados brutos em `baseline/escala.csv`, figura em `baseline/escalabilidade.png`
(`baseline/escala.py`).

**Estes números são da configuração corrigida** — `psAg111.txt`, k = 13,63 Å⁻¹,
`lnum=20` (ver "O k errado"). **Não são comparáveis** com os de qualquer versão
anterior deste arquivo, que foram medidos em k = 16,00 com `lnum=10`. O custo da
correção, medido no V2: 116,06 → 135,65 s em `np=1`, 46,24 → 49,84 s em `np=4`.

| `np` | V0 | V2 | ganho | speedup V0 | speedup V2 |
|---:|---:|---:|---:|---:|---:|
| 1 | 146,49 s | 135,65 s | 1,08× | 1,00× | 1,00× |
| 2 | 87,53 s | 73,91 s | 1,18× | 1,67× | 1,84× |
| 4 | **64,26 s** | **49,84 s** | 1,29× | **2,28×** | **2,72×** |
| 6 | 70,97 s | 53,28 s | 1,33× | 2,06× | 2,55× |
| 8 | 69,06 s | 51,89 s | 1,33× | 2,12× | 2,61× |
| 12 | 69,70 s | 49,86 s | **1,40×** | 2,10× | 2,72× |

**`factors = 0.6710 0.8642` nas 24 rodadas**, nas duas versões e em todos os
`np`. É a validação mais forte que temos: o resultado não depende nem da versão
nem do número de ranks.

Três leituras:

- **O ótimo é `np=4`, não `np=6`.** O `np=4` ganha do `np=6` nos **quatro
  pareamentos independentes** (V0 rep1 64,26 vs 76,77; V0 rep2 64,54 vs 70,97;
  V2 rep1 49,84 vs 53,28; V2 rep2 51,83 vs 55,81), com margem de 7% a 16%. O
  `np=12` empata com o `np=4` no mínimo (49,86 vs 49,84) mas com dispersão muito
  pior no V0 — 15,7% entre repetições contra 0,4% do `np=4`. **Use `-np 4`.**
- **O ganho do V2 cresce com o `np`**: 1,08× em `np=1` contra 1,40× em `np=12`.
  É o argumento central: quanto mais o laço paralelo encolhe, mais o preparo
  serial domina — e é exatamente ele que o V2 cortou. **Otimizar o serial é o que
  permite usar máquina grande.**
- **O gargalo mudou de lugar, e isto redireciona o próximo passo.** O preparo
  serial do V2 é 15,5 s em 135,65 s de `np=1`, ou seja 11,4% ⇒ teto de Amdahl
  **8,8×**. O medido é 2,72×. A diferença não está mais no serial: o laço
  paralelo sozinho vai de ~120 s (`np=1`) para 35,6 s (`np=6`, cronometrado no
  piloto), **3,4× com 6 ranks — 56% de eficiência**. Continuar cortando serial
  rende pouco agora; o alvo passou a ser por que o laço perde 44%.
  *(Derivado de um `np=6` instrumentado mais o total de `np=1`, não de uma
  medição dedicada por `np`. Confirmar antes de agir.)*

### Qualidade destes números

A dispersão entre repetições **cresce com o `np`**, o que faz sentido com 12
threads sobre 6 núcleos físicos:

| `np` | 1 | 2 | 4 | 6 | 8 | 12 |
|---|---:|---:|---:|---:|---:|---:|
| V0 | 1,0% | 2,7% | 0,4% | 8,2% | 8,8% | 15,7% |
| V2 | 0,4% | 1,8% | 4,0% | 4,7% | 7,0% | 1,6% |

Até `np=4` a reprodutibilidade é de fração de por cento — o `V0 np=1` deu 147,93
e 146,49 s. Compare com a campanha suja anterior, onde o mesmo ponto ia de 71,4
para 107,5 s: **o resfriamento inicial e as pausas de 45 s resolveram a deriva.**
Nos pontos de `np` alto a dispersão é real e a curva ali deve ser lida com essa
margem, não como valor exato.

### O `np=1` que não terminava — era bug, não contenção

Uma rodada de `np=1` do V0 passou **897 s sem terminar**, e este arquivo dizia
que era **contenção de CPU** (897/150 = 6,0×, "exatamente o que se espera de um
processo disputando 6 núcleos"). **Esse diagnóstico estava errado**, e a razão
por que sobreviveu é que a aritmética fechava bonito.

Em 04/08/2026 o travamento reapareceu, com a máquina **verificadamente vazia**
(`load 1.00` = só o processo, 5,6 GB livres, `VmSwap: 0`, `read_bytes: 832`).
Ficou 52 min em `R` a 100% de CPU sem escrever uma linha. Diagnóstico real, com
o `rip` colhido por `gdb -p` e conferido no *disassembly*: laço infinito em
**`phase.cpp:320`**, `while (xc-xb>90.0) xc-=180.0f`.

A cadeia, com a parte provada separada da inferida:

1. **Provado.** `phase.cpp:164` zera as 256×61 posições de `phasea`, então valor
   *de dentro* do array é da tabela (graus, ±180) ou zero — com ambos limitados o
   laço termina em ≤2 iterações, **sempre**. Logo o índice tinha de estar fora, e
   só existe um candidato: `(i-1)*61+j+1` com `i=0`, que dá −60..−51.
2. **Provado.** `phase.cpp:312-313` permite `i=0`: o guarda `if (i<0) i=0` é
   código morto, porque `i` nunca sai negativo daquele `for`. A intenção era
   `if (i<1) i=1`.
3. **Provado.** `rip = makephase+0x110`, um laço de 5 instruções
   (`subss`/`movaps`/`subss`/`comiss 90.0`/`ja`) — é a linha 320 compilada.
   Base PIE confere: `0x5ac4142fd5d0 − 0x2a5d0 = 0x5ac4142d3000`.
4. **Provado.** Com `xb` grande negativo, `xc` decrementa 180 por iteração até a
   precisão do `float` saturar (~3e9, onde `xc-180 == xc`) e **nunca termina**.
5. **Inferido.** `phasea` é `float*`; `phasea[-60]` cai ~240 bytes antes do bloco
   do `malloc`. Página nova do kernel vem zerada (caso comum ⇒ roda normal, e daí
   `factors` sair sempre igual); reuso de chunk sujo ⇒ valor grande ⇒ trava. A
   variação vem do Open MPI, que aloca em 3 threads com tempo indeterminado. Não
   foi possível ler a memória do processo travado sem root.

O gatilho era de configuração e está corrigido — ver "O k errado". Com
`psAg111.txt` o `i=0` deixa de ser alcançado no caso normal (13,63 fica com 8,6
Å⁻¹ de margem abaixo do início da tabela, em vez de zero), e o `np=1` rodou três
vezes na campanha fria sem repetir. **Mas o bug continua no código**: a correção
mínima é `if (i<1) i=1` na linha 313, e os quatro `while` de 320-325 seguem sem
limite de iteração.

Três erros de método que geraram isso, todos fáceis de repetir:

1. **Aceitar a explicação cuja aritmética fecha.** 897/150 = 6,0× com 6 núcleos é
   coincidência convincente, e foi ela que fechou a investigação cedo. Número
   redondo não é evidência causal.
2. **Não amostrar a pilha do processo vivo.** Custa 5 s (`gdb -p PID -batch -ex
   "bt 25"`), decide a questão, e a evidência morre com o processo. Precisa de
   `sudo` neste sistema: `ptrace_scope=1` só permite rastrear descendente.
3. **Tomar buffer de saída por progresso.** O `stdout` redirecionado é
   bufferizado em bloco, então o arquivo fica atrás do programa. Aqui isso
   *escondeu* o travamento de verdade: as "24 linhas paradas em 0,00%" foram
   descartadas como artefato de buffer quando eram sintoma real.

Vale para o `precutable` e para qualquer medição futura: **medir com a máquina
quieta, reproduzir antes de concluir, e amostrar a pilha antes de matar.**

## Perfil do `precutable` — o próximo alvo

Instrumentado em 6 pontos (`mscdrunc.cpp`, mesma macro `MSCDT`). Curva idêntica.

| trecho | tempo | |
|---|---:|---|
| coeficientes vibracionais | 0,001 s | |
| matrizes de rotação | 0,003 s | |
| expansão esférica + cortes | 0,058 s | |
| `alltrievent(1,kmin)` | 1,362 s | elementos de matriz de espalhamento |
| **laço do pathcut** | **8,615 s** | **82% do `precutable`** |
| estatísticas + alocações | 0,440 s | |
| **total** | **10,479 s** | |

O laço é `mscdrunc.cpp:404-462`: `for m=2..msorder { for ib { for ic { for ia`,
ou seja **7 × 247³ = 105 milhões de iterações**. Por (ib,ic) ele reduz sobre `ia`
— `xa` por máximo, `cxa` por soma — e no fim grava `bsum[ib*natoms+ic]`.

Três coisas visíveis na leitura, em ordem de ganho esperado:

1. **A ordem dos laços é a pior possível.** O índice é
   `id = ia*natoms² + ib*natoms + ic`, mas `ia` é o laço **mais interno**: cada
   passo de `ia` salta `natoms²` = 244 KB dentro de `tevenadd` (60 MB), e o mesmo
   vale para `tevendim[id]`. São 105 milhões de acessos praticamente aleatórios.
   Com `ia` no laço externo o acesso vira sequencial.

   **E a troca é segura**: `ia` é o índice da redução, `(ib,ic)` é o da saída.
   Basta manter acumuladores `xa[]`/`cxa[]` de tamanho `natoms²` (~0,7 MB) e
   adiar a escrita de `bsum`. A soma de `cxa` para cada `(ib,ic)` continua na
   mesma ordem de `ia`, então **o resultado em ponto flutuante não muda** — que
   é a condição para a curva continuar idêntica. `xa` é máximo, indiferente à
   ordem. E `bsum` não é lido dentro do laço (só `asum`, que é a cópia da
   iteração anterior de `m`), então adiar a escrita não cria dependência.

2. **`cabs(asum[ia*natoms+ib])` aparece 5 vezes com o mesmo argumento** (linhas
   421-425), uma raiz quadrada cada. O compilador provavelmente já faz a
   eliminação de subexpressão comum com `-O3`; medir antes de mexer.

3. **`xb`..`xf` são calculados sempre, mas usados numa cascata de `else if`.**
   Com `raorder=4` (o valor do `Cov0.txt`), o primeiro teste é sobre `xf`; se ele
   passa, os outros quatro produtos foram jogados fora. Calcular sob demanda
   economiza ~4 multiplicações por iteração, 105 milhões de vezes.

Nada disso é paralelismo — é a mesma família de correção algorítmica do V1/V2.

## V3 — inverter a ordem dos laços do pathcut: **piorou, foi revertido**

Resultado negativo, medido e mantido no registro porque a intuição que ele
derruba é forte e voltaria a aparecer.

**A hipótese.** O laço do pathcut (`mscdrunc.cpp:404-465`) é
`for m { for ib { for ic { for ia`, e o índice é `id = ia*natoms² + ib*natoms + ic`.
Com `ia` no laço **mais interno**, cada passo salta `natoms²` = 244 KB dentro de
`tevenadd` e `tevendim` (60 MB cada). Parecia o caso clássico de ordem de laço
errada: inverter para `(ia,ib,ic)` deixaria o acesso sequencial em `ic`.

**A transformação é válida.** `ia` é o índice da redução e `(ib,ic)` o da saída,
então basta manter acumuladores de tamanho `natoms²` e adiar a escrita de `bsum`.
Para cada `(ib,ic)` a soma de `cxa` continua percorrendo `ia` em ordem crescente,
logo a aritmética de ponto flutuante é a mesma; `xa` é máximo, indiferente à
ordem; e `bsum` não é lido dentro do laço (só `asum`, cópia do `m` anterior).
**A curva saiu idêntica byte a byte** — a transformação está correta.

**Mas ficou mais lento.** Comparando medições feitas próximas uma da outra, com
a máquina no mesmo estado — que é a única comparação defensável aqui:

| laço do pathcut | |
|---|---:|
| V3 (`ia` externo) | 10,930 s |
| revertido (`ia` interno) | 9,669 s |

**13% mais lento**, e o `precutable` inteiro de 11,55 s para 12,85 s.

⚠️ **Não compare com os 8,615 s da tabela do V2**: aquilo foi medido ~2 h antes,
com a máquina fria. Depois de carga contínua ela entrega ~15% menos — o mesmo
`-np 6` que dava 48,4 s passou a dar 56,5 s **sem mudança de código**. Os números
de fase só são comparáveis entre si dentro da mesma janela de medição.

**Por quê.** Na versão original `xa` e `cxa` vivem em **registrador** durante todo
o laço de `ia`. Com `ia` por fora eles viram acumuladores em memória, e cada uma
das 105 milhões de iterações passa a fazer load+store. O custo disso superou o
ganho de localidade — ou seja, **o laço não era dominado pelo passo de 244 KB**,
ao contrário do que eu supus. Os acessos que realmente pesam são provavelmente os
dependentes de dado (`tevenpar[k*10+…]`, `tevenelem[memadd+…]`, com
`k=tevenadd[id]`), que a inversão não muda.

**Revertido.** Ficou um comentário no código avisando para não repetir a
tentativa sem medir.

**A lição para o próximo passo:** parar de otimizar por hipótese de leitura. O
`precutable` precisa de **profiling de verdade** — contadores de cache miss e de
stall (`perf stat -e cache-misses,LLC-load-misses`) — para dizer onde o tempo vai
antes de qualquer nova transformação. Foi o segundo palpite meu que a medição
derrubou (o primeiro foi ter estimado a fase C igual à A, quando era o dobro).

## Depois do V2

Nesta ordem, pelo que as medições mostram:

1. **`precutable` (9,7 s)** é o maior item serial, com folga — sozinho é 74% do
   preparo. Perfil acima; o alvo é o laço do pathcut.
2. **Paralelizar o que sobrar.** As três fases são paralelizáveis: A por `ia`
   com hash local por thread e união no fim; B por balde (independentes), com
   prefix-sum para a compactação; C é perfeitamente paralela (cada trio escreve
   um elemento distinto de `tevenadd`), bastando resolver o "representante"
   `tevenpar[j*10+4,7,8,9]` pelo menor `(ia,ib,ic)`, que é o que a ordem serial
   já produz. **OpenMP, não MPI**: tudo já está no espaço do rank 0.
3. **GPU.** A fase A é 15 milhões de cálculos de assinatura independentes mais
   sondas de hash — caso de livro-texto. Fazer só depois de 1 e 2: o ganho
   algorítmico é maior que o de hardware, e é o mesmo trabalho necessário para
   portar.

Medido e descartado: `--mca mpi_yield_when_idle 1` **não** ajuda (69,87 s contra
70,42 s, dentro do ruído). Com `-np 6` numa máquina de 6 núcleos físicos / 12
threads os ranks ociosos cabem nos hyperthreads e não roubam do rank 0. A
observação do `README.md` sobre ranks ociosos vale para `-np 10`.

## Limite de escala achado no caminho

`tevenadd` é indexado por `ia*natoms*natoms+ib*natoms+ic` em `int`. Acima de
**natoms ≈ 1290** isso estoura o int de 32 bits e corrompe silenciosamente. É
pré-existente, não foi introduzido aqui, e ainda não foi corrigido.
