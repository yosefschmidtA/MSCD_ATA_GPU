# MSCDATA — mapa do projeto

**MSCD 1.37** (Van Hove / LBNL, 1997-98): difração e dicroísmo de fotoelétrons.
C++ antigo, 69 arquivos, ~19 mil linhas, paralelo por MPI. Não é o MSCD de
fábrica: tem uma extensão **ATA** da UNICAMP (`Mscdrun::ATAevenelem`,
`mscdrunc.cpp:45`), desligada no `Cov0.txt` (`ATA=0`).

Repositório: <https://github.com/yosefschmidtA/MSCD_ATA_GPU>. **O git é do
usuário** — não rode `add`, `commit`, `tag` nem `push`; entregue o comando pronto.

**Antes de abrir qualquer `.cpp`, leia estes dois** — eles existem justamente
para economizar leitura de código:

- **`README.md`** — caderno de laboratório. A fatoração caro/barato da eq. (46)
  que organiza o programa inteiro, o grafo das equações, onde está (e onde não
  está) o paralelismo, e as medições.
- **`EQUACOES.md`** — as ~40 equações do manual (`MANUAL_MSCD_CEA.pdf`) mapeadas
  para arquivo:linha, nos dois sentidos.

O único executável usado é o **`randmscd_parallel`**. O `makefile` constrói outros
doze; estão compilados e funcionando, e não interessam.

## Compilar e rodar

```bash
make randmscd_parallel CPPFLAGS="-O3 -std=c++98 -w -fpermissive"
mpirun --use-hwthread-cpus -np 6 randmscd_parallel Cov0.txt
```

As duas flags são obrigatórias, não preferência:

- **`-fpermissive`** — `fcomplex.h:46` declara `friend Fcomplex polar(float,float=0)`.
  Argumento padrão em `friend` que não é definição: ilegal no padrão, aceito pelos
  compiladores de 1998, recusado pelo g++ 13. É o **único** erro em 69 arquivos.
- **`--use-hwthread-cpus`** — o Open MPI conta núcleos **físicos** (6 num
  i5-13420H), não threads (12), e desde a série 3.x recusa superalocar em vez de
  aceitar calado como a 1.10 fazia. Sem ela, `-np 10` morre em "not enough slots".

Referência de saída correta: `factors = 0.6724 0.8647`, curva em
`saida1Co-alterado-alexandre.txt`.

**Use `-np 6`, não `-np 10`.** Medido em 04/08/2026 com build instrumentado e
saída redirecionada para arquivo (imprimir no terminal custa ~10 s):

| `-np` | preparo serial | laço dos 779 pontos | total |
|------:|---------------:|--------------------:|------:|
| 1 | 15,6 s | 84,2 s | 99,8 s |
| 2 | 14,8 s | 41,9 s | 56,9 s |
| 6 | 16,1 s | 25,0 s | **41,4 s** |
| 10 | 22,0 s | 21,9 s | 44,0 s |

`-np 10` **perde** para `-np 6`: os 9 ranks ociosos giram em espera ocupada
(92–98% de CPU com 15 MB residentes, antes de receberem o job de 190 MB) e roubam
ciclos do rank 0 durante o preparo serial. Fração serial 15,6% ⇒ teto de Amdahl
6,4×. Detalhes e método no `README.md`.

*(Uma versão anterior deste arquivo dizia `-np 10` = 41,4 s e `-np 6` = 45,2 s.
Está errado — foi medido com saída no terminal e sem separar as fases.)*

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
