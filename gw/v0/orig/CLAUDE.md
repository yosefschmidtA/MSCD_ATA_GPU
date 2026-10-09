# MSCDATA — mapa do projeto

**MSCD 1.37** (Van Hove / LBNL, 1997-98): difração e dicroísmo de fotoelétrons.
C++ antigo, 69 arquivos, ~19 mil linhas, paralelo por MPI. **Não é repo git** —
faça backup antes de mexer.

O único executável usado é o **`randmscd_parallel`**. O `makefile` constrói outros
doze; estão compilados e funcionando, e não interessam.

## Compilar e rodar

```bash
make randmscd_parallel CPPFLAGS="-O3 -std=c++98 -w -fpermissive"
mpirun --use-hwthread-cpus -np 10 randmscd_parallel Cov0.txt
```

As duas flags são obrigatórias, não preferência:

- **`-fpermissive`** — `fcomplex.h:46` declara `friend Fcomplex polar(float,float=0)`.
  Argumento padrão em `friend` que não é definição: ilegal no padrão, aceito pelos
  compiladores de 1998, recusado pelo g++ 13. É o **único** erro em 69 arquivos.
- **`--use-hwthread-cpus`** — o Open MPI conta núcleos **físicos** (6 num
  i5-13420H), não threads (12), e desde a série 3.x recusa superalocar em vez de
  aceitar calado como a 1.10 fazia. Sem ela, `-np 10` morre em "not enough slots".

Medido em 04/08/2026: `-np 10` = **41,4 s**, `-np 6` = 45,2 s. Os 4 ranks extras
rendem 8% — escala quase linear até 6 e satura. Referência de saída correta:
`factors = 0.6724 0.8647`, resultado em `saida1Co-alterado-alexandre.txt`.

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

## Backups de 04/08/2026 (nada apagado)

`backup-original-20260804.tar.gz` (binários e `.o` originais, `tar xzf` restaura)
e `.obj-antigos/` (os 40 `.o` da máquina antiga).
