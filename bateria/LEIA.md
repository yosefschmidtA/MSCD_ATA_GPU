# Bateria de desempenho, original contra GPU

Mede o executável original (o do primeiro commit, só com o limite de átomos
aumentado de 300 para 1250) e a versão com GPU no `Cov0.txt`, com raio e
profundidade crescentes e `-np` de 1 a 12.

```bash
./bateria/bateria.sh             # roda, ou retoma de onde parou
./bateria/bateria.sh resumo      # mostra a tabela do que já foi medido
touch bateria/PAUSAR             # pausa ao fim da rodada atual
```

Ctrl+C para na hora, e a rodada que estava em andamento é refeita na próxima
vez. Pode desligar o computador a qualquer momento, porque cada rodada
terminada já está gravada em `bateria/resultados.csv`.

Os dois executáveis ficam congelados em `bateria/bin/` (md5 em
`bin/md5.txt`), então recompilar o projeto no meio não muda o que é medido.

## A grade

| raio / profundidade | átomos |
|---|---:|
| 8 / 12 | 135 |
| 10 / 15 | 253 |
| 12 / 18 | 424 |
| 14 / 21 | 643 |
| 16 / 24 | 932 |
| 17 / 26 | 1112 |

## O que o script pula sozinho, e por quê

- **Memória.** O original guarda quatro tabelas `natoms³` e cada processo
  recebe uma cópia, então a memória cresce com o número de processos. Antes
  de cada rodada o script estima o que ela vai usar e pula se não couber
  (status `PULADO_MEMORIA`). O original com 932 átomos já não cabe nos 15 GB
  nem com um processo.
- **Tempo.** O Reanalyzing do original pode não terminar. Cada rodada tem no
  máximo `TEMPO_MAX` segundos (padrão 10 horas). Na primeira que estoura
  (`TEMPO_ESGOTADO`), o resto daquele tamanho e todos os maiores daquela
  versão são pulados. O que trava é o preparo, que roda só no processo 0, então
  outro `-np` do mesmo tamanho não terminaria também.
- **VRAM.** Com 1112 átomos o `pathcut` não cabe nos 8 GB da placa, e o script
  liga `MSCD_PATHCUTCPU=1` sozinho. O resultado é o mesmo.

## Ajustes

Tudo por variável de ambiente, por exemplo

```bash
VERSOES=gpu ./bateria/bateria.sh
NPS="1 2 4 8 12" TEMPO_MAX=3600 ./bateria/bateria.sh
TAMANHOS="8:12 10:15" ./bateria/bateria.sh
```

Para medidas limpas, deixe a máquina parada durante a bateria e o protetor de
tela desligado (`omarchy toggle screensaver`).

## Em outra máquina (Ubuntu com placa NVIDIA)

Os executáveis de `bin/` não vão para o git, porque um binário compilado aqui
não roda em outra distribuição. Depois do `git clone` (sem `--depth`, o
script precisa do commit `d4c8408`)

```bash
./bateria/instalar_ubuntu.sh
```

instala compilador, Open MPI, driver e CUDA, compila os dois executáveis a
partir do git e roda 135 átomos comparando com `saidas/gpu_r8p12_np1.out`
deste PC. Se precisar instalar o driver, ele pede para reiniciar e rodar de
novo. No fim imprime o comando da bateria, que é

```bash
MAQUINA=nome VERSOES=gpu ./bateria/bateria.sh
```

Com `MAQUINA` os resultados vão para `resultados_nome.csv` e `saidas_nome/`,
sem misturar com os deste PC. Use nome curto, porque o MSCD corta caminho de
arquivo em 100 caracteres.

Conferido em 09/10/2026 neste PC. O original reconstruído do git tem o mesmo
md5 do `bin/randmscd_original`, e o de GPU dá a mesma curva byte a byte.

### Grade de centenas (`TAMANHOS=centenas`)

Raios escolhidos para dar perto de 100, 200, ... 1100 átomos, com a
profundidade em 1,5 vez o raio. Contados com o original em 09/10/2026.

| raio / prof | átomos | | raio / prof | átomos |
|---|---:|---|---|---:|
| 7 / 10,5 | 96 | | 14,25 / 21,375 | 686 |
| 9,25 / 13,875 | 210 | | 15,25 / 22,875 | 785 |
| 10,75 / 16,125 | 304 | | 15,75 / 23,625 | 890 |
| 11,75 / 17,625 | 394 | | 16,5 / 24,75 | 1016 |
| 13 / 19,5 | 508 | | 17 / 26 | 1112 |
| 13,75 / 20,625 | 592 | | | |

```bash
MAQUINA=nome VERSOES=gpu TAMANHOS=centenas ./bateria/bateria.sh
```

## Com mais de um processo, o limite é ~500 átomos (ERRO_134)

Com `-np` maior que 1, o processo 0 empacota o job inteiro numa mensagem só
para mandar aos outros (`Mscdjob::sendjobs`, `mscdjob.cpp:108`). O tamanho
dela é somado num `int` em `Mscdrun::getlength` (`mscdrun.cpp:168`), e as
tabelas de trios entram com 16 bytes por trio. Com 508 átomos isso passa de
2³¹ bytes, o número fica negativo, o `new` lança `std::bad_alloc` e o
programa aborta com código 134. Reproduzido em 09/10/2026 com 508 átomos e
`-np 2`, com 11 GB livres para uma rodada de ~7 GB, então não é falta de RAM.

Vale para o original e para a GPU, porque o código é o mesmo. Com 424 átomos
ainda passa, com 508 já não. Os `ERRO_134` da bateria a partir de ~500 átomos
com `np>1` são esse limite, e o `np=1` não é afetado. Para passar dele seria
preciso trocar o `int` por 64 bits na serialização e dividir o `MPI_Send` em
pedaços menores que 2 GB.
