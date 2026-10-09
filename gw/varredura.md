# Varredura da Fase 6

Lista de todos os fontes do randmscd_gpu. Marcar cada um quando lido inteiro.
Caso de teste principal gw/c963.in (963 atomos). Referencias bit a bit em
gw/cov_refgpu.out, gw/iron_refgpu.out e gw/c963_old.out.

Perfil de partida (963 atomos, randmscd_gpu_timer, depois do symdblvert com
hash): total 70 s. symtrivert 18,3 s, estatisticas do precutable 18,2 s,
allrotation 10,1 s, pathcut 6,2 s, alltrievent 1,7 s, esfericas 1,3 s,
laco dos pontos 22,6 s (build instrumentado, sem encadeamento).

## Feito antes da lista
- symdblvert (mscdrunb_not_reanalize.cpp). Busca linear O(pares x distintos)
  trocada por hash com a mesma ordem de insercao. 963 atomos 237 s -> 0,04 s.
  Chave MSCD_DBLSERIAL=1.

## Etapas feitas (963 atomos, producao, np=1)
| etapa | onde | chave | 963 atomos |
|---|---|---|---:|
| partida | | | 310 s |
| symdblvert com hash | mscdrunb_not_reanalize.cpp | MSCD_DBLSERIAL | 70 s (perfil) |
| estatisticas do precutable em paralelo | mscdrunc.cpp | MSCD_STATSERIAL | 62,9 s |
| allrotation em paralelo | mscdrunc.cpp | MSCD_ROTSERIAL | 63,7 s (swap) |
| sem talpha/tgamma (rotacao sob demanda) | mscdrunc.cpp, mscdgpu.cu | MSCD_ROTFULL | 51,0 s, pico 13,2 -> 8,2 GB |
| dedup por tabela direta | mscdrunb_not_reanalize.cpp | MSCD_DEDUPHASH | 38,3 s |
| estatisticas numa varredura contigua | mscdrunc.cpp | MSCD_STATSERIAL | 28,8 s |
| pathcut com registro por trio | mscdgpu.cu | (sem chave, mesmo kernel) | 25,5 s |
| dedup com o tevenadd guardando a entrada (uma varredura a menos) + paginas grandes | mscdrunb_not_reanalize.cpp | MSCD_NOHUGE | ~26 s (ruido) |
| allevendetec e alldblevent so' no conjunto U (0,9% das posicoes) | mscdgpu.cu | MSCD_FULLMAT | 19,8 s (ferro 10,1 s) |
| euler pre-calculado por chamada, somas finais inline | mscdrunc.cpp, mscdrund.cpp | MSCD_EULERFULL | 18,3 s |
| fexpix e produtos complexos inline nos valores do bloco final | mscdrunc.cpp | (sem chave) | 18,2 s (sem ganho, mantido) |
| round() inline na dedup, sem zerar tevendim no host | mscdrunb_not_reanalize.cpp, mscdrunc.cpp | MSCD_ZEROALL | 16,9 s |
| descer do tevendim so' as linhas usadas | mscdgpu.cu, mscdrunc.cpp | MSCD_DIMFULL | 15,7 s |

## Arquivos

Dois metodos. "Lido" quer dizer as funcoes que rodam no preparo ou no laco
lidas inteiras. "Mapa + perfil" quer dizer que gw/laços.py listou todos os
lacos aninhados do arquivo e o perfil de 963 atomos mostrou que o arquivo
inteiro cabe no tempo nao medido (1,1 s de 15,5 s, dos quais 0,18 s de
pos-processamento e o resto partida do MPI e liberacao de memoria).

| arquivo | metodo | o que foi feito ou descartado |
|---|---|---|
| mscdrunc.cpp | lido | estatisticas, allrotation, talpha/tgamma, euler, inline no bloco final e nas cadeias do alltrievent, zeragem do tevendim, descida parcial. precutable "maximo do pemeven" (0,14 s) mantido. |
| mscdrunb_not_reanalize.cpp | lido | symdblvert com hash, dedup por tabela direta com tevenadd guardando a entrada, round inline, paginas grandes. |
| mscdrund.cpp | lido | somas finais inline. O laco dos pontos e o summation ja tinham sido refeitos na Fase 5. |
| mscdgpu.cu | lido | conjunto U, registro do pathcut, descida parcial do tevendim. k_sum_step e k_sum_prod antigos ficaram no arquivo sem uso (avisos do nvcc). |
| rotamat.cpp | lido (makerotation, rotelem, rotharma, termination) | nada novo. makecurve roda uma vez. |
| msfuncs.cpp | lido (Hankel, Expix) | Expix::fexpix copiado inline no bloco final. makecurve roda uma vez. |
| cartesia.cpp | lido | euler recalculava seno e cosseno a cada chamada, resolvido por pre-calculo no chamador. |
| fcomplex.cpp | lido | operadores fora de linha. Copiados inline nos tres lacos quentes, sem mexer no arquivo. |
| userutil.cpp | lido (round, confine) e mapa | round copiado inline na dedup. O resto e' Textout, saida de texto. |
| radmat.cpp | lido (famphase) e mapa | famphase tem busca linear numa tabela de dezenas de energias, 2 vezes por bloco de emissor. Desprezivel. |
| meanpath.cpp | lido (finvpath) e mapa | cache por energia, energia fixa. Nada. |
| vibrate.cpp | lido (fvibmsrd) e mapa | interpolacao pura. Nada. |
| phase.cpp | mapa + perfil, fsinexpa lido | makephase roda uma vez por energia. O bug do phase.cpp:312 (CLAUDE.md, Armadilhas) continua sem gatilho e sem correcao, porque corrigir muda a fisica no caso que hoje nao ocorre. |
| mscdruna.cpp | lido (makeatoms) e mapa | busca linear com tolerancia 0,1 A sobre <=1250 atomos, milissegundos e dependente de ordem. Mantida. readparameter roda uma vez (0,006 s ate o fim da leitura). |
| mscdrune.cpp | mapa + perfil | fitphotoemission e savecurve, uma vez por job. |
| mscdrun.cpp | mapa + perfil | exportacao para MPI (so' np>1) e assistant. |
| mscdjob.cpp | mapa + perfil | leitura da lista de jobs. |
| mscdmain.cpp | mapa | sem lacos aninhados. |
| pdinten.cpp | mapa + perfil | leitura da intensidade experimental, uma vez. |
| pdintena.cpp | mapa + perfil | chicalc e savecurve depois do laco, 0,18 s no total. |
| pdchifit.cpp | mapa | so' com nfit>0 (ajuste). Nao roda nos casos de teste. |
| curvefit.cpp | mapa | so' com nfit>0. Nao roda nos casos de teste. |
| polation.cpp | mapa | splines usadas na leitura. |
| userinfo.cpp | mapa | leitura de parametros. |
| jobtime.cpp | mapa | sem lacos aninhados. |
| userCluster.cpp | mapa | camada de MPI, sem lacos. |

## O que sobra (963 atomos, 15,5 s)
dedup 4,1 s (calculo da assinatura de 893 milhoes de trios, limitado por
conta), pathcut 2,4 s (subida de 3,6 GB e leituras aleatorias no kernel),
alltrievent 1,3 a 1,6 s, laco 5,1 s (bloco final na CPU com geometria em
double). Nenhum desses e' busca linear.
