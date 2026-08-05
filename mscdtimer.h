#ifndef MSCDTIMER_H
#define MSCDTIMER_H

/* Cronometro de fases, so' existe quando se compila com -DMSCDTIMER.
   Sem a macro nao sobra instrucao nenhuma no binario -- e' por isso que a
   instrumentacao pode ficar no fonte sem contaminar a producao.
   Escreve em stderr de proposito: stdout e o flogout vao para arquivos que
   entram na comparacao de regressao. */

#ifdef MSCDTIMER

#include <time.h>
#include <stdio.h>

static double mscdt_now()
{ struct timespec t;
  clock_gettime(CLOCK_MONOTONIC,&t);
  return (double)t.tv_sec+1.0e-9*(double)t.tv_nsec;
}

#define MSCDT_DECL double mscdt_p=mscdt_now(), mscdt_c
#define MSCDT(nome) do { mscdt_c=mscdt_now(); \
    fprintf(stderr,"[timer] %-26s %8.3f s\n",(nome),mscdt_c-mscdt_p); \
    fflush(stderr); mscdt_p=mscdt_c; } while (0)

/* Acumuladores, para regiao chamada muitas vezes: summation roda 3895 vezes,
   e um fprintf por chamada custaria mais que a regiao medida. Acumula em
   memoria e imprime uma vez no fim.
   Os acumuladores sao static: cada unidade de traducao tem os seus. Marcar e
   relatar tem de ficar no mesmo .cpp -- hoje e' o mscdrund.cpp. */

#define MSCDT_NACC 12

static double mscdt_acc[MSCDT_NACC]={0};
static long mscdt_cnt[MSCDT_NACC]={0};
static const char *mscdt_nam[MSCDT_NACC]={0};

#define MSCDT_ADECL double mscdt_a=mscdt_now(), mscdt_b
#define MSCDT_AMARK do { mscdt_a=mscdt_now(); } while (0)
#define MSCDT_A(slot,nome) do { mscdt_b=mscdt_now(); \
    mscdt_acc[slot]+=mscdt_b-mscdt_a; mscdt_cnt[slot]+=1; \
    mscdt_nam[slot]=(nome); mscdt_a=mscdt_b; } while (0)

#define MSCDT_AREPORT(rotulo,mype) do { int _i; double _t=0.0; \
    for (_i=0;_i<MSCDT_NACC;++_i) if (mscdt_nam[_i]) _t+=mscdt_acc[_i]; \
    fprintf(stderr,"[acc rank %d] %s  soma dos trechos %8.3f s\n", \
      (mype),(rotulo),_t); \
    for (_i=0;_i<MSCDT_NACC;++_i) if (mscdt_nam[_i]) \
      fprintf(stderr,"[acc rank %d]   %-30s %8.3f s  %6.2f%%  %9ld chamadas\n", \
        (mype),mscdt_nam[_i],mscdt_acc[_i], \
        (_t>0.0?100.0*mscdt_acc[_i]/_t:0.0),mscdt_cnt[_i]); \
    fflush(stderr); } while (0)

#else

#define MSCDT_DECL
#define MSCDT(nome)
#define MSCDT_ADECL
#define MSCDT_AMARK
#define MSCDT_A(slot,nome)
#define MSCDT_AREPORT(rotulo,mype)

#endif
#endif
