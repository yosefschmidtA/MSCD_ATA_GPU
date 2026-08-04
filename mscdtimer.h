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

#else

#define MSCDT_DECL
#define MSCDT(nome)

#endif
#endif
