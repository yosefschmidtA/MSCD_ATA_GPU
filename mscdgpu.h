/*
----------------------------------------------------------------------
  Port de GPU do MSCD -- Fase 1: alldblevent / evenelem.
  Interface POD entre o C++98 de 1998 e o nvcc.

  Por que uma interface POD e nao os tipos do programa: fcomplex.h:46
  declara "friend Fcomplex polar(float,float=0)" -- argumento padrao num
  friend que nao e' definicao. E' ilegal no padrao, e' o unico motivo do
  -fpermissive no g++, e o front-end de device do nvcc recusa. Nenhum
  cabecalho do programa entra no .cu.
----------------------------------------------------------------------
*/

#ifndef __MSCDGPU_H
#define __MSCDGPU_H

/* Mesmo layout de Fcomplex (fcomplex.h:31) e de float2. */
typedef struct { float re,im; } Gcplx;

/* Tabelas constantes na corrida inteira: sobem para a placa uma vez.
   Valem porque kmin==kmax (a energia nunca muda), verificado no Cov0.txt. */
typedef struct
{ const float *patom;      int natoms;   /* natoms*12 */
  const float *devenpar;   int ndbleven; /* ndbleven*7 */
  const int *devenadd;     int msorder;
  const float *rotmata;    const float *rotmatc;
  int rlnum,lamdum,betanum;
  const Gcplx *hankmat_a;  int handata,halnum,hacmnum;
  const Gcplx *hankarg_b;                /* foto do cache de hankb */
  const Gcplx *phasec;     int pclnum;   /* por especie */
  int nkind;
  const Gcplx *cexpix;     int exndata,exmdata;
  const float *thermat;    int thernum;
  float therstep,mweight,tdebye,tsample;
  const float *aweight;
  int radim,raorder;   /* alnum vai por especie, em mscdgpu_set_alnum */
} Gconst;

#ifdef __cplusplus
extern "C" {
#endif

/* Sobe as tabelas constantes. Chamar uma vez. 0 = ok. */
int mscdgpu_setup(const Gconst *k);

/* Um ponto: calcula devenelem[ndbleven*radim] para este xdetec.
   xc = meanpath->finvpath(akin), constante, calculado no host. */
int mscdgpu_alldblevent(float akin,const float *xdetec,float xc);
int mscdgpu_get_devenelem(Gcplx *out);
int mscdgpu_allevendetec(float akin,const float *xdetec,float xc,Gcplx *devendetec_out);

int mscdgpu_setup_summation(
    const int *tevencut, const int *tevendim, const int *tevenadd, 
    const float *tevenpar, const float *talpha, const float *tgamma,
    int ntrieven, int ntrielem, const float *patom, int msorder);

int mscdgpu_summation(float akin, const Gcplx *tevenelem, Gcplx *asum_host, const float *patom);

void mscdgpu_teardown(void);
const char *mscdgpu_lasterror(void);

#ifdef __cplusplus
}
#endif

#endif
