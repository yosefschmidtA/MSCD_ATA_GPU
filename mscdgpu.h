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
/* As duas metades do summation: launch enfileira tudo na placa e volta na
   hora; finish espera e copia as linhas dos emissores para asum_host. */
int mscdgpu_summation_launch(float akin, const Gcplx *tevenelem);
int mscdgpu_summation_finish(Gcplx *asum_host);
/* Chamadas encadeadas: launch2 devolve o slot, finish2 espera so' ele. */
int mscdgpu_summation_launch2(float akin, const Gcplx *tevenelem, int *slot);
int mscdgpu_summation_finish2(int slot, Gcplx *asum_host);

/* pathcut do precutable na placa (08/10/2026). Thread por (ib,ic), laco de
   ia em serie, mesma ordem da CPU. O pow do passo m=2 fica no host: a libm
   da CUDA nao garante os bits da glibc. bsum tem natoms^2 complexos. */
int mscdgpu_pathcut_begin(int natoms,int msorder,int raorder,float pathcut,
  const float *patom,const int *tevenadd,const float *tevenpar,int ntrieven,
  const Gcplx *tevenelem,int ntrielem,const Gcplx *bsum);
int mscdgpu_pathcut_step(int m,float *xa_m2);
int mscdgpu_pathcut_setbsum(const Gcplx *bsum);
int mscdgpu_pathcut_end(int *tevencut,int *tevendim);
/* Fase 6: rowsonly!=0 desce do tevendim so' as linhas (ia,ib) com algum
   tevencut, as unicas que o host le no modo GPU com np=1. */
int mscdgpu_pathcut_end2(int *tevencut,int *tevendim,int rowsonly);

/* Fase 6: sem talpha/tgamma, o setup pede a rotacao de cada trio que
   precisa dela (evedim>1) por esta funcao. */
typedef void (*mscdgpu_rotfn)(int ia,int ib,int ic,float *alpha,float *gamma);
void mscdgpu_set_rotfn(mscdgpu_rotfn f);

void mscdgpu_teardown(void);
/* So' para o build -DMSCDTIMER: espera a placa, para os acumuladores
   medirem o kernel e nao o lancamento assincrono. */
void mscdgpu_sync(void);
const char *mscdgpu_lasterror(void);

#ifdef __cplusplus
}
#endif

#endif
