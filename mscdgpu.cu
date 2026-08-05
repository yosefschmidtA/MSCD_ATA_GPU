/*
----------------------------------------------------------------------
  Port de GPU do MSCD -- Fase 1: alldblevent / evenelem (57% do laco).
  Ver PLANO_CUDA.md. Interface em mscdgpu.h.

  A ARITMETICA E' COPIA FIEL DA CPU, na mesma ordem de operacoes.
  Nao e' preciosismo: enquanto o kernel bate bit a bit com o
  mscdrunc.cpp, qualquer divergencia nas fases seguintes tem uma causa
  so'. Otimizar a ordem (fatorar cxb*cxc para fora do laco de ma, por
  exemplo) economiza ~metade das multiplicacoes complexas e pode ser
  feito depois -- mas ai a validacao vira "dentro de 1e-4" e para de
  apontar o dedo para o culpado.

  Onde a CPU calcula em double e arredonda para float (sqrt, atan2,
  acos, cos, sin, exp -- todos vem de <math.h> com argumento promovido),
  o device tambem usa double. FP64 e' 1/64 nesta placa, mas sao ~6
  chamadas por par contra ~5400 flops float: medido irrelevante, e e' o
  que mantem o beta discreto igual ao da CPU. Ver "O risco do beta".
----------------------------------------------------------------------
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include "mscdgpu.h"

/* Escrito exatamente como nas funcoes da CPU. NAO substituir por um literal
   float: radian divide ka, que chega a ~14000 graus, e o resultado vira
   INDICE de cexpix em fexpix -- 1 ulp aqui pula uma entrada da tabela. Os
   dois compiladores fazem 3.14159265/180.0 no mesmo double IEEE e arredondam
   para float do mesmo jeito, entao a constante e' bit a bit a mesma. */
#define RADIANF ((float)(3.14159265/180.0))

static char g_err[512]={0};
static int  g_ready=0;
static int  g_dbg=0;   /* so' o modo validate liga: sao 650 KB por ponto */

#define CK(call) do { cudaError_t e_=(call); if (e_!=cudaSuccess) { \
  snprintf(g_err,sizeof(g_err),"%s:%d %s: %s",__FILE__,__LINE__,#call, \
    cudaGetErrorString(e_)); return 1; } } while (0)

/* ---- aritmetica complexa, mesma semantica de fcomplex.h ---- */
__host__ __device__ static inline Gcplx cmul(Gcplx a,Gcplx b)
{ Gcplx r; r.re=a.re*b.re-a.im*b.im; r.im=a.re*b.im+a.im*b.re; return r; }
__host__ __device__ static inline Gcplx csmul(float s,Gcplx a)
{ Gcplx r; r.re=s*a.re; r.im=s*a.im; return r; }
__host__ __device__ static inline Gcplx cadd(Gcplx a,Gcplx b)
{ Gcplx r; r.re=a.re+b.re; r.im=a.im+b.im; return r; }
__host__ __device__ static inline Gcplx cneg(Gcplx a)
{ Gcplx r; r.re=-a.re; r.im=-a.im; return r; }

/* ---- estado do device ---- */
typedef struct
{ float *patom,*devenpar,*rotmata,*rotmatc,*thermat,*aweight;
  Gcplx *hankmat_a,*hankarg_b,*phasec,*cexpix;
  float *pairgeo;          /* 3 floats por par: cosa,sina,phia */
  int   *pairkind;         /* akind-1 */
  Gcplx *devenelem;
  int   *alnum;            /* por especie */
  float *dbg;
} Dev;

static Dev D;
static Gconst K;

/* Os cinco ma usados por alldblevent sao 0..4 e, com mb=0,
   getkelem(ma,0)==getkharm(ma,0)==ma*(ma+1) -- sempre >=0, entao o ramo
   de sinal de rotharma (kelem<0) nunca dispara aqui. */
__constant__ int c_kel[5]={0,2,6,12,20};

/* As 15 saidas vem de 9 evenelem distintos: (ma,na) e qual multiplicador
   complexo (cxa..cxi) cada saida usa. Tabela em mscdrunc.cpp:815-895. */
__constant__ int c_ma[9] ={0,1,0,2,1,3,0,2,4};
__constant__ int c_na[9] ={0,0,1,0,1,0,2,1,0};

/* ---------------- Expix::fexpix (msfuncs.cpp:306) ---------------- */
__device__ static inline Gcplx d_fexpix(const Gcplx *cexpix,int ndata,
  int mdata,float xa)
{ int k;
  /* A reducao e' por subtracao repetida, nao fmodf: ka/radian chega a
     ~14000 graus e as ~39 subtracoes em float nao dao o mesmo resultado
     que um resto exato. Trocar isto muda o indice k. */
  while (xa>180.0) xa-=360.0f;
  while (xa<-180.0) xa+=360.0f;
  k=(int)(mdata+(ndata-1.0)*xa/360.0+0.5);
  if (k<0) k=0; else if (k>ndata-1) k=ndata-1;
  return cexpix[k];
}

/* ------------- Vibration::fvibmsrd (vibrate.cpp:149) ------------- */
__device__ static inline float d_fvibmsrd(const float *thermat,int thernum,
  float therstep,float mweight,float tdebye,float tsample,
  float bondlength,float aweight)
{ int i; float y,value;
  if (bondlength<1.0/therstep/thernum) y=thernum-2.0f;
  else y=1.0f/bondlength/therstep;
  if (y<0.0) y=0.0f; else if (y>thernum-2.0) y=thernum-2.0f;
  i=(int)y;
  if ((mweight>0.1)&&(aweight>0.1)&&(tdebye>0.1)&&(tsample>0.1))
  { value=thermat[i]+(y-(float)i)*(thermat[i+1]-thermat[i]);
    value=value*(float)sqrt((double)(mweight/aweight));
  }
  else value=0.0f;
  return value;
}

/* ---- precomputo da perna (ia->ib), constante em todos os pontos ----
   Em alldblevent o terceiro ponto e' patom[ib]+xdetec, entao o vetor
   (patomc-patomb) E' xdetec: igual para todos os pares num ponto. So' a
   perna (patomb-patoma) varia por par, e ela nao depende do detector.
   O onerotation da CPU refaz este sqrt/atan2 779 vezes por par. */
__global__ static void k_pairgeo(const float *patom,const float *devenpar,
  int ndbleven,float *pairgeo,int *pairkind)
{ int j=blockIdx.x*blockDim.x+threadIdx.x;
  if (j>=ndbleven) return;
  int ia=(int)devenpar[j*7+5], ib=(int)devenpar[j*7+6];
  float xa=patom[ib*12]-patom[ia*12];
  float ya=patom[ib*12+1]-patom[ia*12+1];
  float za=patom[ib*12+2]-patom[ia*12+2];
  float lenga=(float)sqrt((double)(xa*xa+ya*ya+za*za));
  float cosa=za/lenga, sina, phia;
  if (fabs((double)xa)+fabs((double)ya)<1.0e-5) phia=0.0f;
  else phia=(float)atan2((double)ya,(double)xa);
  /* mscdrunc.cpp:694 e' sqrt(1.0-cosa*cosa): cosa*cosa e' FLOAT, so' a
     subtracao promove. Calcular o quadrado em double muda zc no 8o digito,
     e como beta=acos(zc) vira indice inteiro de rotmata -- com derivada
     1/sqrt(1-zc^2), que explode perto de zc=+-1 -- isso vira erro de 1 grau.
     Foi exatamente este bug na primeira versao: maxrel 7,25. */
  if (fabs((double)cosa)>=1.0) sina=0.0f;
  else sina=(float)sqrt(1.0-(double)(cosa*cosa));
  pairgeo[j*3]=cosa; pairgeo[j*3+1]=sina; pairgeo[j*3+2]=phia;
  pairkind[j]=(int)devenpar[j*7+4]-1;
}

/* --------------------------- o kernel --------------------------- */
struct Kargs
{ const float *devenpar,*pairgeo,*rotmata,*rotmatc,*thermat,*aweight;
  const Gcplx *hankmat_a,*hankarg_b,*phasec,*cexpix;
  const int *pairkind,*alnum;
  Gcplx *out;
  float *dbg;              /* 4 por par: beta,gamma,xb,row -- so' depuracao */
  int ndbleven,rlnum,lamdum,betanum,handata,halnum,hacmnum,pclnum;
  int exndata,exmdata,thernum,radim,raorder;
  float therstep,mweight,tdebye,tsample;
  float akin,xc,cosb,sinb,phib;
};

__global__ static void k_alldblevent(Kargs a)
{ int j=blockIdx.x*blockDim.x+threadIdx.x;
  if (j>=a.ndbleven) return;

  /* --- onerotation (mscdrunc.cpp:672), so' beta e gamma: alpha nao e'
     consumido por alldblevent, entao os dois atan2 dele saem. --- */
  float cosa=a.pairgeo[j*3], sina=a.pairgeo[j*3+1], phia=a.pairgeo[j*3+2];
  float cosb=a.cosb, sinb=a.sinb;
  float phiab=phia-a.phib;
  float cosphi=(float)cos((double)phiab), sinphi=(float)sin((double)phiab);
  float xd_=sinb*cosa*cosphi-cosb*sina, yd_=sinb*sinphi;
  float zc=cosb*cosa+sinb*sina*cosphi;
  float beta,gamma;
  if ((zc>0.99999)&&(cosb>=0.0))      { beta=0.0f;   gamma=0.0f; }
  else if (zc>0.99999)                { beta=0.0f;   gamma=0.0f; }
  else if ((zc<-0.99999)&&(cosb>=0.0)){ beta=180.0f; gamma=0.0f; }
  else if (zc<-0.99999)               { beta=180.0f; gamma=0.0f; }
  else
  { beta=(float)acos((double)zc)/RADIANF;
    if (fabs((double)xd_)+fabs((double)yd_)<1.0e-5) gamma=0.0f;
    else gamma=(float)atan2((double)yd_,(double)xd_)/RADIANF;
  }
  while (gamma>180.0) gamma-=360.0f;
  while (gamma<-180.0) gamma+=360.0f;

  /* O beta vira INDICE DISCRETO de rotmata (linha 2*beta) -- por isso o
     acos e' em double: 1 ulp aqui e' 1 grau la'. Ver "O risco do beta". */
  beta=(float)floor(((double)beta*10.0+0.5)/10.0);

  float dist=a.devenpar[j*7+3];
  int   kind=a.pairkind[j];
  float ka=a.akin*dist;
  float cosbeta=(float)cos((double)(beta*RADIANF));
  float xdv=d_fvibmsrd(a.thermat,a.thernum,a.therstep,a.mweight,a.tdebye,
    a.tsample,dist,a.aweight[kind]);
  /* mscdrunc.cpp:795 e' akin*akin*xd*(1.0-cosbeta): a associacao a' esquerda
     faz akin*akin*xd ser calculado em FLOAT e so' o ultimo fator promover
     para double. Fazer tudo em double aqui daria outro ultimo bit. */
  float w=a.akin*a.akin*xdv;
  float xb=(float)exp(-0.5*(double)dist*(double)a.xc
    -(double)w*(1.0-(double)cosbeta))/ka;

  float t=ka/RADIANF;
  Gcplx cxa=csmul(xb,d_fexpix(a.cexpix,a.exndata,a.exmdata,t));
  Gcplx cxb_,cxc_,cxd_,cxe_,cxf_,cxg_,cxh_,cxi_;
  cxb_=csmul(xb,d_fexpix(a.cexpix,a.exndata,a.exmdata,t-gamma));
  cxc_=csmul(xb,d_fexpix(a.cexpix,a.exndata,a.exmdata,t+gamma));
  cxd_=csmul(xb,d_fexpix(a.cexpix,a.exndata,a.exmdata,t-gamma-gamma));
  cxe_=csmul(xb,d_fexpix(a.cexpix,a.exndata,a.exmdata,t+gamma+gamma));
  cxf_=csmul(xb,d_fexpix(a.cexpix,a.exndata,a.exmdata,t-gamma-gamma-gamma));
  cxg_=csmul(xb,d_fexpix(a.cexpix,a.exndata,a.exmdata,t+gamma+gamma+gamma));
  cxh_=csmul(xb,d_fexpix(a.cexpix,a.exndata,a.exmdata,
    t-gamma-gamma-gamma-gamma));
  cxi_=csmul(xb,d_fexpix(a.cexpix,a.exndata,a.exmdata,
    t+gamma+gamma+gamma+gamma));

  float vka;
  if (a.raorder<0) vka=0.0f; else vka=1.0f/ka;

  /* Linha de rotmata. V6 mostrou que beta e' inteiro aqui em 30.726.810 de
     30.726.810 chamadas, entao a interpolacao de makerotation e' copia de
     linha e da' para indexar direto.
     UNICA divergencia conhecida com a CPU: em beta=180 a CPU satura i em 359
     e calcula a+1,0*(b-a), que nao e' bit a bit b; aqui le-se a linha 360
     direto. Sao 346 pares em 30,7 milhoes, e o valor certo e' o daqui. */
  int row=(int)((double)fabs((double)beta)*(a.betanum-1.0)/180.0);
  if (row<0) row=0; else if (row>a.betanum-1) row=a.betanum-1;
  const float *rota=a.rotmata+(size_t)row*a.rlnum*a.lamdum;

  /* hanka: vka varia por par, entao a interpolacao de 2 pontos de hankmat
     e' feita inline -- o cache de 100 complexos de fhankelfaca existe na
     CPU so' para amortizar as 15 chamadas, e aqui nao paga o espaco. */
  float ya=(float)(4.0*(a.handata-1.0)*(double)vka);
  int hi=(int)ya;
  if (hi<0) hi=0; else if (hi>a.handata-2) hi=a.handata-2;
  float hf=ya-(float)hi;
  const Gcplx *hrow=a.hankmat_a+(size_t)hi*a.halnum*a.hacmnum;
  const Gcplx *hnxt=hrow+(size_t)a.halnum*a.hacmnum;

  const Gcplx *phc=a.phasec+(size_t)kind*(a.pclnum+1);
  /* rotharma devolve 0 para al>=lnum, entao as parcelas alem das tabelas
     somam zero na CPU. Cortar o laco da' o mesmo resultado e evita ler fora
     de hankmat, que fhankelfaca nao protege. Com lnum=20 nada muda. */
  int alnum=a.alnum[kind];
  if (alnum>a.rlnum) alnum=a.rlnum;
  if (alnum>a.halnum) alnum=a.halnum;

  Gcplx acc[9];
  for (int e=0;e<9;++e) { acc[e].re=0.0f; acc[e].im=0.0f; }

  for (int al=0;al<alnum;++al)
  { Gcplx cpb=phc[al+1];                       /* fsinexpa: akin constante */
    Gcplx hb=a.hankarg_b[al*a.hacmnum];        /* kb==0 em todas as 15 */
    Gcplx ha[3];
    for (int m=0;m<3;++m)
    { Gcplx p=hrow[al*a.hacmnum+m], q=hnxt[al*a.hacmnum+m];
      ha[m].re=p.re+hf*(q.re-p.re);
      ha[m].im=p.im+hf*(q.im-p.im);
    }
    for (int e=0;e<9;++e)
    { int kel=c_kel[c_ma[e]];
      float rot=rota[al*a.lamdum+kel]*a.rotmatc[al*a.lamdum+kel];
      /* mesma ordem da CPU: ((xa*cxb)*cxc)*cxd */
      Gcplx v=cmul(cmul(csmul(rot,cpb),hb),ha[c_na[e]]);
      acc[e]=cadd(acc[e],v);
    }
  }

  if (a.dbg)
  { a.dbg[j*4]=beta; a.dbg[j*4+1]=gamma; a.dbg[j*4+2]=xb;
    a.dbg[j*4+3]=(float)row;
  }

  /* radim encolhe com raorder (mscdruna.cpp:472): 15/10/6/3/1 para 4/3/2/1/0.
     Os guardas sao os mesmos de mscdrunc.cpp:824-895 -- sem eles, qualquer
     raorder<4 escreveria fora de devenelem. */
  Gcplx *o=a.out+(size_t)j*a.radim;
  o[0] =cmul(cxa,acc[0]);
  if (a.raorder>0)
  { o[1] =cmul(cxb_,acc[1]);   o[2] =cmul(cneg(cxc_),acc[1]);
  }
  if (a.raorder>1)
  { o[3] =cmul(cxa,acc[2]);
    o[4] =cmul(cxd_,acc[3]);   o[5] =cmul(cxe_,acc[3]);
  }
  if (a.raorder>2)
  { o[6] =cmul(cxb_,acc[4]);   o[7] =cmul(cneg(cxc_),acc[4]);
    o[8] =cmul(cxf_,acc[5]);   o[9] =cmul(cneg(cxg_),acc[5]);
  }
  if (a.raorder>3)
  { o[10]=cmul(cxa,acc[6]);
    o[11]=cmul(cxd_,acc[7]);   o[12]=cmul(cxe_,acc[7]);
    o[13]=cmul(cxh_,acc[8]);   o[14]=cmul(cxi_,acc[8]);
  }
}

/* ------------------------- lado do host ------------------------- */

static int upf(float **d,const float *h,size_t n)
{ CK(cudaMalloc((void**)d,n*sizeof(float)));
  CK(cudaMemcpy(*d,h,n*sizeof(float),cudaMemcpyHostToDevice)); return 0; }
static int upc(Gcplx **d,const Gcplx *h,size_t n)
{ CK(cudaMalloc((void**)d,n*sizeof(Gcplx)));
  CK(cudaMemcpy(*d,h,n*sizeof(Gcplx),cudaMemcpyHostToDevice)); return 0; }

extern "C" int mscdgpu_setup(const Gconst *k)
{ if (g_ready) return 0;
  K=*k; memset(&D,0,sizeof(D));

  if (upf(&D.patom,k->patom,(size_t)k->natoms*12)) return 1;
  if (upf(&D.devenpar,k->devenpar,(size_t)k->ndbleven*7)) return 1;
  if (upf(&D.rotmata,k->rotmata,
      (size_t)k->betanum*k->rlnum*k->lamdum)) return 1;
  if (upf(&D.rotmatc,k->rotmatc,(size_t)k->rlnum*k->lamdum)) return 1;
  if (upf(&D.thermat,k->thermat,(size_t)k->thernum)) return 1;
  if (upf(&D.aweight,k->aweight,(size_t)k->nkind)) return 1;
  if (upc(&D.hankmat_a,k->hankmat_a,
      (size_t)k->handata*k->halnum*k->hacmnum)) return 1;
  if (upc(&D.hankarg_b,k->hankarg_b,
      (size_t)k->halnum*k->hacmnum)) return 1;
  if (upc(&D.phasec,k->phasec,(size_t)k->nkind*(k->pclnum+1))) return 1;
  if (upc(&D.cexpix,k->cexpix,(size_t)k->exndata)) return 1;

  CK(cudaMalloc((void**)&D.pairgeo,(size_t)k->ndbleven*3*sizeof(float)));
  CK(cudaMalloc((void**)&D.pairkind,(size_t)k->ndbleven*sizeof(int)));
  CK(cudaMalloc((void**)&D.devenelem,
    (size_t)k->ndbleven*k->radim*sizeof(Gcplx)));
  CK(cudaMalloc((void**)&D.alnum,(size_t)k->nkind*sizeof(int)));
  CK(cudaMalloc((void**)&D.dbg,(size_t)k->ndbleven*4*sizeof(float)));

  int nb=(k->ndbleven+255)/256;
  k_pairgeo<<<nb,256>>>(D.patom,D.devenpar,k->ndbleven,D.pairgeo,D.pairkind);
  CK(cudaGetLastError());
  CK(cudaDeviceSynchronize());
  g_ready=1; return 0;
}

extern "C" int mscdgpu_set_alnum(const int *alnum,int nkind)
{ CK(cudaMemcpy(D.alnum,alnum,(size_t)nkind*sizeof(int),
    cudaMemcpyHostToDevice)); return 0; }

/* hankb e' refotografado a cada ponto porque o cache de fhankelfaca e'
   keyed no argumento e outras funcoes do laco mexem nele. 100 complexos. */
extern "C" int mscdgpu_set_hankb(const Gcplx *h)
{ CK(cudaMemcpy(D.hankarg_b,h,(size_t)K.halnum*K.hacmnum*sizeof(Gcplx),
    cudaMemcpyHostToDevice)); return 0; }

extern "C" int mscdgpu_alldblevent(float akin,const float *xdetec,float xc,
  Gcplx *out)
{ Kargs a;
  if (!g_ready) { snprintf(g_err,sizeof(g_err),"setup nao chamado"); return 1; }

  /* (patomc-patomb) == xdetec para todos os pares: uma vez por ponto. */
  float xb=xdetec[0],yb=xdetec[1],zb=xdetec[2];
  double lengb=sqrt((double)(xb*xb+yb*yb+zb*zb));
  if (lengb<1.0e-5) { snprintf(g_err,sizeof(g_err),"lengb=0"); return 1; }
  float cosb=(float)(zb/lengb),sinb,phib;
  if (fabs((double)xb)+fabs((double)yb)<1.0e-5) phib=0.0f;
  else phib=(float)atan2((double)yb,(double)xb);
  if (fabs((double)cosb)>=1.0) sinb=0.0f;
  else sinb=(float)sqrt(1.0-(double)(cosb*cosb));   /* float, ver k_pairgeo */

  a.devenpar=D.devenpar; a.pairgeo=D.pairgeo; a.rotmata=D.rotmata;
  a.rotmatc=D.rotmatc;   a.thermat=D.thermat; a.aweight=D.aweight;
  a.hankmat_a=D.hankmat_a; a.hankarg_b=D.hankarg_b; a.phasec=D.phasec;
  a.cexpix=D.cexpix; a.pairkind=D.pairkind; a.alnum=D.alnum;
  a.out=D.devenelem;
  a.ndbleven=K.ndbleven; a.rlnum=K.rlnum; a.lamdum=K.lamdum;
  a.betanum=K.betanum;   a.handata=K.handata; a.halnum=K.halnum;
  a.hacmnum=K.hacmnum;   a.pclnum=K.pclnum;
  a.exndata=K.exndata;   a.exmdata=K.exmdata; a.thernum=K.thernum;
  a.radim=K.radim;       a.raorder=K.raorder;
  a.therstep=K.therstep; a.mweight=K.mweight; a.tdebye=K.tdebye;
  a.tsample=K.tsample;
  a.akin=akin; a.xc=xc; a.cosb=cosb; a.sinb=sinb; a.phib=phib;
  a.dbg=g_dbg?D.dbg:NULL;

  int nb=(K.ndbleven+127)/128;
  k_alldblevent<<<nb,128>>>(a);
  CK(cudaGetLastError());
  CK(cudaMemcpy(out,D.devenelem,
    (size_t)K.ndbleven*K.radim*sizeof(Gcplx),cudaMemcpyDeviceToHost));
  return 0;
}

extern "C" void mscdgpu_set_dbg(int on) { g_dbg=on; }

extern "C" int mscdgpu_get_dbg(float *out)
{ CK(cudaMemcpy(out,D.dbg,(size_t)K.ndbleven*4*sizeof(float),
    cudaMemcpyDeviceToHost)); return 0; }

extern "C" void mscdgpu_teardown(void)
{ if (!g_ready) return;
  cudaFree(D.patom); cudaFree(D.devenpar); cudaFree(D.rotmata);
  cudaFree(D.rotmatc); cudaFree(D.thermat); cudaFree(D.aweight);
  cudaFree(D.hankmat_a); cudaFree(D.hankarg_b); cudaFree(D.phasec);
  cudaFree(D.cexpix); cudaFree(D.pairgeo); cudaFree(D.pairkind);
  cudaFree(D.devenelem); cudaFree(D.alnum);
  memset(&D,0,sizeof(D)); g_ready=0;
}

extern "C" const char *mscdgpu_lasterror(void) { return g_err; }
