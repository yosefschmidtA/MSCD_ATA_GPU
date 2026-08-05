#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "cartesia.h"
#include "phase.h"
#include "radmat.h"
#include "meanpath.h"
#include "pdinten.h"
#include "fcomplex.h"
#include "vibrate.h"
#include "msfuncs.h"
#include "rotamat.h"
#include "pdchifit.h"
#include "jobtime.h"
#include "userutil.h"
#include "mscdrun.h"
#include "mscdtimer.h"

//beta unit: degree
Fcomplex Mscdrun::evenelem(int akind,int ma,int na,int mb,int nb,
  float akin,float vka,float vkb,float beta)
{ int al,ka,kb,kelem,kharm,alnum;
  float xa;
  Fcomplex cxa,cxb,cxc,cxd;

  if (error==0)
  { ka=iabs(na); kb=iabs(mb)+iabs(nb);
    kelem=evenmat->getkelem(ma,mb);
    kharm=evenmat->getkharm(ma,mb);
    alnum=phaseshift[akind-1].getalnum(akin);
    cxa=0.0f;
    for (al=0;al<alnum;++al)
    { xa=evenmat->rotharma(al,kelem,kharm,beta);
      cxb=phaseshift[akind-1].fsinexpa(akin,al);
      cxc=hankb->fhankelfaca(al,kb,vkb);
      cxd=hanka->fhankelfaca(al,ka,vka);
      cxa+=xa*cxb*cxc*cxd;
    }
  }
  return(cxa);
} //end of Mscdrun::evenelem

//beta unit: degree
Fcomplex Mscdrun::ATAevenelem(int akind,int ma,int na,int mb,int nb,
  float akin,float vka,float vkb,float beta, float ATAweight1,
  float ATAweight2)
{ int al,ka,kb,kelem,kharm,alnum;
  float xa;
  Fcomplex cxa,cxb,cxb1,cxb2,cxc,cxd;

  if (error==0)
  { if (akind==3)
    {  ka=iabs(na); kb=iabs(mb)+iabs(nb);
       kelem=evenmat->getkelem(ma,mb);
       kharm=evenmat->getkharm(ma,mb);
       alnum=phaseshift[akind-1].getalnum(akin);
       cxa=0.0f;
       for (al=0;al<alnum;++al)
       {  xa=evenmat->rotharma(al,kelem,kharm,beta);
          cxb1=phaseshift[akind-3].fsinexpa(akin,al);
          cxb2=phaseshift[akind-2].fsinexpa(akin,al);
          cxb=(1.0 - ATAweight1)*cxb1 + ATAweight1*(cxb2);
          cxc=hankb->fhankelfaca(al,kb,vkb);
          cxd=hanka->fhankelfaca(al,ka,vka);
          cxa+=xa*cxb*cxc*cxd;
       }
    }
    else
    {  if (akind==4)
       {  ka=iabs(na); kb=iabs(mb)+iabs(nb);
          kelem=evenmat->getkelem(ma,mb);
          kharm=evenmat->getkharm(ma,mb);
          alnum=phaseshift[akind-1].getalnum(akin);
          cxa=0.0f;
          for (al=0;al<alnum;++al)
          {  xa=evenmat->rotharma(al,kelem,kharm,beta);
             cxb1=phaseshift[akind-4].fsinexpa(akin,al);
             cxb2=phaseshift[akind-3].fsinexpa(akin,al);
             cxb=(1.0 - ATAweight2)*cxb1 + ATAweight2*(cxb2);
             cxc=hankb->fhankelfaca(al,kb,vkb);
             cxd=hanka->fhankelfaca(al,ka,vka);
             cxa+=xa*cxb*cxc*cxd;
          }
       }
       else
       {  ka=iabs(na); kb=iabs(mb)+iabs(nb);
          kelem=evenmat->getkelem(ma,mb);
          kharm=evenmat->getkharm(ma,mb);
          alnum=phaseshift[akind-1].getalnum(akin);
          cxa=0.0f;
          for (al=0;al<alnum;++al)
          { xa=evenmat->rotharma(al,kelem,kharm,beta);
            cxb=phaseshift[akind-1].fsinexpa(akin,al);
            cxc=hankb->fhankelfaca(al,kb,vkb);
            cxd=hanka->fhankelfaca(al,ka,vka);
            cxa+=xa*cxb*cxc*cxd;
          }
       }
    }
  }
  return(cxa);
} //end of Mscdrun::ATAevenelem

//beta unit: degree
Fcomplex Mscdrun::evenbelem(int alf,int ma,int na,int mb,int nb,
  float vkb,float beta)
{ float xb;
  Fcomplex cxa,cxb,cxc,cxd;

  xb=(float)na;
  if (error==0)
  { xb=termmat->termination(alf,ma,mb,beta);
    cxc=hankb->fhankelfac(alf,iabs(mb)+iabs(nb),vkb);
    cxa=xb*cxc;
  }
  return(cxa);
} //end of Mscdrun::evenbelem

/*
----------------------------------------------------------------------
 calculate all the multiple scattering event

p   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14
ma  0   1  -1   0   2  -2   1  -1   3  -3   0   2  -2   4  -4
na  0   0   0   1   0   0   1   1   0   0   2   1   1   0   0

q   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14
mb  0   1  -1   0   2  -2   1  -1   3  -3   0   2  -2   4  -4
nb  0   0   0   1   0   0   1   1   0   0   2   1   1   0   0

----------------------------------------------------------------------
*/
int Mscdrun::alltrievent(int forcut,float akin)
{ int j,k,p,q,t,pt,eledim,peledim,akind,memadd,pmemadd,ma,na,mb,nb,
    smear;
  float xa,xb,xc,xd,xe,xf,cosbeta,ka,vka,vkb,beta;
  Fcomplex cxa,cxb,cxc,cvalue;
  int lamda[64];
  Fcomplex *tempeven;
  const float radian=(float)(3.14159265/180.0);

  tempeven=NULL;
  if ((error==0)&&(msorder>1)&&(!tevenelem)) error=901;
  if (error==0)
  { tempeven=new Fcomplex [radim*radim*2];
    if (!tempeven) error=102;
  }
  for (j=0;(error==0)&&(j<32);++j)
  { if ((j==0)||(j==3)||(j==10)) ma=0;
    else if (j<3) ma=3-j*2;
    else if (j<6) ma=18-j*4;
    else if (j<8) ma=13-j*2;
    else if (j<10) ma=51-j*6;
    else if (j<13) ma=46-j*4;
    else if (j<15) ma=108-j*8;
    else ma=0;

    if (j==10) na=2;
    else if ((j==3)||(j==6)||(j==7)||(j==11)||(j==12)) na=1;
    else na=0;

    lamda[j]=ma; lamda[32+j]=na;
  }

  for (j=0;(error==0)&&(j<ntrieven);++j)
  { if ((j>0)&&(tevenpar[j*10+1]==tevenpar[(j-1)*10+1])&&
      (tevenpar[j*10+2]==tevenpar[(j-1)*10+2])&&
      (tevenpar[j*10+3]==tevenpar[(j-1)*10+3])&&
      (tevenpar[j*10+4]==tevenpar[(j-1)*10+4])&&
      (tevenpar[j*10]>=tevenpar[(j-1)*10])&&
      (tevenpar[j*10+5]<=tevenpar[(j-1)*10+5]))
      smear=1;
    else smear=0;
    xa=tevenpar[j*10]; ka=akin*xa; cosbeta=tevenpar[j*10+3];
    akind=(int)tevenpar[j*10+4]; eledim=(int)tevenpar[j*10+5];
    memadd=(int)tevenpar[j*10+6];
    if (j>0)
    { peledim=(int)tevenpar[(j-1)*10+5];
      pmemadd=(int)tevenpar[(j-1)*10+6];
    }
    else peledim=pmemadd=0;
    xc=meanpath->finvpath(akin);
    xd=vibrate->fvibmsrd(xa,aweight[akind-1]);
    xb=(float)exp(-0.5*xa*xc-akin*akin*xd*(1.0-cosbeta))/ka;
    cxa=xb*expix->fexpix(ka/radian);
    if (smear==0)
    { cxb=0.0f; cxc=cxa; beta=(float)acos(cosbeta)/radian;
      beta=(float)floor(beta*10+0.5)*0.1f;
    }
    else
    { cxb=cxa/cxc; cxc=cxa; beta=0.0f;
    }
    if (raorder<0) vka=vkb=0.0f;
    else
    { vka=tevenpar[j*10+1]/akin; vkb=tevenpar[j*10+2]/akin;
    }

    for (p=0;p<eledim;++p)
    { for (q=0;q<eledim;++q)
      { t=q+p*eledim; pt=pmemadd+q+p*peledim;
        ma=lamda[p]; na=lamda[32+p];
        mb=lamda[q]; nb=lamda[32+q];
        k=(ma+mb)&1;
        if ((ma==0)&&(mb<0)&&(k==0))
          tempeven[t]=tempeven[t-1];
        else if ((ma==0)&&(mb<0))
          tempeven[t]=-tempeven[t-1];
        else if ((ma<0)&&(mb==0)&&(k==0))
          tempeven[t]=tempeven[t-eledim];
        else if ((ma<0)&&(mb==0))
          tempeven[t]=-tempeven[t-eledim];
        else if ((ma<0)&&(mb>0)&&(k==0))
          tempeven[t]=tempeven[t-eledim+1];
        else if ((ma<0)&&(mb>0))
          tempeven[t]=-tempeven[t-eledim+1];
        else if ((ma<0)&&(mb<0)&&(k==0))
          tempeven[t]=tempeven[t-eledim-1];
        else if ((ma<0)&&(mb<0))
          tempeven[t]=-tempeven[t-eledim-1];
        else if (smear==0)
        { if (ATA==1)
          {  cvalue=ATAevenelem(akind,ma,na,mb,nb,akin,vka,vkb,beta,ATAweight1,                                 ATAweight2);
             tempeven[t]=cxa*cvalue;
          }
          else
          {
             cvalue=evenelem(akind,ma,na,mb,nb,akin,vka,vkb,beta);
             tempeven[t]=cxa*cvalue;
          }
        }
        else if (forcut==0) tempeven[t]=cxb*tevenelem[pt];
        else tempeven[t]=cxb*tempeven[t+radim*radim];
      }
    }

    if (forcut==0)
    { for (p=0;p<eledim;++p)
      { for (q=0;q<eledim;++q)
        { t=q+p*eledim;
          tevenelem[memadd+t]=tempeven[t];
        }
      }
    }
    else
    { xa=cabs(tempeven[0]); xb=xc=xd=xe=0.0f;
      for (p=0;p<eledim;++p)
      { for (q=0;q<eledim;++q)
        { t=q+p*eledim;
          ma=lamda[p]; na=lamda[32+p];
          mb=lamda[q]; nb=lamda[32+q];

          if (ma<0) continue;
          else xf=cabs(tempeven[t]);

          for (k=0;k<ma+na;++k) xf*=vka;
          for (k=0;k<nb;++k) xf*=vkb;
          if ((p<1)&&(q<1)) continue;
          else if ((p<3)&&(q<3))
          { if (xb<xf) xb=xf;
          }
          else if ((p<6)&&(q<6))
          { if (xc<xf) xc=xf;
          }
          else if ((p<10)&&(q<10))
          { if (xd<xf) xd=xf;
          }
          else if (xe<xf) xe=xf;

          tempeven[t+radim*radim]=tempeven[t];
        }
      }
      if (xd<xe) xd=xe; if (xc<xd) xc=xd;
      if (xb<xc) xb=xc; if (xa<xb) xa=xb;
      tevenelem[memadd]=tempeven[0]; tevenelem[memadd+1]=xa;
      tevenelem[memadd+2]=xb; tevenelem[memadd+3]=xc;
      tevenelem[memadd+4]=xd; tevenelem[memadd+5]=xe;
    }
  }
  if (tempeven) delete [] tempeven;

  return(error);
} //end of Mscdrun::alltrievent

int Mscdrun::precutable()
{ int ia,ib,ic,id,j,k,m,memadd,eledim;
  float xa,xb,xc,xd,xe,xf,akin,pemeven,centa,centb,centc,centd;
  float *stat;
  Fcomplex *asum,*bsum;
  Fcomplex cxa;

  MSCDT_DECL;
  Textout conout;
  if (error==0)
  { if (dispmode>4)
      conout.string("calculating vibrational coefficients ...",0,16+1);
    if ((displog>0)&&(flogout))
      flogout->string("calculating vibrational coefficients ...",
        0,16+1);
    if (vibrate) delete vibrate; vibrate=new Vibration;
    if (error==0) vibrate->init();
    if (error==0)
      error=vibrate->loadparameter(density,mweight,tdebye,tsample);
  }
  MSCDT("  precut vibracionais");
  if (error==0)
  { if (dispmode>4)
      conout.string("calculating rotational matrices ...",0,1);
    if ((displog>0)&&(flogout))
      flogout->string("calculating rotational matrices ...",0,1);
    if (evenmat) delete evenmat; evenmat=new Rotamat;
    if (termmat) delete termmat;
    if (raorder>linitial) k=raorder; else k=linitial;
    termmat=new Rotamat;
    if (error==0) evenmat->init(lnum,raorder);
    if (error==0) termmat->init(linitial+2,k);
    if (error==0) error=evenmat->makecurve();
    if (error==0) error=termmat->makecurve();
  }
  MSCDT("  precut rotacionais");
  if (error==0)
  { if (dispmode>4)
      conout.string("calculating spherical expension coefficients ...",
        0,2);
    if ((displog>0)&&(flogout)) flogout->string(
      "calculating spherical expension coefficients ...",0,2);
    if (expix) delete expix;
    if (hanka) delete hanka; if (hankb) delete hankb;
    expix=new Expix;
    if (error==0) expix->init(3601);
    if (error==0) error=expix->makecurve();
    k=raorder+1; if (k<1) k=1;
    if (lnum>linitial+2) j=lnum; else j=linitial+2;
    hanka=new Hankel; hankb=new Hankel;
    if (error==0) hanka->init(101,j,k);
    if (error==0) hankb->init(101,j,k);
    if (error==0) error=hanka->makecurve();
    if (error==0) error=hankb->makecurve();
  }

  asum=bsum=0; stat=NULL;
  ntrielem=ntrieven*6;
  if (error==0) stat=new float [21*21];
  if ((error==0)&&(!stat)) error=102;
  if ((error==0)&&(msorder>1))
  { if (tevencut) delete [] tevencut; if (tevendim) delete [] tevendim;
    if (tevenelem) delete [] tevenelem;
    tevencut=new int [msorder*natoms*natoms];
    tevendim=new int [natoms*natoms*natoms];
    tevenelem=new Fcomplex [ntrielem];
    asum=new Fcomplex [radim*natoms*natoms];
    bsum=new Fcomplex [radim*natoms*natoms];
    if ((!tevencut)||(!tevendim)||(!tevenelem)||(!asum)||(!bsum))
      error=102;
  }
  if ((error==0)&&(msorder>1))
  { for (j=0;j<msorder*natoms*natoms;++j) tevencut[j]=0;
    for (j=0;j<natoms*natoms*natoms;++j) tevendim[j]=0;
  }

  for (m=0;(error==0)&&(m<ntrieven);++m)
  { k=(int)tevenpar[m*10+4];
    if ((k<1)||(k>katoms)) error=621;
    else
    { tevenpar[m*10+5]=(float)radim; tevenpar[m*10+6]=(float)m*6;
    }
  }
  if ((msorder<2)||(pathcut<1.0e-10)) pathcut=0.0f;
  else if (error==0)
  { if (dispmode>4)
      conout.string("setting up pre-cut table ...",0,2);
    if ((displog>0)&&(flogout))
      flogout->string("setting up pre-cut table ...",0,2);
    akin=kmin;
    MSCDT("  precut esfericas+cortes");
    error=alltrievent(1,akin);
    MSCDT("  precut alltrievent");
  }
  pemeven=0.0f;
  if ((error==0)&&(pathcut>0.0))
  {
#pragma omp parallel for collapse(2) schedule(static) \
    private(ia,id,k,eledim,memadd,xb) reduction(max:pemeven)
    for (ib=0;ib<natoms;++ib)
    { for (ic=0;ic<natoms;++ic)
      { if (ic==ib) continue;
        for (ia=0;ia<natoms;++ia)
        { if ((ia==ib)||(patom[ia*12+7]==0.0)) continue;
          id=ia*natoms*natoms+ib*natoms+ic;
          k=tevenadd[id];
          eledim=(int)tevenpar[k*10+5]; memadd=(int)tevenpar[k*10+6];
          if (eledim<1) continue;
          xb=real(tevenelem[memadd+1]);
          if (pemeven<xb) pemeven=xb;
        }
      }
    }
  }

  if ((pathcut>0.0)&&(pemeven<1.0e-10)) pemeven=1.0f;
  if ((error==0)&&(pathcut>0.0))
  { for (ib=0;ib<natoms;++ib)
    { for (ic=0;ic<natoms;++ic) bsum[ib*natoms+ic]=1.0f/pemeven;
    }
  }
  MSCDT("  precut maximo do pemeven");
  /* A ordem (ib,ic,ia) parece ruim -- o indice e' id=ia*natoms^2+..., entao
     cada passo de ia salta natoms^2 = 244 KB dentro de tevenadd/tevendim.
     MEDIDO: inverter para (ia,ib,ic) deixa 27% MAIS LENTO (8,6 s -> 10,9 s).
     O motivo e' que xa e cxa vivem em registrador durante todo o laco de ia;
     com ia por fora eles viram acumuladores em memoria, e as 105 milhoes de
     iteracoes passam a fazer load+store. Nao inverta sem medir de novo.
     Detalhes em OTIMIZACAO.md, secao "V3". */
  for (m=2;(error==0)&&(m<=msorder);++m)
  { for (ib=0;ib<natoms;++ib)
    { for (ic=0;ic<natoms;++ic) asum[ib*natoms+ic]=bsum[ib*natoms+ic];
    }
    /* O laco de ia FICA sequencial dentro da thread de proposito: cxa e' uma
       soma em float, e so' assim a ordem das parcelas continua a mesma da
       versao serial -- que e' a condicao para a curva sair identica bit a bit.
       Reduzir cxa entre threads mudaria o ultimo bit, e a saida deste laco e'
       uma MASCARA discreta (xb>pathcut): um bit vira um trio que entra ou sai
       da conta, nao um erro pequeno.
       A corrida em tevencut e' real e benigna: todos os ic de um mesmo (ia,ib)
       escrevem o mesmo slot, mas so' 0->1, e o valor nao decide nada aqui. */
#pragma omp parallel for collapse(2) schedule(static) \
    private(ia,id,k,memadd,xa,xb,xc,xd,xe,xf,cxa)
    for (ib=0;ib<natoms;++ib)
    { for (ic=0;ic<natoms;++ic)
      { if (ic==ib) continue;
        xa=0.0f; cxa=0.0f;
        for (ia=0;ia<natoms;++ia)
        { if ((ia==ib)||((m==2)&&(patom[ia*12+7]==0.0))) continue;
          id=ia*natoms*natoms+ib*natoms+ic;
          if (pathcut>1.0e-10)
          { k=tevenadd[id];
            if ((tevenpar[k*10+1]>tevenpar[k*10+2])&&(m>2))
              k=tevenadd[ic*natoms*natoms+ib*natoms+ia];
            memadd=(int)tevenpar[k*10+6];
            xb=cabs(asum[ia*natoms+ib])*real(tevenelem[memadd+1]);
            xc=cabs(asum[ia*natoms+ib])*real(tevenelem[memadd+2]);
            xd=cabs(asum[ia*natoms+ib])*real(tevenelem[memadd+3]);
            xe=cabs(asum[ia*natoms+ib])*real(tevenelem[memadd+4]);
            xf=cabs(asum[ia*natoms+ib])*real(tevenelem[memadd+5]);
            if ((m==2)&&(xa<xb)) xa=xb;
            else if (m>2) cxa+=asum[ib*natoms+ic]*tevenelem[memadd];
          }
          else
          { xb=pathcut+1.0f;
            if (raorder>0) xc=pathcut+1.0f; else xc=0.0f;
            if (raorder>1) xd=pathcut+1.0f; else xd=0.0f;
            if (raorder>2) xe=pathcut+1.0f; else xe=0.0f;
            if (raorder>3) xf=pathcut+1.0f; else xf=0.0f;
          }
          if ((raorder>3)&&(xf>pathcut))
          { if (tevencut[(m-1)*natoms*natoms+ia*natoms+ib]==0)
              tevencut[(m-1)*natoms*natoms+ia*natoms+ib]=1;
            if (m<=8) tevendim[id]+=(15<<((m-2)*4));
          }
          else if ((raorder>2)&&(xe>pathcut))
          { if (tevencut[(m-1)*natoms*natoms+ia*natoms+ib]==0)
              tevencut[(m-1)*natoms*natoms+ia*natoms+ib]=1;
            if (m<=8) tevendim[id]+=(10<<((m-2)*4));
          }
          else if ((raorder>1)&&(xd>pathcut))
          { if (tevencut[(m-1)*natoms*natoms+ia*natoms+ib]==0)
              tevencut[(m-1)*natoms*natoms+ia*natoms+ib]=1;
            if (m<=8) tevendim[id]+=(6<<((m-2)*4));
          }
          else if ((raorder>0)&&(xc>pathcut))
          { if (tevencut[(m-1)*natoms*natoms+ia*natoms+ib]==0)
              tevencut[(m-1)*natoms*natoms+ia*natoms+ib]=1;
            if (m<=8) tevendim[id]+=(3<<((m-2)*4));

          }
          else if (xb>pathcut)
          { if (tevencut[(m-1)*natoms*natoms+ia*natoms+ib]==0)
              tevencut[(m-1)*natoms*natoms+ia*natoms+ib]=1;
            if (m<=8) tevendim[id]+=(1<<((m-2)*4));
          }
        }
        if (m==2) bsum[ib*natoms+ic]=(float)pow((double)xa,0.75);
        else bsum[ib*natoms+ic]=cxa;
      }
    }
  }

  MSCDT("  precut pathcut msorder x natoms^3");
  for (j=0;(error==0)&&(j<ntrieven);++j) tevenpar[j*10+5]=0.0f;
  for (j=0;(error==0)&&(j<21*21);++j) stat[j]=0.0f;
  for (m=2;(error==0)&&(m<=msorder);++m)
  { for (ib=0;ib<natoms;++ib)
    { for (ic=0;ic<natoms;++ic)
      { if (ic==ib) continue;
        for (ia=0;ia<natoms;++ia)
        { if ((ia==ib)||((m==2)&&(patom[ia*12+7]==0.0))) continue;
          id=ia*natoms*natoms+ib*natoms+ic;
          k=tevenadd[id]; eledim=tevendim[id];
          if ((sizeint<4)&&(m>5)) eledim>>=12;
          else if (m>8) eledim>>24;
          else eledim>>=(m-2)*4;
          eledim&=15;
          if (tevencut[(m-1)*natoms*natoms+ia*natoms+ib]!=0)
          { if (tevenpar[k*10+5]<eledim)
              tevenpar[k*10+5]=(float)eledim;
            stat[(m-1)*21+eledim]+=1.0f;
          }
        }
      }
    }
  }
  ntrielem=0;
  for (j=0;(error==0)&&(j<ntrieven);++j)
  { eledim=(int)tevenpar[j*10+5]; tevenpar[j*10+6]=(float)ntrielem;
    ntrielem+=eledim*eledim; stat[msorder*21+eledim]+=1.0f;
  }

  for (m=1;(error==0)&&(m<=msorder);++m)
  { for (j=0;j<6;++j)
    { if ((m==1)&&(raorder>=0)&&(j==raorder+1))
        stat[(m-1)*21+j]=100.0f;
      else if ((m==1)&&(raorder<0)&&(j==1)) stat[(m-1)*21+j]=100.0f;
      else if (m==1) stat[(m-1)*21+j]=0.0f;
      else if ((m==2)&&(j<2)) stat[(m-1)*21+j]*=
        (float)(100.0/eatoms/(natoms-1.0)/(natoms-1.0));
      else if ((m==2)&&(j==2)) stat[(m-1)*21+j]=(float)
        (stat[(m-1)*21+3]*100.0/eatoms/(natoms-1.0)/(natoms-1.0));
      else if ((m==2)&&(j==3)) stat[(m-1)*21+j]=(float)
        (stat[(m-1)*21+6]*100.0/eatoms/(natoms-1.0)/(natoms-1.0));
      else if ((m==2)&&(j==4)) stat[(m-1)*21+j]=(float)
        (stat[(m-1)*21+10]*100.0/eatoms/(natoms-1.0)/(natoms-1.0));
      else if ((m==2)&&(j==5)) stat[(m-1)*21+j]=(float)
        (stat[(m-1)*21+15]*100.0/eatoms/(natoms-1.0)/(natoms-1.0));
      else if (j<2) stat[(m-1)*21+j]*=
        (float)(100.0/natoms/(natoms-1.0)/(natoms-1.0));
      else if (j==2) stat[(m-1)*21+j]=(float)
        (stat[(m-1)*21+3]*100.0/natoms/(natoms-1.0)/(natoms-1.0));
      else if (j==3) stat[(m-1)*21+j]=(float)
        (stat[(m-1)*21+6]*100.0/natoms/(natoms-1.0)/(natoms-1.0));
      else if (j==4) stat[(m-1)*21+j]=(float)
        (stat[(m-1)*21+10]*100.0/natoms/(natoms-1.0)/(natoms-1.0));
      else if (j==5) stat[(m-1)*21+j]=(float)
        (stat[(m-1)*21+15]*100.0/natoms/(natoms-1.0)/(natoms-1.0));
    }
  }
  for (m=1;(error==0)&&(m<=msorder);++m)
  { xa=0.0f;
    for (j=1;j<6;++j) xa+=stat[(m-1)*21+j];
    if (xa>100.0) stat[(m-1)*21]=0.0f;
    else stat[(m-1)*21]=100.0f-xa;
  }

  for (j=0;(error==0)&&(j<6);++j)
  { if ((ntrieven>0)&&(j<2))
      stat[msorder*21+j]*=(float)(100.0/ntrieven);
    else if ((ntrieven>0)&&(j==2))
      stat[msorder*21+j]=(float)(stat[msorder*21+3]*100.0/ntrieven);
    else if ((ntrieven>0)&&(j==3))
      stat[msorder*21+j]=(float)(stat[msorder*21+6]*100.0/ntrieven);
    else if ((ntrieven>0)&&(j==4))
      stat[msorder*21+j]=(float)(stat[msorder*21+10]*100.0/ntrieven);
    else if ((ntrieven>0)&&(j==5))
      stat[msorder*21+j]=(float)(stat[msorder*21+15]*100.0/ntrieven);
  }
  if ((error==0)&&(msorder>2)&&(natoms>1))
    centa=(float)(ntrieven*100.0/natoms/(natoms-1.0)/(natoms-1.0));
  else if ((error==0)&&(msorder==2)&&(natoms>1)&&(eatoms>0))
    centa=(float)(ntrieven*100.0/eatoms/(natoms-1.0)/(natoms-1.0));
  else centa=0.0f;
  if ((error==0)&&(msorder>1)&&(natoms>1))
    centb=(float)(ndbleven*100.0/natoms/(natoms-1.0));
  else if ((error==0)&&(natoms>1)&&(eatoms>0))
    centb=(float)(ndbleven*100.0/eatoms/(natoms-1.0));
  else centb=0.0f;
  if ((error==0)&&(ntrieven>0)&&(radim>0))
  { centc=(float)(nscorse*100.0/ntrieven);
    centd=(float)(ntrielem*100.0/ntrieven/radim/radim);
  }
  else centc=centd=0.0f;

  if ((error==0)&&(dispmode>2))
  { xa=getmemory();
    if (xa>0.0) xa*=1.0e-6f;
    if ((xa>0.0)&&(sizeint<4)) xa+=0.2f;
    else if (xa>0.0) xa+=0.8f;
    if ((numpe>1)&&(xa>0.0)) xa*=(float)(numpe*2.0);
    if (xa>0.0)
    { conout.string("This job allocated ",0,16);
      conout.floating(xa,6.1f,256);
      conout.string(" megabytes memory in ");
      conout.integer(numpe,4);
      if (numpe<2) conout.string(" processor",0,2);
      else conout.string(" processors",0,2);
    }
  }
  if ((error==0)&&(dispmode>4))
  { conout.floating(centa,8.2f,256); conout.floating(centb,8.2f,256);
    conout.floating(centc,8.2f,256); conout.floating(centd,8.2f,256);
    conout.string("   trieven dbleven corse trielem",0,1);
  }
  if ((error==0)&&(dispmode>4)&&(msorder>1))
  { for (m=2;m<=msorder;++m)
    { if (m>8) continue;
      for (j=0;j<6;++j) conout.floating(stat[(m-1)*21+j],7.1f,256);
      conout.string("   event percent (msorder=");
      conout.integer(m); conout.string(")",0,1);
    }
  }
  if ((error==0)&&(dispmode>4))
  { for (j=0;j<6;++j) conout.floating(stat[msorder*21+j],7.1f,256);
    conout.string("   unique event percent",0,1);
    waitenter();
  }

  if ((error==0)&&(displog>0)&&(flogout))
  { xa=getmemory();
    if (xa>0.0) xa*=1.0e-6f;
    if ((xa>0.0)&&(sizeint<4)) xa+=0.2f;
    else if (xa>0.0) xa+=0.8f;
    if (xa>0.0)
    { flogout->string("This job allocated ",0,16);
      flogout->floating(xa,6.1f,256);
      flogout->string(" megabytes memory in ");
      flogout->integer(numpe,4);
      if (numpe<2) flogout->string(" processor",0,2);
      else flogout->string(" processors",0,2);
    }
  }
  if ((error==0)&&(displog>0)&&(flogout))
  { for (m=2;m<=msorder;++m)
    { if (m>8) continue;
      for (j=0;j<6;++j) flogout->floating(stat[(m-1)*21+j],7.1f,256);
      flogout->string("   event percent (msorder=");
      flogout->integer(m); flogout->string(")",0,1);
    }
    for (j=0;j<6;++j) flogout->floating(stat[msorder*21+j],7.1f,256);
    flogout->string("   unique event percent",0,1);
  }
  if ((error==0)&&(displog>0)&&(flogout)) error=flogout->geterror();

  if (stat) delete [] stat;
  if (asum) delete [] asum; if (bsum) delete [] bsum;
  if (pdnum) delete [] pdnum;
  if (tevenelem) delete [] tevenelem;
  if (devenelem) delete [] devenelem;
  if (devendetec) delete [] devendetec;
  if (talpha) delete [] talpha; if (tgamma) delete [] tgamma;
  if (fithist) delete [] fithist;
  if (error==0)
  { pdnum=new int [npoint];
    tevenelem=new Fcomplex [ntrielem];
    devenelem=new Fcomplex [ndbleven*radim];
    devendetec=new Fcomplex [natoms*natoms*radim];
    if ((!pdnum)||(!tevenelem)||(!devenelem)||(!devendetec)) error=102;
  }
  if ((error==0)&&(msorder>1)&&(raorder>0))
  { talpha=new float [natoms*natoms*natoms];
    tgamma=new float [natoms*natoms*natoms];
    if ((!talpha)||(!tgamma)) error=102;
  }
  if ((error==0)&&(nfit>0))
  { fithist=new float [(trymax+50)*(nfit+10)];
    if (!fithist) error=102;
  }

  MSCDT("  precut estatisticas+alocacoes");
  return(error);
} //end of Mscdrun::precutable

// alpha, beta and gamma use unit of degree
int Mscdrun::onerotation(float *patoma,float *patomb,float *patomc,
  float *alpha,float *beta,float *gamma)
{ float xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,lenga,lengb,sina,cosa,
    phia,sinb,cosb,phib,phiab,sinphi,cosphi,aalpha,abeta,agamma;
  const float radian=(float)(3.14159265/180.0);

  if (error==0)
  { xa=patomb[0]-patoma[0]; ya=patomb[1]-patoma[1];
    za=patomb[2]-patoma[2];
    xb=patomc[0]-patomb[0]; yb=patomc[1]-patomb[1];
    zb=patomc[2]-patomb[2];
    lenga=(float)sqrt(xa*xa+ya*ya+za*za);
    lengb=(float)sqrt(xb*xb+yb*yb+zb*zb);
    if ((lenga<1.0e-5)||(lengb<1.0e-5)) error=605;
  }
  else xa=ya=za=lenga=xb=yb=zb=lengb=0.0f;

  if (error==0)
  { cosa=za/lenga;
    if (fabs(xa)+fabs(ya)<1.0e-5) phia=0.0f;
    else phia=(float)atan2(ya,xa);
    if (fabs(cosa)>=1.0) sina=0.0f;
    else sina=(float)sqrt(1.0-cosa*cosa);

    cosb=zb/lengb;
    if (fabs(xb)+fabs(yb)<1.0e-5) phib=0.0f;
    else phib=(float)atan2(yb,xb);
    if (fabs(cosb)>=1.0) sinb=0.0f;
    else sinb=(float)sqrt(1.0-cosb*cosb);

    phiab=phia-phib;
    cosphi=(float)cos(phiab); sinphi=(float)sin(phiab);

    xc=sinb*cosa-cosb*sina*cosphi; yc=-sina*sinphi;
    xd=sinb*cosa*cosphi-cosb*sina; yd=sinb*sinphi;
    zc=cosb*cosa+sinb*sina*cosphi;

    if ((zc>0.99999)&&(cosb>=0.0))
    { aalpha=phiab/radian; abeta=agamma=0.0f;
    }
    else if (zc>0.99999)
    { aalpha=-phiab/radian; abeta=agamma=0.0f;
    }
    else if ((zc<-0.99999)&&(cosb>=0.0))
    { aalpha=phiab/radian; abeta=180.0f; agamma=0.0f;
    }
    else if (zc<-0.99999)
    { aalpha=-phiab/radian; abeta=180.0f; agamma=0.0f;
    }
    else
    { if (fabs(xc)+fabs(yc)<1.0e-5) aalpha=0.0f;
      else aalpha=(float)atan2(yc,xc)/radian;
      abeta=(float)acos(zc)/radian;
      if (fabs(xd)+fabs(yd)<1.0e-5) agamma=0.0f;
      else agamma=(float)atan2(yd,xd)/radian;
    }
  }
  else aalpha=abeta=agamma=0.0f;
  if (error==0)
  { while (aalpha>180.0) aalpha-=360.0f;
    while (aalpha<-180.0) aalpha+=360.0f;
    while (agamma>180.0) agamma-=360.0f;
    while (agamma<-180.0) agamma+=360.0f;
    *alpha=aalpha; *beta=abeta; *gamma=agamma;
  }

  return(error);
} //end of Mscdrun::onerotation

int Mscdrun::allrotation()
{ int ia,ib,ic,id,k,eledim;
  float emiter,alpha,beta,gamma;

  if ((error==0)&&(msorder>1)&&(raorder>0)&&((!talpha)||(!tgamma)))
    error=901;
  for (ia=0;(error==0)&&(ia<natoms);++ia)
  { emiter=patom[ia*12+7];
    if ((msorder<2)||(raorder<1)||((msorder==2)&&(emiter==0))) continue;
    for (ib=0;(error==0)&&(ib<natoms);++ib)
    { if (ib==ia) continue;
      for (ic=0;(error==0)&&(ic<natoms);++ic)
      { if (ic==ib) continue;
        id=ia*natoms*natoms+ib*natoms+ic;
        k=tevenadd[id]; eledim=(int)tevenpar[k*10+5];
        if (eledim<2) alpha=gamma=0.0f;
        else error=onerotation(patom+ia*12,patom+ib*12,patom+ic*12,
          &alpha,&beta,&gamma);
        talpha[id]=alpha; tgamma[id]=gamma;
      }
    }
  }

  return(error);
} //end of Mscdrun::allrotation

#ifdef MSCDGPU
#include <stdlib.h>
#include <string.h>
#include "mscdgpu.h"
extern "C" int mscdgpu_set_alnum(const int *alnum,int nkind);
extern "C" int mscdgpu_set_hankb(const Gcplx *h);
extern "C" int mscdgpu_get_dbg(float *out);
extern "C" void mscdgpu_set_dbg(int on);

/* Fase 1 do port (PLANO_CUDA.md). validate!=0: roda a GPU DEPOIS do
   alldblevent da CPU e compara par a par, sem tocar em devenelem -- e' a
   validacao isolada que o plano manda fazer antes de mexer em outro bloco.
   validate==0: a GPU escreve devenelem e a CPU nao roda.

   Ligado por MSCD_GPU=1 (substitui) ou MSCD_GPU=validate (compara). */
int Mscdrun::gpudblevent(float akin,float *xdetec,int validate)
{ static int ready=0;
  static Gcplx *host=NULL;
  static double wmaxabs=0.0,wmaxrel=0.0;
  static long wcount=0,wbig=0;
  int i;

  if (!ready)
  { Gconst k;
    int al[8];
    /* phase.cpp:28 aloca phasec com 61 entradas fixas, independente do lnum
       do arquivo -- e o lnum DIFERE entre as especies (psAg111 e psl9), entao
       o passo tem de ser o da alocacao, nao o do arquivo. */
    int pl=60;
    static Gcplx *pc=NULL;

    for (i=0;i<katoms;++i)
    { phaseshift[i].fsinexpa(akin,0);      /* forca makesinexp */
      al[i]=phaseshift[i].getalnum(akin);
    }
    pc=new Gcplx [katoms*(pl+1)];
    for (i=0;i<katoms;++i)
      memcpy(pc+i*(pl+1),phaseshift[i].gpu_phasec(),
        (pl+1)*sizeof(Gcplx));

    k.patom=patom;             k.natoms=natoms;
    k.devenpar=devenpar;       k.ndbleven=ndbleven;
    k.rotmata=evenmat->gpu_rotmata();
    k.rotmatc=evenmat->gpu_rotmatc();
    k.rlnum=evenmat->gpu_lnum();
    k.lamdum=evenmat->gpu_lamdum();
    k.betanum=evenmat->gpu_betanum();
    k.hankmat_a=(const Gcplx *)hanka->gpu_hankmat();
    k.handata=hanka->gpu_ndata();
    k.halnum=hanka->gpu_lnum();
    k.hacmnum=hanka->gpu_cmnum();
    k.hankarg_b=(const Gcplx *)hankb->gpu_hankarg();
    k.phasec=pc;               k.pclnum=pl;
    k.cexpix=(const Gcplx *)expix->gpu_cexpix();
    k.exndata=expix->gpu_ndata();
    k.exmdata=expix->gpu_mdata();
    k.thermat=vibrate->gpu_thermat();
    k.thernum=vibrate->gpu_thernum();
    k.therstep=vibrate->gpu_therstep();
    k.mweight=vibrate->gpu_mweight();
    k.tdebye=vibrate->gpu_tdebye();
    k.tsample=vibrate->gpu_tsample();
    k.aweight=aweight;         /* mscdrun.cpp:76 aloca 4 */
    k.nkind=(katoms<4)?katoms:4;
    k.radim=radim; k.raorder=raorder;

    if (mscdgpu_setup(&k))
    { std::cerr<<"GPU setup: "<<mscdgpu_lasterror()<<"\n"; return 1; }
    mscdgpu_set_dbg(validate);
    mscdgpu_set_alnum(al,katoms);
    host=new Gcplx [(long)ndbleven*radim];
    ready=1;
    std::cerr<<"GPU: pronta ("<<ndbleven<<" pares, radim="<<radim
      <<", alnum="<<al[0]<<", rlnum="<<k.rlnum<<", lamdum="<<k.lamdum
      <<", betanum="<<k.betanum<<", halnum="<<k.halnum
      <<", handata="<<k.handata<<", cmnum="<<k.hacmnum
      <<", nkind="<<k.nkind<<")\n";
  }

  /* hankb: alldblevent chama sempre com vkb=0, mas o cache e' keyed no
     argumento e outras funcoes do laco o mexem. Em vez de reproduzir esse
     automato no device, forca-se o mesmo estado e tira-se a foto. */
  hankb->fhankelfaca(0,0,0.0f);
  mscdgpu_set_hankb((const Gcplx *)hankb->gpu_hankarg());

  if (mscdgpu_alldblevent(akin,xdetec,meanpath->finvpath(akin),host))
  { std::cerr<<"GPU: "<<mscdgpu_lasterror()<<"\n"; return 1; }

  if (!validate)
  { memcpy(devenelem,host,(long)ndbleven*radim*sizeof(Gcplx));
    return error;
  }

  /* Diagnostico do primeiro ponto: para os piores pares, recalcula na CPU as
     grandezas que o kernel derivou e imprime lado a lado. beta e row sao
     INDICES -- se divergirem, o erro nao e' de arredondamento. */
  { static int once=0;
    if (validate&&!once)
    { once=1;
      float *dbg=new float [(long)ndbleven*4];
      if (!mscdgpu_get_dbg(dbg))
      { const Gcplx *c=(const Gcplx *)devenelem;
        int shown=0,j;
        for (j=0;(j<ndbleven)&&(shown<12);++j)
        { double worst=0.0; int e;
          for (e=0;e<radim;++e)
          { long ix=(long)j*radim+e;
            double dr=(double)host[ix].re-c[ix].re;
            double di=(double)host[ix].im-c[ix].im;
            double da=sqrt(dr*dr+di*di);
            double m=sqrt((double)c[ix].re*c[ix].re+(double)c[ix].im*c[ix].im);
            if ((m>1.0e-12)&&(da/m>worst)) worst=da/m;
          }
          if (worst>1.0e-3)
          { int ia=(int)devenpar[j*7+5], ib=(int)devenpar[j*7+6];
            float ydetec[3],alpha,beta,gamma,cxa;
            ydetec[0]=patom[ib*12]+xdetec[0];
            ydetec[1]=patom[ib*12+1]+xdetec[1];
            ydetec[2]=patom[ib*12+2]+xdetec[2];
            onerotation(patom+ia*12,patom+ib*12,ydetec,&alpha,&beta,&gamma);
            cxa=beta;
            beta=(float)floor((beta*10.0+0.5)/10.0);
            std::cerr<<"DBG j="<<j<<" rel="<<worst
              <<" | CPU beta_bruto="<<cxa<<" beta="<<beta
              <<" gamma="<<gamma
              <<" | GPU beta="<<dbg[j*4]<<" gamma="<<dbg[j*4+1]
              <<" row="<<dbg[j*4+3]
              <<" | dbeta="<<(dbg[j*4]-beta)
              <<" dgamma="<<(dbg[j*4+1]-gamma)<<"\n";
            ++shown;
          }
        }
      }
      delete [] dbg;
    }
  }

  { long n=(long)ndbleven*radim,idx;
    const Gcplx *c=(const Gcplx *)devenelem;
    for (idx=0;idx<n;++idx)
    { double dr=(double)host[idx].re-c[idx].re;
      double di=(double)host[idx].im-c[idx].im;
      double da=sqrt(dr*dr+di*di);
      double m=sqrt((double)c[idx].re*c[idx].re+(double)c[idx].im*c[idx].im);
      ++wcount;
      if (da>wmaxabs) wmaxabs=da;
      if (m>1.0e-12)
      { double rel=da/m;
        if (rel>wmaxrel) wmaxrel=rel;
        if (rel>1.0e-3) ++wbig;
      }
    }
    std::cerr<<"GPUVAL n="<<wcount<<" maxabs="<<wmaxabs
      <<" maxrel="<<wmaxrel<<" rel>1e-3="<<wbig<<"\n";
  }
  return error;
}
#endif

#ifdef BETATEST
extern int betatest_indbl;
#endif

int Mscdrun::alldblevent(float akin,float *xdetec)
{ int ia,ib,j,akind;
  float xa,xb,xc,xd,ka,cosbeta,vka,vkb,alpha,beta,gamma;
  float ydetec[3];
  Fcomplex cxa,cxb,cxc,cxd,cxe,cxf,cxg,cxh,cxi,cvalue;
  const float radian=(float)(3.14159265/180.0);

#ifdef BETATEST
  betatest_indbl=1;
#endif
  if ((error==0.0)&&(msorder>0)&&(!devenelem)) error=901;
  for (j=0;(error==0)&&(j<ndbleven);++j)
  { ia=(int)devenpar[j*7+5]; ib=(int)devenpar[j*7+6];
    ydetec[0]=patom[ib*12]+xdetec[0];
    ydetec[1]=patom[ib*12+1]+xdetec[1];
    ydetec[2]=patom[ib*12+2]+xdetec[2];
    error=onerotation(patom+ia*12,patom+ib*12,ydetec,
      &alpha,&beta,&gamma);
    if (error!=0) break;
    beta=(float)floor((beta*10.0+0.5)/10.0);
    xa=devenpar[j*7+3]; akind=(int)devenpar[j*7+4]; ka=akin*xa;
    cosbeta=(float)cos(beta*radian);
    xc=meanpath->finvpath(akin);
    xd=vibrate->fvibmsrd(xa,aweight[akind-1]);
    xb=(float)exp(-0.5*xa*xc-akin*akin*xd*(1.0-cosbeta))/ka;
    cxa=xb*expix->fexpix(ka/radian);
    if (raorder>0)
    { cxb=xb*expix->fexpix(ka/radian-gamma);
      cxc=xb*expix->fexpix(ka/radian+gamma);
    }
    if (raorder>1)
    { cxd=xb*expix->fexpix(ka/radian-gamma-gamma);
      cxe=xb*expix->fexpix(ka/radian+gamma+gamma);
    }
    if (raorder>2)
    { cxf=xb*expix->fexpix(ka/radian-gamma-gamma-gamma);
      cxg=xb*expix->fexpix(ka/radian+gamma+gamma+gamma);
    }
    if (raorder>3)
    { cxh=xb*expix->fexpix(ka/radian-gamma-gamma-gamma-gamma);
      cxi=xb*expix->fexpix(ka/radian+gamma+gamma+gamma+gamma);
    }
    if (raorder<0) vka=0.0f; else vka=1.0f/ka;
    vkb=0.0f;
    if (ATA==1)
    { cvalue=ATAevenelem(akind,0,0,0,0,akin,vka,vkb,beta,ATAweight1,                                 ATAweight2);
      devenelem[j*radim]=cxa*cvalue;
    }
    else
    {
      cvalue=evenelem(akind,0,0,0,0,akin,vka,vkb,beta);
      devenelem[j*radim]=cxa*cvalue;
    }
    if (raorder>0)
    { if (ATA==1)
      { cvalue=ATAevenelem(akind,1,0,0,0,akin,vka,vkb,beta,ATAweight1,                                 ATAweight2);
        devenelem[j*radim+1]=cxb*cvalue;
        devenelem[j*radim+2]=-cxc*cvalue;
      }
      else
      {
        cvalue=evenelem(akind,1,0,0,0,akin,vka,vkb,beta);
        devenelem[j*radim+1]=cxb*cvalue;
        devenelem[j*radim+2]=-cxc*cvalue;
      }
    }
    if (raorder>1)
    { if (ATA==1)
      {  cvalue=ATAevenelem(akind,0,1,0,0,akin,vka,vkb,beta,ATAweight1,                                 ATAweight2);
         devenelem[j*radim+3]=cxa*cvalue;
         cvalue=evenelem(akind,2,0,0,0,akin,vka,vkb,beta);
         devenelem[j*radim+4]=cxd*cvalue;
         devenelem[j*radim+5]=cxe*cvalue;
      }
      else
      {
         cvalue=evenelem(akind,0,1,0,0,akin,vka,vkb,beta);
         devenelem[j*radim+3]=cxa*cvalue;
         cvalue=evenelem(akind,2,0,0,0,akin,vka,vkb,beta);
         devenelem[j*radim+4]=cxd*cvalue;
         devenelem[j*radim+5]=cxe*cvalue;
      }
    }
    if (raorder>2)
    { if (ATA==1)
      {  cvalue=ATAevenelem(akind,1,1,0,0,akin,vka,vkb,beta,ATAweight1,                                 ATAweight2);
         devenelem[j*radim+6]=cxb*cvalue;
         devenelem[j*radim+7]=-cxc*cvalue;
         cvalue=ATAevenelem(akind,3,0,0,0,akin,vka,vkb,beta,ATAweight1,                                 ATAweight2);
         devenelem[j*radim+8]=cxf*cvalue;
         devenelem[j*radim+9]=-cxg*cvalue;
      }
      else
      {
         cvalue=evenelem(akind,1,1,0,0,akin,vka,vkb,beta);
         devenelem[j*radim+6]=cxb*cvalue;
         devenelem[j*radim+7]=-cxc*cvalue;
         cvalue=evenelem(akind,3,0,0,0,akin,vka,vkb,beta);
         devenelem[j*radim+8]=cxf*cvalue;
         devenelem[j*radim+9]=-cxg*cvalue;
      }
    }
    if (raorder>3)
    { if (ATA==1)
      { cvalue=ATAevenelem(akind,0,2,0,0,akin,vka,vkb,beta,ATAweight1,                                 ATAweight2);
        devenelem[j*radim+10]=cxa*cvalue;
        cvalue=ATAevenelem(akind,2,1,0,0,akin,vka,vkb,beta,ATAweight1,                                 ATAweight2);
        devenelem[j*radim+11]=cxd*cvalue;
        devenelem[j*radim+12]=cxe*cvalue;
        cvalue=ATAevenelem(akind,4,0,0,0,akin,vka,vkb,beta,ATAweight1,                                 ATAweight2);
        devenelem[j*radim+13]=cxh*cvalue;
        devenelem[j*radim+14]=cxi*cvalue;
      }
      else
      {
        cvalue=evenelem(akind,0,2,0,0,akin,vka,vkb,beta);
        devenelem[j*radim+10]=cxa*cvalue;
        cvalue=evenelem(akind,2,1,0,0,akin,vka,vkb,beta);
        devenelem[j*radim+11]=cxd*cvalue;
        devenelem[j*radim+12]=cxe*cvalue;
        cvalue=evenelem(akind,4,0,0,0,akin,vka,vkb,beta);
        devenelem[j*radim+13]=cxh*cvalue;
        devenelem[j*radim+14]=cxi*cvalue;
      }
    }
  }
#ifdef BETATEST
  betatest_indbl=0;
#endif
  return(error);
} //end of Mscdrun::alldblevent

int Mscdrun::allevendetec(float akin,float *xdetec)
{ int ia,ib,j,k;
  float xa,xc,xb,yb,zb,xd,yd,zd,cosd,ka,emiter;
  Fcomplex cxa;
  const float radian=(float)(3.14159265/180.0);

  xd=xdetec[0]; yd=xdetec[1]; zd=cosd=xdetec[2];
  if ((error==0)&&(msorder>0)&&((!devendetec)||(cosd<1.0e-5)||
    (cosd>1.001))) error=901;
  else if ((error==0)&&(msorder>0))
  { for (ib=0;ib<natoms;++ib)
    { xb=patom[ib*12]; yb=patom[ib*12+1]; zb=patom[ib*12+2];
      ka=-akin*(xb*xd+yb*yd+zb*zd); xc=meanpath->finvpath(akin);
      xa=(float)exp(0.5*xc*zb/cosd); cxa=xa*expix->fexpix(ka/radian);
      for (ia=0;ia<natoms;++ia)
      { if (ia==ib) continue;
        emiter=patom[ia*12+7];
        if ((msorder==1)&&(emiter==0)) continue;
        k=devenadd[ia*natoms+ib];
        for (j=0;j<radim;++j) devendetec[ia*natoms*radim+ib*radim+j]=
          devenelem[k*radim+j]*cxa;
      }
    }
  }
  return(error);
} //end of Mscdrun::allevendetec

int Mscdrun::onevenemit(int ia,int ib,int alf,int am,float akin,
  float *xdetec,float *polaron,Fcomplex *aemitelem)
{ float xa,xb,yb,zb,lengb,kb,vkb,alpha,beta,gamma,emiter;
  float xatoma[3],xatomb[3],xatomc[3];
  Fcomplex cxa,cxb,cxc,cxd,cxe,cxf,cxg,cxh,cxi,cvalue;

  emiter=patom[ia*12+7];
  if ((error==0)&&(msorder>0)&&(emiter!=0))
  { xb=patom[ib*12]-patom[ia*12]; yb=patom[ib*12+1]-patom[ia*12+1];
    zb=patom[ib*12+2]-patom[ia*12+2];
    lengb=(float)sqrt(xb*xb+yb*yb+zb*zb);
    kb=akin*lengb;
    if ((raorder<0)||(kb<1.0e-5)) vkb=0.0f; else vkb=1.0f/kb;
    if (lengb<1.0e-5)
    { alpha=0.0f; beta=0.0f; gamma=180.0f;
    }
    else
    { xatomc[0]=patom[ib*12]+xdetec[0];
      xatomc[1]=patom[ib*12+1]+xdetec[1];
      xatomc[2]=patom[ib*12+2]+xdetec[2];
      error=onerotation(patom+ia*12,patom+ib*12,xatomc,
        &alpha,&beta,&gamma);
      xa=gamma;
      xatoma[0]=xatoma[1]=xatoma[2]=0.0f;
      Cartesia pbond(xb,yb,zb);
      pbond=pbond.euler(0.0f,polaron[3],180.0f-polaron[4]);
      pbond.getcoordinates(xatomb,xatomb+1,xatomb+2);
      pbond.loadcoordinates(xb+xdetec[0],yb+xdetec[1],zb+xdetec[2]);
      pbond=pbond.euler(0.0f,polaron[3],180.0f-polaron[4]);
      pbond.getcoordinates(xatomc,xatomc+1,xatomc+2);
      error=onerotation(xatoma,xatomb,xatomc,&alpha,&beta,&gamma);
      alpha=gamma-xa;

      pbond.loadcoordinates(xb,yb,zb);
      pbond=pbond.euler(0.0f,polaron[3],180.0f-polaron[4]);
      beta=pbond.theta(); gamma=180.0f-pbond.phi();
    }
    cxa=expix->fexpix(-gamma*am);
    cxb=expix->fexpix(-alpha-gamma*am);
    cxc=expix->fexpix(alpha-gamma*am);
    cxd=expix->fexpix(-alpha-alpha-gamma*am);
    cxe=expix->fexpix(alpha+alpha-gamma*am);
    cxf=expix->fexpix(-alpha-alpha-alpha-gamma*am);
    cxg=expix->fexpix(alpha+alpha+alpha-gamma*am);
    cxh=expix->fexpix(-alpha-alpha-alpha-alpha-gamma*am);
    cxi=expix->fexpix(alpha+alpha+alpha+alpha-gamma*am);
    cvalue=evenbelem(alf,am,0,0,0,vkb,beta);
    aemitelem[0]=cxa*cvalue;
    if ((ia!=ib)&&(raorder>0))
    { cvalue=evenbelem(alf,am,0,1,0,vkb,beta);
      aemitelem[1]=cxb*cvalue;
      cvalue=evenbelem(alf,am,0,-1,0,vkb,beta);
      aemitelem[2]=cxc*cvalue;
    }
    if ((ia!=ib)&&(raorder>1))
    { cvalue=evenbelem(alf,am,0,0,1,vkb,beta);
      aemitelem[3]=cxa*cvalue;
      cvalue=evenbelem(alf,am,0,2,0,vkb,beta);
      aemitelem[4]=cxd*cvalue;
      cvalue=evenbelem(alf,am,0,-2,0,vkb,beta);
      aemitelem[5]=cxe*cvalue;
    }
    if ((ia!=ib)&&(raorder>2))
    { cvalue=evenbelem(alf,am,0,1,1,vkb,beta);
      aemitelem[6]=cxb*cvalue;
      cvalue=evenbelem(alf,am,0,-1,1,vkb,beta);
      aemitelem[7]=cxc*cvalue;
      cvalue=evenbelem(alf,am,0,3,0,vkb,beta);
      aemitelem[8]=cxf*cvalue;
      cvalue=evenbelem(alf,am,0,-3,0,vkb,beta);
      aemitelem[9]=cxg*cvalue;
    }
    if ((ia!=ib)&&(raorder>3))
    { cvalue=evenbelem(alf,am,0,0,2,vkb,beta);
      aemitelem[10]=cxa*cvalue;
      cvalue=evenbelem(alf,am,0,2,1,vkb,beta);
      aemitelem[11]=cxd*cvalue;
      cvalue=evenbelem(alf,am,0,-2,1,vkb,beta);
      aemitelem[12]=cxe*cvalue;
      cvalue=evenbelem(alf,am,0,4,0,vkb,beta);
      aemitelem[13]=cxh*cvalue;
      cvalue=evenbelem(alf,am,0,-4,0,vkb,beta);
      aemitelem[14]=cxi*cvalue;
    }
  }
  return(error);
} //end of Mscdrun::onevenemit

Fcomplex Mscdrun::onemidetec(float akin,int ie,int alf,int am,
  float *xdetec,float *polaron)
{ float xa,ya,za,xd,yd,zd,cosd,ka,xe,xf,atheta,aphi,gamma,emiter;
  Fcomplex cxa,cvalue;
  const float radian=(float)(3.14159265/180.0);

  xd=xdetec[0]; yd=xdetec[1]; zd=cosd=xdetec[2];
  emiter=patom[ie*12+7];
  if ((error==0)&&(emiter!=0)&&(cosd>1.0e-5))
  { xa=patom[ie*12]; ya=patom[ie*12+1]; za=patom[ie*12+2];
    Cartesia patomb(xd,yd,zd);
    Cartesia patomc=patomb.euler(0.0f,polaron[3],180.0f-polaron[4]);
    atheta=patomc.theta(); aphi=patomc.phi();
    xe=meanpath->finvpath(akin);
    xf=(float)exp(0.5*xe*za/cosd); ka=-akin*(xa*xd+ya*yd+za*zd);
    gamma=180.0f-aphi;
    cxa=xf*expix->fexpix(-gamma*am+ka/radian);
    cvalue=evenbelem(alf,am,0,0,0,0.0f,atheta);
    cvalue=cxa*cvalue;
  }
  else cvalue=0.0f;
  return(cvalue);
} //end of Mscdrun::onemidetec

