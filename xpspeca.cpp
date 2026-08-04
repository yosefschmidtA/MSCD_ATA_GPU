#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

#include "curvefit.h"
#include "userutil.h"
#include "xpspec.h"

/*
----------------------------------------------------------------------
emission  0        1       2     3
          ephoton  theta   phi   reserved

          4                    5               6            7
          peak-binding-energy  peak-intensity  peak-height  peak-width

          8                    9               10           11
          peak-binding-energy  peak-intensity  peak-height  peak-width

          12                   13              14           15
          peak-binding-energy  peak-intensity  peak-height  peak-width

----------------------------------------------------------------------
*/
int Spectrum::peakfit(int inpeak,float *epeak,int backassociate,
  int peakassociate)
{ int i,j,ka,kb,kc,kd,ke,ndata,fitmath,trymax,trynum,jobcount;
  float xa,xb,xc,xd,xe,xf,xg,xh,xs,tolerance,reliable,ephoton;
  float *afit,*safit,*abest,*ybest;

  afit=safit=abest=ybest=NULL; jobcount=0;
  Textout conout;
  if (error==0)
  { if (inpeak<0) npeak=0;
    else if (inpeak>3) npeak=3;
    else npeak=inpeak;
    nfit=3+npeak*3;
    afit=new float [ncurve*nfit]; safit=new float [ncurve*nfit];
    abest=new float [ncurve*nfit];
    ybest=new float [1024];
    if ((!afit)||(!safit)||(!abest)||(!ybest)) error=102;
  }

  if ((error==0)&&(npeak>0))
  { numint=(npeak+1)*4;
    if (emission) delete [] emission;
    emission=new float [ncurve*numint];
    if (!emission) error=102;
  }

  //first pass of fitting
  fitmath=2; trymax=2000; tolerance=1.0e-10f;
  if (error==0)
  { for (i=0;i<ncurve;++i)
    { for (j=0;j<nfit;++j) afit[i*nfit+j]=0.0f;
      ephoton=curvepar[i*8];
      ka=(int)curvepar[i*8+3];
      kb=(int)curvepar[i*8+4];
      ndata=kb-ka;
      if (ndata>256) ndata=256;

      if (ndata<1) kd=1;
      else if (ndata<5) kd=ndata;
      else kd=5;

      xa=xb=xc=xs=0.0f;
      for (j=0;j<kd;++j) xa+=xpscount[ka+j];
      xa/=(float)kd;
      for (j=ndata-1;j>=ndata-kd;--j) xb+=xpscount[ka+j];
      xb/=(float)kd;
      for (j=1;j<ndata;++j)
        xs+=xpscount[ka+j]*(ebind[ka+j]-ebind[ka+j-1]);
      xc=xs-xa*(ebind[ka+ndata-1]-ebind[ka]);
      if ((xc>-1.0e-20)&&(xc<1.0e-20)) xc=1.0f;
      if (nfit>0) afit[i*nfit]=xa;
      if (nfit>1) afit[i*nfit+1]=(xb-xa)/xc;
      if (nfit>2) afit[i*nfit+2]=0.0f;

      xc=ebind[ka+ndata-1]-ebind[ka];
      if (xc<1.0e-20) xc=1.0f;
      xc=(xb-xa)/xc;
      xd=0.0f;
      for (kc=ka+2;kc<kb-2;++kc)
      { xg=xa+xc*(ebind[kc]-ebind[ka]);
        xf=(xpscount[kc-2]+xpscount[kc-1]+xpscount[kc]+xpscount[kc+1]+
          xpscount[kc+2])*0.2f-xg;
        for (j=1;j<npeak;++j)
        { xe=epeak[j]-epeak[0]+ebind[kc];
          kd=ka;
          while ((kd<kb)&&(ebind[kd]<xe)) ++kd;
          if ((kd>ka+2)&&(kd<kb-2))
          { xg=xa+xc*(ebind[kd]-ebind[ka]);
            xf+=(xpscount[kd-2]+xpscount[kd-1]+xpscount[kd]+
              xpscount[kd+1]+xpscount[kd+2])*0.2f-xg;
          }
          else xf=0.0f;
        }
        if (xd<xf)
        { xd=xf; ke=kc;
        }
      }
      for (j=0;j<npeak;++j)
      { xf=ebind[ke]+epeak[j]-epeak[0];
        kd=ka;
        while ((xf>ebind[kd])&&(kd<kb-1)) ++kd;
        xf-=ebind[ka];
        xd=xpscount[kd]-xa-xc*xf;
        if (xd<1.0e-10) xd=1.0e-10f;
        xe=ebind[ka+ndata-1]-ebind[ka];
        if (xe>1.0e-10) xe=50.0f/xe/xe;
        xe*=(float)(npeak*npeak);
        afit[i*nfit+j*3+3]=xd; afit[i*nfit+j*3+4]=xe;
        afit[i*nfit+j*3+5]=xf;
        if (j>0) afit[i*nfit+j*3+5]-=afit[i*nfit+5];
      }
      if (npeak>1)
      { xc=xd=xe=0.0f;
        if (npeak>0) xc=afit[i*nfit+3];
        if (npeak>1) xd=afit[i*nfit+6];
        if (npeak>2) xe=afit[i*nfit+9];
        if ((xe>xc)&&(xe>xd)) kd=2;
        else if ((xd>xc)&&(xd>xe)) kd=1;
        else kd=0;
        if (kd>npeak-1) kd=npeak-1;
        if (kd<0) kd=0;
        xe=afit[i*nfit+kd*3+3];
        for (j=0;j<npeak;++j)
        { if (j==kd) continue;
          else if (j==0) xc=afit[i*nfit+kd*3+5];
          else if (kd==0) xc=afit[i*nfit+j*3+5];
          else xc=afit[i*nfit+j*3+5]-afit[i*nfit+kd*3+5];
          xd=afit[i*nfit+j*3+4]*xc*xc;
          if ((xd>0.0)&&(xd<60.0))
            xe+=afit[i*nfit+j*3]*(float)exp(-xc);
        }
        if (xe>1.0e-10) xe=afit[i*nfit+kd*3+3]/xe;
        else xe=1.0f;
        for (j=0;j<npeak;++j) afit[i*nfit+j*3+3]*=xe;
      }
      xa=afit[i*nfit]; xb=0.0f;
      for (j=0;j<npeak;++j)
        if (xb<afit[i*nfit+j*3+3]) xb=afit[i*nfit+j*3+3];

      //first pass of fitting
      for (j=0;j<nfit;++j)
      { if (nfit<4) safit[i*nfit+j]=0.01f;
        else if (j<3) safit[i*nfit+j]=0.0f;
        else if ((j%3)==0) safit[i*nfit+j]=0.01f;
        else safit[i*nfit+j]=0.01f;
      }

      if ((error==0)&&(nfit>1))
      { Curvefit xpscurve;
        xpscurve.init(ndata,nfit);
        error=xpscurve.loadcurve(fitmath,trymax,tolerance,ebind+ka,
          xpscount+ka,afit+i*nfit,safit+i*nfit,fitxpspeak);
        if (error==0)
          error=xpscurve.dofit(ybest,abest+i*nfit,&reliable,&trynum);
      }
      if (error==0)
      { curvepar[i*8+5]=reliable;
        for (j=0;j<nfit;++j)
        { afit[i*nfit+j]=abest[i*nfit+j];
          if (((j==1)||(j==2))&&(abest[i*nfit]>1.0e-20))
            abest[i*nfit+j]/=abest[i*nfit];
          else if (j==5) abest[i*nfit+j]+=ebind[ka];
        }
      }
      kd=(i+1)*21/ncurve-jobcount;
      if ((error==0)&&(kd>0))
      { conout.charfill('-',kd); conout.flush(); jobcount+=kd;
      }
    }
  }

  //second pass of fitting
  fitmath=2; trymax=2000; tolerance=1.0e-10f;
  if ((ncurve>3)&&(error==0))
  { for (j=1;(j<3);++j)
    { xa=abest[j]; xb=abest[nfit+j];
      for (i=2;i<ncurve-2;++i)
      { xc=abest[i*nfit+j];
        xd=xa+xb+xc+abest[(i+1)*nfit+j]+abest[(i+2)*nfit+j];
        xe=xa;
        if (xe>xb) xe=xb;
        if (xe>xc) xe=xc;
        if (xe>abest[(i+1)*nfit+j]) xe=abest[(i+1)*nfit+j];
        if (xe>abest[(i+2)*nfit+j]) xe=abest[(i+2)*nfit+j];
        xd-=xe;
        xe=xa;
        if (xe<xb) xe=xb;
        if (xe<xc) xe=xc;
        if (xe<abest[(i+1)*nfit+j]) xe=abest[(i+1)*nfit+j];
        if (xe<abest[(i+2)*nfit+j]) xe=abest[(i+2)*nfit+j];
        xd-=xe;
        xa=xb; xb=xc;
        abest[i*nfit+j]=xd/3.0f;
      }
      abest[nfit+j]=(abest[j]+abest[nfit+j]+abest[2*nfit+j])/3.0f;
      abest[(ncurve-2)*nfit+j]=(abest[(ncurve-3)*nfit+j]+
        abest[(ncurve-2)*nfit+j]+abest[(ncurve-1)*nfit+j])/3.0f;
    }

    for (j=3;j<nfit;++j)
    { xa=xb=abest[j]; xc=0.0f;
      for (i=0;i<ncurve;++i)
      { if (xa>abest[i*nfit+j]) xa=abest[i*nfit+j];
        if (xb<abest[i*nfit+j]) xb=abest[i*nfit+j];
        xc+=abest[i*nfit+j];
      }
      if (ncurve>0) xc/=(float)ncurve;
      xa=xc-xa; xb=xb-xc;
      if (xa<xb) xa=xb;
      xg=xc-0.4f*xa; xh=xc+0.4f*xa;
      xe=0.0f; kd=0;
      for (i=0;i<ncurve;++i)
      { xf=abest[i*nfit+j];
        if ((xf>xg)&&(xf<xh))
        { ++kd; xe+=xf;
        }
      }
      if (kd>0) xe=xe/(float)kd;
      xg+=xe-xc; xh+=xe-xc;
      kd=0; xa=xb=xc=xd=0.0f;
      for (i=0;i<ncurve;++i)
      { xe=curvepar[i*8]; xf=abest[i*nfit+j];
        if ((xf>xg)&&(xf<xh))
        { ++kd;
          xa+=xe; xb+=xe*xe;
          xc+=xf; xd+=xe*xf;
        }
      }
      xe=xb*xc-xa*xd; xf=(float)kd*xd-xa*xc;
      xg=(float)kd*xb-xa*xa;
      if ((kd>1)&&(xg>1.0e-20))
      { xe/=xg; xf/=xg;
        for (i=0;i<ncurve;++i)
        { xg=curvepar[i*8]; xh=xe+xf*xg;
          if ((j%3)!=0) abest[i*nfit+j]=xh;
          else if (abest[i*nfit+j]<0.0) abest[i*nfit+j]=xh;
        }
      }
    }
  }
  if (error==0)
  { for (i=0;i<ncurve;++i)
    { ephoton=curvepar[i*8];
      ka=(int)curvepar[i*8+3];
      kb=(int)curvepar[i*8+4];
      ndata=kb-ka;
      if (ndata>256) ndata=256;
      if (npeak>1)
      { kc=1; xb=abest[i*nfit+8];
        for (j=2;j<npeak;++j)
        { if (xb>abest[i*nfit+j*3+5])
          { xb=abest[i*nfit+j*3+5];
            kc=j;
          }
        }
        if (xb>0.0) kc=0;
      }
      else kc=0;
      if (npeak>0)
      { xa=abest[i*nfit+kc*3+5];
        if (kc>0) xa+=abest[i*nfit+5];
        xb=abest[i*nfit+kc*3+4];
        if (xb>1.0e-5) xb=(float)sqrt(2.0/xb);
      }
      else xa=xb=0.0f;
      xe=xa-2.0f*xb;
      kd=ka;
      while ((kd<kb)&&(ebind[kd]<xe)) ++kd;
      kd=kd-ka-1;
      if (kd<5) kd=5;
      xa=0.0f;
      for (j=0;j<kd;++j) xa+=xpscount[ka+j];
      xa/=(float)kd;
      afit[i*nfit]=xa;
      for (j=3;j<nfit;++j)
      { if ((j==1)||(j==2))
          afit[i*nfit+j]=abest[i*nfit]*abest[i*nfit+j];
        else if (j==5)
          afit[i*nfit+j]=abest[i*nfit+j]-ebind[ka];
        else
          afit[i*nfit+j]=abest[i*nfit+j];
      }
      xa=xb=xc=xd=0.0f;
      for (j=1;j<npeak;++j)
      { if (xa>afit[i*nfit+j*3+5]) xa=afit[i*nfit+j*3+5];
        if (xb<afit[i*nfit+j*3+5]) xb=afit[i*nfit+j*3+5];
      }
      for (j=0;j<npeak;++j)
      { xe=afit[i*nfit+j*3+4];
        if (xe>1.0e-5) xe=(float)sqrt(2.0/xe);
        if (xc<xe) xc=xe;
        xd+=xe;
      }
      if ((xa>0.0)&&(xb>0.0)) xa=xb;
      else if ((xa<0.0)&&(xb<0.0)) xa=-xa;
      else xa=xb-xa;
      xa+=xd; xc*=2.0f;
      if (xa<xc) xa=xc;
      xa*=2.5f;
      if (xa<8.0) xa=8.0f;
      xb=ebind[ka+ndata-1]-ebind[ka];

      //second pass of fitting
      for (j=0;j<nfit;++j)
      { if (j<1) safit[i*nfit+j]=0.0f;
        else if (j==1) safit[i*nfit+j]=0.01f;
        else if ((j==2)&&(xa<xb)) safit[i*nfit+j]=0.0f;
        else if (j==2) safit[i*nfit+j]=0.0f;
        else if ((j%3)!=2) safit[i*nfit+j]=0.01f;
        else if (j!=5) safit[i*nfit+j]=0.01f;
        else if (afit[i*nfit+j]<0.1) safit[i*nfit+j]=1.0f;
        else safit[i*nfit+j]=0.01f/afit[i*nfit+j];
      }
    }
  }

  if (error==0)
  { for (i=0;(error==0)&&(i<ncurve);++i)
    { ephoton=curvepar[i*8];
      ka=(int)curvepar[i*8+3];
      kb=(int)curvepar[i*8+4];
      ndata=kb-ka;
      if (ndata>256) ndata=256;

      if ((error==0)&&(nfit>1))
      { Curvefit xpscurve;
        xpscurve.init(ndata,nfit);
        error=xpscurve.loadcurve(fitmath,trymax,tolerance,ebind+ka,
          xpscount+ka,afit+i*nfit,safit+i*nfit,fitxpspeak);
        if (error==0)
          error=xpscurve.dofit(ybest,abest+i*nfit,&reliable,&trynum);
      }
      if (error==0)
      { curvepar[i*8+5]=reliable;
        for (j=0;j<nfit;++j)
        { afit[i*nfit+j]=abest[i*nfit+j];
          if (((j==1)||(j==2))&&(abest[i*nfit]>1.0e-20))
            abest[i*nfit+j]/=abest[i*nfit];
          else if (j==5) abest[i*nfit+j]+=ebind[ka];
        }
      }
      kd=21+(i+1)*21/ncurve-jobcount;
      if ((error==0)&&(kd>0))
      { conout.charfill('-',kd); conout.flush(); jobcount+=kd;
      }
    }
  }

  //third pass of fitting
  fitmath=2; trymax=2000; tolerance=1.0e-10f;
  if ((ncurve>3)&&(error==0))
  { for (j=1;j<3;++j)
    { xa=abest[j]; xb=abest[nfit+j];
      for (i=2;i<ncurve-2;++i)
      { xc=abest[i*nfit+j];
        xd=xa+xb+xc+abest[(i+1)*nfit+j]+abest[(i+2)*nfit+j];
        xe=xa;
        if (xe>xb) xe=xb;
        if (xe>xc) xe=xc;
        if (xe>abest[(i+1)*nfit+j]) xe=abest[(i+1)*nfit+j];
        if (xe>abest[(i+2)*nfit+j]) xe=abest[(i+2)*nfit+j];
        xd-=xe;
        xe=xa;
        if (xe<xb) xe=xb;
        if (xe<xc) xe=xc;
        if (xe<abest[(i+1)*nfit+j]) xe=abest[(i+1)*nfit+j];
        if (xe<abest[(i+2)*nfit+j]) xe=abest[(i+2)*nfit+j];
        xd-=xe;
        xa=xb; xb=xc;
        abest[i*nfit+j]=xd/3.0f;
      }
      abest[nfit+j]=(abest[j]+abest[nfit+j]+abest[2*nfit+j])/3.0f;
      abest[(ncurve-2)*nfit+j]=(abest[(ncurve-3)*nfit+j]+
        abest[(ncurve-2)*nfit+j]+abest[(ncurve-1)*nfit+j])/3.0f;
    }

    for (j=3;j<nfit;++j)
    { xa=xb=abest[j]; xc=0.0f;
      for (i=0;i<ncurve;++i)
      { if (xa>abest[i*nfit+j]) xa=abest[i*nfit+j];
        if (xb<abest[i*nfit+j]) xb=abest[i*nfit+j];
        xc+=abest[i*nfit+j];
      }
      if (ncurve>0) xc/=(float)ncurve;
      xa=xc-xa; xb=xb-xc;
      if (xa<xb) xa=xb;
      xg=xc-0.4f*xa; xh=xc+0.4f*xa;
      xe=0.0f; kd=0;
      for (i=0;i<ncurve;++i)
      { xf=abest[i*nfit+j];
        if ((xf>xg)&&(xf<xh))
        { ++kd; xe+=xf;
        }
      }
      if (kd>0) xe=xe/(float)kd;
      xg+=xe-xc; xh+=xe-xc;
      kd=0; xa=xb=xc=xd=0.0f;
      for (i=0;i<ncurve;++i)
      { xe=curvepar[i*8]; xf=abest[i*nfit+j];
        if ((xf>xg)&&(xf<xh))
        { ++kd;
          xa+=xe; xb+=xe*xe;
          xc+=xf; xd+=xe*xf;
        }
      }
      xe=xb*xc-xa*xd; xf=(float)kd*xd-xa*xc;
      xg=(float)kd*xb-xa*xa;
      if ((kd>1)&&(xg>1.0e-20))
      { xe/=xg; xf/=xg;
        for (i=0;i<ncurve;++i)
        { xg=curvepar[i*8]; xh=xe+xf*xg;
          if ((j%3)!=0) abest[i*nfit+j]=xh;
          else if (abest[i*nfit+j]<0.0) abest[i*nfit+j]=xh;
        }
      }
    }
  }
  if (error==0)
  { for (i=0;i<ncurve;++i)
    { ephoton=curvepar[i*8];
      ka=(int)curvepar[i*8+3];
      kb=(int)curvepar[i*8+4];
      ndata=kb-ka;
      if (ndata>256) ndata=256;
      if (npeak>1)
      { kc=1; xb=abest[i*nfit+8];
        for (j=2;j<npeak;++j)
        { if (xb>abest[i*nfit+j*3+5])
          { xb=abest[i*nfit+j*3+5];
            kc=j;
          }
        }
        if (xb>0.0) kc=0;
      }
      else kc=0;
      if (npeak>0)
      { xa=abest[i*nfit+kc*3+5];
        if (kc>0) xa+=abest[i*nfit+5];
        xb=abest[i*nfit+kc*3+4];
        if (xb>1.0e-5) xb=(float)sqrt(2.0/xb);
      }
      else xa=xb=0.0f;
      xe=xa-2.0f*xb;
      kd=ka;
      while ((kd<kb)&&(ebind[kd]<xe)) ++kd;
      kd=kd-ka-1;
      if (kd<5) kd=5;
      xa=0.0f;
      for (j=0;j<kd;++j) xa+=xpscount[ka+j];
      xa/=(float)kd;
      afit[i*nfit]=abest[i*nfit];
      for (j=1;j<nfit;++j)
      { if ((j==1)||(j==2))
          afit[i*nfit+j]=abest[i*nfit]*abest[i*nfit+j];
        else if (j==5)
          afit[i*nfit+j]=abest[i*nfit+j]-ebind[ka];
        else
          afit[i*nfit+j]=abest[i*nfit+j];
      }
      xa=xb=xc=xd=0.0f;
      for (j=1;j<npeak;++j)
      { if (xa>afit[i*nfit+j*3+5]) xa=afit[i*nfit+j*3+5];
        if (xb<afit[i*nfit+j*3+5]) xb=afit[i*nfit+j*3+5];
      }
      for (j=0;j<npeak;++j)
      { xe=afit[i*nfit+j*3+4];
        if (xe>1.0e-5) xe=(float)sqrt(2.0/xe);
        if (xc<xe) xc=xe;
        xd+=xe;
      }
      if ((xa>0.0)&&(xb>0.0)) xa=xb;
      else if ((xa<0.0)&&(xb<0.0)) xa=-xa;
      else xa=xb-xa;
      xa+=xd; xc*=2.0f;
      if (xa<xc) xa=xc;
      xa*=2.5f;
      if (xa<8.0) xa=8.0f;
      xb=ebind[ka+ndata-1]-ebind[ka];

      //third pass of fitting
      for (j=0;j<nfit;++j)
      { if (j<1) safit[i*nfit+j]=0.0f;
        else if ((j<2)&&(backassociate==0)) safit[i*nfit+j]=0.01f;
        else if (j<2) safit[i*nfit+j]=0.0f;
        else if ((j==2)&&(xa<xb)) safit[i*nfit+j]=0.0f;
        else if (j==2) safit[i*nfit+j]=0.0f;
        else if ((j%3)==0) safit[i*nfit+j]=0.01f;
        else if (((j%3)==1)&&(peakassociate==0)) safit[i*nfit+j]=0.01f;
        else if ((j%3)==1) safit[i*nfit+j]=0.0f;
        else if ((j>6)&&(peakassociate==0)) safit[i*nfit+j]=0.01f;
        else if (j>6) safit[i*nfit+j]=0.0f;
        else if (afit[i*nfit+j]<0.1) safit[i*nfit+j]=1.0f;
        else safit[i*nfit+j]=0.01f/afit[i*nfit+j];
      }
    }
  }

  if (error==0)
  { for (i=0;i<ncurve;++i)
    { ephoton=curvepar[i*8];
      ka=(int)curvepar[i*8+3];
      kb=(int)curvepar[i*8+4];
      ndata=kb-ka;
      if (ndata>256) ndata=256;

      if ((error==0)&&(nfit>1))
      { Curvefit xpscurve;
        xpscurve.init(ndata,nfit);
        error=xpscurve.loadcurve(fitmath,trymax,tolerance,ebind+ka,
          xpscount+ka,afit+i*nfit,safit+i*nfit,fitxpspeak);
        if (error==0)
          error=xpscurve.dofit(ybest,abest+i*nfit,&reliable,&trynum);
      }
      if (error==0)
      { curvepar[i*8+5]=reliable;
        for (j=0;j<nfit;++j)
        { afit[i*nfit+j]=abest[i*nfit+j];
          if (((j==1)||(j==2))&&(abest[i*nfit]>1.0e-20))
            abest[i*nfit+j]/=abest[i*nfit];
          else if (j==5) abest[i*nfit+j]+=ebind[ka];
        }
      }
      kd=42+(i+1)*22/ncurve-jobcount;
      if ((error==0)&&(kd>0))
      { conout.charfill('-',kd); conout.flush(); jobcount+=kd;
      }
    }
  }

  if ((error==0)&&(ncurve>0)&&(nfit>0))
  { if (bestfit) delete [] bestfit;
    bestfit=new float [ncurve*nfit];
    if (!bestfit) error=102;
    else
    { for (i=0;i<ncurve;++i)
      { for (j=0;j<nfit;++j)
        { bestfit[i*nfit+j]=afit[i*nfit+j];
          if (((j==1)||(j==2))&&(afit[i*nfit]>1.0e-20))
            bestfit[i*nfit+j]/=afit[i*nfit];
          if ((j>6)&&((j%3)==2)) bestfit[i*nfit+j]+=bestfit[i*nfit+5];
        }
      }
    }
  }
  if ((error==0)&&(npeak>0)&&(emission)&&(bestfit))
  { for (i=0;i<ncurve;++i)
    { ephoton=curvepar[i*8];
      ka=(int)curvepar[i*8+3];
      kb=(int)curvepar[i*8+4];
      ndata=kb-ka;
      if (ndata>256) ndata=256;
      for (j=0;j<3;++j) emission[i*numint+j]=curvepar[i*8+j];
      emission[i*numint+3]=curvepar[i*8+5];
      for (j=0;j<npeak;++j)
      { xa=bestfit[i*nfit+j*3+5];
        xc=bestfit[i*nfit+j*3+3];
        xd=bestfit[i*nfit+j*3+4];
        xe=bestfit[i*nfit];
        xa+=ebind[ka];
        if (xd>1.0e-10) xd=(float)(sqrt(2.0/(double)xd));
        else xd=0.0f;
        if (xe>1.0e-10) xc/=xe;
        else xc=0.0f;
        xb=1.25331f*xc*xd;

        emission[i*numint+j*4+4]=xa;
        emission[i*numint+j*4+5]=xb;
        emission[i*numint+j*4+6]=xc;
        emission[i*numint+j*4+7]=xd;
      }
    }
    for (j=0;j<npeak;++j)
    { xa=0.0f;
      for (i=0;i<ncurve;++i)
      { xd=bestfit[i*nfit+j*3+4];
        if (xd>1.0e-10) xd=(float)(sqrt(2.0/(double)xd));
        else xd=0.0f;
        xa+=xd;
      }
      if (ncurve>0) xa/=(float)ncurve;
      xa*=1.25331f;
      for (i=0;i<ncurve;++i) emission[i*numint+j*4+6]*=xa;
    }
    conout.string("-",0,1);
  }

  if (afit) delete [] afit; if (safit) delete [] safit;
  if (abest) delete []abest; if (ybest) delete [] ybest;

  return(error);
} //end of Spectrum::peakfit

int Spectrum::savepeak(char *filename,char *usermessage)
{ int i,j,begrow,linenum,scanuni;
  float xa,xb,ephoton;

  datatype=scanuni=0;
  if ((error==0)&&(emission)&&(ncurve>1))
  { if (emission[0]!=emission[(ncurve-1)*numint]) datatype=851;
    else if (emission[1]!=emission[(ncurve-1)*numint+1]) datatype=852;
    else if (emission[2]!=emission[(ncurve-1)*numint+2]) datatype=853;
    else datatype=851;
  }
  else if (error==0) datatype=851;
  if (error==0)
  { scanuni=datatype%10;
    if (scanuni>3) scanuni=1;
  }

  if ((error==0)&&(emission))
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { i=0; begrow=0; linenum=ncurve;
      while (usermessage[i]!='\0')
      { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
        ++i;
      }
      begrow+=9;
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(linenum,6); fileout.space(5);
      if (linenum>0)
        fileout.string("datakind beginning-row linenumbers",0,1);
      else fileout.string("datakind beginning-row multi-curves",0,1);
      fileout.string(usermessage,0,1);
      fileout.string("   xps spectrum peak intensities");
      if (scanuni==1)
      { fileout.string("   photoemission energy scan curves",0,16+1);
        fileout.string(
          "     parameters: ncurve npoint npeak theta phi",
          0,1);
      }
      else if (scanuni==2)
      { fileout.string("   photoemission polar scan curves",0,16+1);
        fileout.string(
          "     parameters: ncurve npoint npeak wavenum phi",
          0,1);
      }
      else if (scanuni==3)
      { fileout.string("   photoemission azimuthal scan curves",
          0,16+1);
        fileout.string(
          "     parameters: ncurve npoint npeak wavenum theta",
          0,1);
      }
      else
      { fileout.string("   photoelectron diffraction curves",0,16+1);
        fileout.string(
          "     parameters: ncurve npoint npeak wavenum",
          0,1);
      }
      fileout.string("     columns: ");
      for (j=0;j<npeak;++j)
      { if (scanuni==1) fileout.string(" wavenum");
        else if (scanuni==2) fileout.string(" theta");
        else if (scanuni==3) fileout.string(" phi");
        if ((npeak<3)||(j==0)) fileout.string(" intensity");
        else fileout.string(" int");
      }
      fileout.string(" reliability",0,1);
      fileout.integer(1,8);
      fileout.integer(ncurve,8);
      fileout.string("       ncurve npoint",0,1);
      if (scanuni!=1)
      { xa=emission[0]-emission[4];
        if (xa>1.0e-5) xa=0.512f*(float)sqrt(xa);
        else xa=0.0f;
      }
      else xa=emission[1];
      if (scanuni!=3) xb=emission[2];
      else xb=emission[1];
      fileout.integer(1,8);
      fileout.integer(ncurve,8);
      fileout.integer(npeak,8);
      if (scanuni!=1) fileout.floating(xa,8.2f,256);
      else fileout.floating(xa,8.1f,256);
      fileout.floating(xb,8.1f,256);
      fileout.string("       ------------------",0,1);
      for (i=0;i<ncurve;++i)
      { for (j=0;j<npeak;++j)
        { if (scanuni==1)
          { xa=emission[i*numint]-emission[i*numint+j*4+4];
            if (xa>1.0e-5) xa=0.512f*(float)sqrt(xa);
            else xa=0.0f;
          }
          else xa=emission[i*numint+scanuni-1];
          if (scanuni==1) fileout.floating(xa,8.2f,256);
          else fileout.floating(xa,8.1f,256);
          fileout.floating(emission[i*numint+j*4+5],10.4f,256);
        }
        fileout.floating(emission[i*numint+3],10.4f,256+1);
      }
      if (error==0) error=fileout.geterror();
    }

    if (error==0)
    { fileout.charfill('-',65,16+1);
      fileout.string("   X-ray photoemission spectrum peaks",0,2);
      fileout.string(
        "     columns: photon-energy(eV) theta phi reliability",0,1);
      fileout.string(
        "              peak-binding(eV) peak-int peak-width");
      if (npeak==1) fileout.newline();
      else
      { fileout.string(" for peaks 1-");
        fileout.integer(npeak,0,1);
      }
      fileout.integer(ncurve,8);
      fileout.integer(npeak,8);
      fileout.integer(npeak*3+4,8);
      fileout.string("         ndata npeak numbers-per-line",0,1);
      for (i=0;i<ncurve;++i)
      { for (j=0;j<numint;++j)
        { if ((j>4)&&((j%4)==2)) continue;
          if ((j>0)&&((j%8)==0)) fileout.space(6,16);
          if (j==0) xa=10.2f;
          else if (j<3) xa=8.1f;
          else if (j==3) xa=10.4f;
          else if ((j%4)==0) xa=10.2f;
          else if ((j%4)==1) xa=10.4f;
          else xa=10.3f;
          fileout.floating(emission[i*numint+j],xa,256);
        }
        fileout.newline();
      }
    }

    if ((error==0)&&(bestfit))
    { fileout.charfill('-',65,16+1);
      fileout.string("best-fit parameters",0,2);
      if (nfit>0) fileout.string("  function: y=a0");
      if (nfit>1) fileout.string("(1+a1 int((ydata-a0)dx)");
      if (nfit>2) fileout.string("+a2(x-x0)2");
      if (nfit>0) fileout.string(")");
      if (nfit>5)
      { fileout.space(14,16);
        fileout.string("+a3exp(-a4(x-x0-a5)2)");
      }
      if (nfit>8) fileout.string("+a6exp(-a7(x-x0-a8)2)");
      if (nfit>11)
      { fileout.space(14,16);
        fileout.string("+a9exp(-a10(x-x0-a11)2)");
      }
      fileout.newline(2);
      fileout.string("  parameters: curve ephoton reliability a0-a");
      fileout.integer(nfit-1,0,1);
      for (i=0;i<ncurve;++i)
      { ephoton=curvepar[i*8];
        fileout.integer(i+1,8);
        fileout.floating(ephoton,10.2f,256);
        fileout.floating(curvepar[i*8+5],10.4f,256);
        for (j=0;j<nfit;++j)
        { if ((j==3)||(j==9)) fileout.space(8,16);
          fileout.floating(bestfit[i*nfit+j],10.4f,256);
        }
        fileout.newline();
      }
    }
    if (error==0) error=fileout.geterror();
  }

  return(error);
} //end of Spectrum::savepeak

int Spectrum::savelist(char *filename,char *usermessage)
{ int i,j,ka,kb,ndata,begrow,linenum;
  float xa,xb,xc,xd,xe,xs,ephoton,atheta,aphi;

  if ((error==0)&&(npeak>0)&&(bestfit))
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { datatype=842;
      i=0; begrow=0;
      if (ncurve<2) linenum=npoint;
      else linenum=0;
      while (usermessage[i]!='\0')
      { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
        ++i;
      }
      if (linenum>0) begrow+=9;
      else begrow+=7;
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(linenum,6); fileout.space(5);
      if (linenum>0)
        fileout.string("datakind beginning-row linenumbers",0,1);
      else fileout.string("datakind beginning-row multi-curves",0,1);
      fileout.string(usermessage,0,1);
      fileout.string("   X-ray photoemission spectrum best-fit",0,1);
      fileout.string("   counts versus photoelectron kinetic energy");
      fileout.string(" (from Fermi level)",0,1);
      fileout.string(
        "     parameters: curve point photon-energy(eV) theta phi",
        0,1);
      fileout.string(
        "     columns: binding-energy counts background");
      if (npeak>1)
      { fileout.string(" peak 1-");
        fileout.integer(npeak);
      }
      fileout.string(" best-fit",0,1);
      fileout.integer(ncurve,6);
      fileout.integer(npoint,6);
      fileout.string("  ncurve npoint",0,1);

      for (i=0;i<ncurve;++i)
      { ephoton=curvepar[i*8];
        atheta=curvepar[i*8+1];
        aphi=curvepar[i*8+2];
        ka=(int)curvepar[i*8+3];
        kb=(int)curvepar[i*8+4];
        ndata=kb-ka;
        if (ndata>256) ndata=256;
        ka=(int)curvepar[i*8+3];
        fileout.integer(i+1,4);
        fileout.integer(ndata,6);
        fileout.floating(ephoton,10.2f,256);
        fileout.floating(atheta,8.1f,256);
        fileout.floating(aphi,8.1f,256);
        fileout.floating(curvepar[i*8+5],8.3f,256);
        fileout.string("   ------------------",0,1);
        xs=0.0f;
        for (j=0;j<ndata;++j)
        { if (j>0)
          { xs+=xpscount[ka+j]*(ebind[ka+j]-ebind[ka+j-1]);
            xd=xs-bestfit[i*nfit]*(ebind[ka+j]-ebind[ka]);
          }
          else xd=0.0f;
          fileout.floating(ebind[ka+j],10.2f,256);
          fileout.floating(xpscount[ka+j],11.3f,256);
          xa=0.0f; xc=ebind[ka+j]-ebind[ka];
          if (nfit>0) xa+=bestfit[i*nfit];
          if (nfit>1) xa+=bestfit[i*nfit]*bestfit[i*nfit+1]*xd;
          if (nfit>2) xa+=bestfit[i*nfit]*bestfit[i*nfit+2]*xc*xc;
          fileout.floating(xa,9.3f,256);
          xe=0.0f;
          for (kb=0;kb<npeak;++kb)
          { xc=ebind[ka+j]-ebind[ka]-bestfit[i*nfit+kb*3+5];
            xd=bestfit[i*nfit+kb*3+4]*xc*xc;
            if (xd<-30.0) xd=1.0f;
            else if (xd>60.0) xd=0.0f;
            else xd=(float)exp(-xd);
            xb=bestfit[i*nfit+kb*3+3]*xd;
            xe+=xb;
            if (npeak>1) fileout.floating(xa+xb,9.3f,256);
          }
          fileout.floating(xa+xe,9.3f,256+1);
        }
      }
      if (error==0) error=fileout.geterror();
    }
    else error+=(datatype<<16);
  }

  return(error);
} //end of Spectrum::savelist

float fitxpspeak(int fitmath,int ndata,int nfit,int *yafit,
  float *afit,float *safit,float *xdata,float *ydata,float *ymod,
  float *dyda)
{ int i,j,k,mfit,error;
  float xc,xd,xs,reliable;
  float za[20];

  if (nfit>20) error=642;
  else
  { error=0; mfit=nfit/3*3;
    for (j=0;j<20;++j) za[j]=0.0f;
    for (k=0;k<nfit;++k)
    { for (j=k+1;j<nfit;++j)
        if (-yafit[j]==k+1) afit[j]=afit[k]*safit[j];
    }

    xs=0.0f;
    for (i=0;i<ndata;++i)
    { if (i>0) xs+=ydata[i]*(xdata[i]-xdata[i-1]);
      for (k=3;k<mfit;k+=3)
      { xc=xdata[i]-xdata[0]-afit[k+2];
        if (k>3) xc-=afit[5];
        xd=afit[k+1]*xc*xc;
        if (xd<-30.0) xd=1.0f;
        else if (xd>60.0) xd=0.0f;
        else xd=(float)exp(-xd);
        za[k]=xc; za[k+1]=xd;
      }
      xd=xdata[i]-xdata[0];

      if (fitmath==1)
      { if (nfit>0) dyda[i]=1.0f-afit[1]*xd;
        if (nfit>1) dyda[ndata+i]=xs-afit[0]*xd;
        if (nfit>2) dyda[2*ndata+i]=xd*xd;
        for (k=3;k<mfit;k+=3)
        { dyda[k*ndata+i]=za[k+1];
          dyda[(k+1)*ndata+i]=-afit[k]*za[k]*za[k]*za[k+1];
          dyda[(k+2)*ndata+i]=2.0f*afit[k]*afit[k+1]*za[k]*za[k+1];
        }
        for (k=6;k<mfit;k+=3)
          dyda[(3+2)*ndata+i]+=dyda[(k+2)*ndata+i];
        for (k=mfit;k<nfit;++k) dyda[k*ndata+i]=0.0f;
        for (k=0;k<nfit;++k)
        { if (yafit[k]>=0)
          { for (j=k+1;j<nfit;++j)
            { if (-yafit[j]==k+1)
                dyda[k*ndata+i]+=dyda[j*ndata+i]*safit[j];
            }
          }
        }
      }
      xc=0.0f;
      if (nfit>0) xc+=afit[0]*(1.0f-afit[1]*xd);
      if (nfit>1) xc+=afit[1]*xs;
      if (nfit>2) xc+=afit[2]*xd*xd;
      for (k=3;k<mfit;k+=3) xc+=afit[k]*za[k+1];
      ymod[i]=xc;
    }
  }
  xc=xd=0.0f;
  if (error==0)
  { for (i=0;i<ndata;++i)
    { xc+=(ymod[i]-ydata[i])*(ymod[i]-ydata[i]);
      xd+=ymod[i]*ymod[i]+ydata[i]*ydata[i];
    }
  }
  if ((error==0)&&(xd>1.0e-10))
  { reliable=xc/xd;
    if (reliable>10.0) reliable=10.0f;
  }
  else if (error==0) reliable=0.0f;
  else reliable=(float)error;

  return(reliable);
} //end of fitxpspeak

