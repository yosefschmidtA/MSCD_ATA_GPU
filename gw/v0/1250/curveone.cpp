#include <iostream>
#include <math.h>

#include "polation.h"
#include "curveone.h"

Curveone::Curveone()
{ error=107; ndata=nfft=nsmooth=nderiv=0; init();
} // end of Curveone::Curveone

void Curveone::init(int indata)
{ if ((error==107)||(error==103))
  { xdata=ydata=ysmooth=yfderiv=ysderiv=xfft=yfft=zfft=NULL;
  }
  if (error==103)
  { if (indata<3) error=701;
    else if (indata>1024) error=702;
    else
    { error=101; ndata=indata; nfft=256; //nfft must be a power of 2
      xdata=new float [ndata]; ydata=new float [ndata];
      ysmooth=new float [ndata];
      yfderiv=new float [ndata]; ysderiv=new float [ndata];
      xfft=new float [nfft]; yfft=new float [nfft];
      zfft=new float [nfft];
      if ((!xdata)||(!ydata)||(!ysmooth)||(!yfderiv)||(!ysderiv)||
        (!xfft)||(!yfft)||(!zfft)) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} // end of Curveone::init

Curveone::~Curveone()
{ if (xdata) delete [] xdata; if (ydata) delete [] ydata;
  if (ysmooth) delete [] ysmooth;
  if (yfderiv) delete [] yfderiv; if (ysderiv) delete [] ysderiv;
  if (xfft) delete [] xfft; if (yfft) delete [] yfft;
  if (zfft) delete [] zfft;
} // end of Curveone::~Curveone

float Curveone::getmemory()
{ float xa;
  if (error==0) xa=(float)(sizeof(Curveone)+ndata*5*sizeof(float)+
    nfft*3*sizeof(float));
  else xa=-1.0f;
  return(xa);
} //end of Curveone::getmemory

int Curveone::loadcurve(float *ixdata,float *iydata)
{ int i;

  if ((error==101)||(error==0))
  { error=0;
    for (i=0;i<ndata;++i)
    { xdata[i]=ixdata[i]; ydata[i]=iydata[i];
    }
    for (i=0;i<ndata-1;++i)
    { if (xdata[i+1]-xdata[i]<1.0e-10)
      { i=ndata; error=707;
      }
    }
    nsmooth=nderiv=0;
  }

  return(error);
} //end of Curveone::loadcurve

int Curveone::smoothcurve(int insmooth)
{ int i;
  float xa;

  Polation aspline; aspline.init(ndata);
  if (error==0) error=aspline.loadcurve(xdata,ydata);
  if (error==0)
  { for (i=0;i<nfft;++i)
    { xa=xdata[0]+(float)i*(xdata[ndata-1]-xdata[0])/float(nfft-1);
      yfft[i]=aspline.fspline((float)xa);
    }
    error=fftsmooth(insmooth);
    for (i=0;i<5;++i) yfft[i]=yfft[5];
    for (i=0;i<5;++i) yfft[nfft-1-i]=yfft[nfft-1-5];
  }

  Polation bspline; bspline.init(nfft);
  if (error==0)
  { for (i=0;i<nfft;++i) xfft[i]=float(i);
    error=bspline.loadcurve(xfft,yfft);
  }
  if (error==0)
  { for (i=0;i<ndata;++i)
    { xa=(float)i*(nfft-1)/float(ndata-1);
      ysmooth[i]=bspline.fspline((float)xa);
    }
  }

  return(error);
} //end of Curveone::smoothcurve

int Curveone::fftsmooth(int insmooth)
{ int i;

  float xa,xb;

  if (error==0)
  { xa=(float)(insmooth*nfft)/(float)(ndata);
    nsmooth=(int)xa;
    if (nsmooth<3) nsmooth=3;
    else if (nsmooth>nfft/2) nsmooth=nfft/2;
    error=realfft(0);
  }
  xb=0.0f;
  if (error==0)
  { xa=4.0f/float(nsmooth*nsmooth);
    for (i=0;i<nfft;++i)
    { zfft[i]=yfft[i];
      yfft[i]=(float)exp((float)(-xa*(i-nfft/2+0.5f)*(i-nfft/2+0.5f)));
      xb+=yfft[i];
    }
    if (xb<1.0e-20) xb=1.0f;
    for (i=0;i<nfft;++i) yfft[i]/=xb;
    if (error==0) error=realfft(0);
  }
  if (error==0)
  { for (i=0;i<nfft/2;++i)
    { xa=yfft[i*2]; xb=yfft[i*2+1];
      xa=(float)sqrt((double)(xa*xa+xb*xb));
      yfft[i*2]=zfft[i*2]*xa;
      yfft[i*2+1]=zfft[i*2+1]*xa;
    }
    error=realfft(1);
  }

  return(error);
} //end of Curveone::fftsmooth

int Curveone::derivative(int insmooth)
{ int i,ka;
  float xa;
  float *tempfft;

  tempfft=NULL;
  if (error==0)
  { tempfft=new float [nfft];
    if (!tempfft) error=102;
  }
  if ((error==0)&&(nsmooth==0)) error=smoothcurve(insmooth);
  if (error==0)
  { nderiv=1;
    xa=2.0f*float(nfft)/float(ndata);
    ka=(int)(xa+0.5);
    if (ka<2) ka=2;
    else if (ka>10) ka=10;
  }
  if (error==0)
  { for (i=0;i<nfft;++i) tempfft[i]=yfft[i];
    xa=2.0f*(xdata[ndata-1]-xdata[0])/float(nfft);
    if (xa<1.0e-20) xa=1.0f;
    else xa=1.0f/xa;
    for (i=ka;i<nfft-ka;++i)
    { yfft[i]=xa*(tempfft[i+1]-tempfft[i-1]);
    }
    for (i=0;i<ka;++i) yfft[i]=yfft[ka];
    for (i=nfft-ka;i<nfft;++i) yfft[i]=yfft[nfft-ka-1];
    if (insmooth>0) error=fftsmooth(insmooth);
  }
  if (error==0)
  { Polation aspline; aspline.init(nfft);
    if (error==0)
    { for (i=0;i<nfft;++i) xfft[i]=float(i);
      error=aspline.loadcurve(xfft,yfft);
    }
    if (error==0)
    { for (i=0;i<ndata;++i)
      { xa=(float)i*(nfft-1)/float(ndata-1);
        yfderiv[i]=aspline.fspline((float)xa);
      }
    }
  }

  if (error==0)
  { xa=(xdata[ndata-1]-xdata[0])/float(nfft);
    xa*=xa;
    if (xa<1.0e-20) xa=1.0f;
    else xa=1.0f/xa;
    for (i=ka;i<nfft-ka;++i)
    { yfft[i]=xa*(tempfft[i+1]+tempfft[i-1]-tempfft[i]*2.0f);
    }
    for (i=0;i<ka;++i) yfft[i]=yfft[ka];
    for (i=nfft-ka;i<nfft;++i) yfft[i]=yfft[nfft-ka-1];
    if (insmooth>0) error=fftsmooth(insmooth);
  }
  if (error==0)
  { Polation aspline; aspline.init(nfft);
    if (error==0)
    { for (i=0;i<nfft;++i) xfft[i]=float(i);
      error=aspline.loadcurve(xfft,yfft);
    }
    if (error==0)
    { for (i=0;i<ndata;++i)
      { xa=(float)i*(nfft-1)/float(ndata-1);
        ysderiv[i]=aspline.fspline((float)xa);
      }
    }
  }

  if (tempfft) delete [] tempfft;
  return(error);
} //end of Curveone::derivative

float Curveone::fsmooth(float xa)
{ int i,ka,kb;
  float ya;

  if ((error==0)&&(nsmooth==0)) error=smoothcurve(1);
  if (error==0)
  { ka=0; kb=ndata-1;
    while (kb-ka>1)
    { i=(kb+ka)/2;
      if (xdata[i] > xa) kb=i;
      else ka=i;
    }
    if (ka > ndata-2) ka=ndata-2;
    if (kb < 1) kb=1;
    if (kb <= ka) kb=ka+1;
    ya=ysmooth[ka]+(xa-xdata[ka])*(ysmooth[kb]-ysmooth[ka])/
      (xdata[kb]-xdata[ka]);
  }
  else ya=0.0f;
  return(ya);
} //end of Curveone::fsmooth

float Curveone::ffderiv(float xa)
{ int i,ka,kb;
  float ya;

  if ((error==0)&&(nsmooth==0)) error=smoothcurve(5);
  if ((error==0)&&(nderiv==0)) error=derivative(5);
  if (error==0)
  { ka=0; kb=ndata-1;
    while (kb-ka>1)
    { i=(kb+ka)/2;
      if (xdata[i] > xa) kb=i;
      else ka=i;
    }
    if (ka > ndata-2) ka=ndata-2;
    if (kb < 1) kb=1;
    if (kb <= ka) kb=ka+1;
    ya=yfderiv[ka]+(xa-xdata[ka])*(yfderiv[kb]-yfderiv[ka])/
      (xdata[kb]-xdata[ka]);
  }
  else ya=0.0f;
  return(ya);
} //end of Curveone::ffderiv

float Curveone::fsderiv(float xa)
{ int i,ka,kb;
  float ya;

  if ((error==0)&&(nsmooth==0)) error=smoothcurve(5);
  if ((error==0)&&(nderiv==0)) error=derivative(5);
  if (error==0)
  { ka=0; kb=ndata-1;
    while (kb-ka>1)
    { i=(kb+ka)/2;
      if (xdata[i] > xa) kb=i;
      else ka=i;
    }
    if (ka > ndata-2) ka=ndata-2;
    if (kb < 1) kb=1;
    if (kb <= ka) kb=ka+1;
    ya=ysderiv[ka]+(xa-xdata[ka])*(ysderiv[kb]-ysderiv[ka])/
      (xdata[kb]-xdata[ka]);
  }
  else ya=0.0f;
  return(ya);
} //end of Curveone::fsderiv

int Curveone::getsdmini(float xa,float *xm,float *ym)
{ int i,j,ka,kb,kc;
  float xb,xc,xd;

  if ((error==0)&&(nsmooth==0)) error=smoothcurve(5);
  if ((error==0)&&(nderiv==0)) error=derivative(5);
  if (error==0)
  { ka=0; kb=ndata-1;
    while (kb-ka>1)
    { i=(kb+ka)/2;
      if (xdata[i] > xa) kb=i;
      else ka=i;
    }
    if (ka > ndata-3) ka=ndata-3;
    if (ka < 2) ka=2;

    kb=kc=ka;
    xb=xc=ysderiv[ka-1]+ysderiv[ka]+ysderiv[ka+1];
    for (i=ka-1;i>=2;--i)
    { xd=ysderiv[i-1]+ysderiv[i]+ysderiv[i+1];
      if (xb>=xd)
      { xb=xd; kb=i;
      }
      else i=0;
    }
    for (i=ka+1;i<ndata-2;++i)
    { xd=ysderiv[i-1]+ysderiv[i]+ysderiv[i+1];
      if (xc>=xd)
      { xc=xd; kc=i;
      }
      else i=ndata;
    }
    if (kb==ka) j=kc;
    else if (kc==ka) j=kb;
    else if ((kc-ka)>(ka-kb)) j=kb;
    else if ((kc-ka)<(ka-kb)) j=kc;
    else if (xb>xc) j=kc;
    else j=kb;
    *xm=xdata[j];
    *ym=-(ysderiv[j-1]+ysderiv[j]+ysderiv[j+1])/3.0f;
  }
  return(error);
} //end of Curveone::getsdmini

int Curveone::realfft(int backward)
{ int i,ka,kb,kc,kd;
  float xa,xb,xc,xd,xe,xf,theta,wr,wi,wpr,wpi,wtemp;

  if (error==0)
  { xe=0.5f;
    theta=2.0f*3.14159265f/(float)nfft;
    if (backward==0)
    { xf=-0.5f;
      error=dofft(backward);
    }
    else
    { xf=0.5f; theta=-theta;
    }
    wtemp=(float)sin(0.5*(double)theta);
    wpr=-2.0f*wtemp*wtemp;
    wpi=(float)sin((double)theta);
    wr=1.0f+wpr; wi=wpi;
    for (i=1;i<nfft/4;++i)
    { ka=i+i; kb=ka+1; kc=nfft+1-kb; kd=kc+1;
      xa=xe*(yfft[ka]+yfft[kc]);
      xb=xe*(yfft[kb]-yfft[kd]);
      xc=-xf*(yfft[kb]+yfft[kd]);
      xd=xf*(yfft[ka]-yfft[kc]);
      yfft[ka]=xa+wr*xc-wi*xd;
      yfft[kb]=xb+wr*xd+wi*xc;
      yfft[kc]=xa-wr*xc+wi*xd;
      yfft[kd]=-xb+wr*xd+wi*xc;
      wtemp=wr;
      wr=wr*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    if (backward==0)
    { xa=yfft[0];
      yfft[0]=yfft[0]+yfft[1];
      yfft[1]=xa-yfft[1];
    }
    else
    { xa=yfft[0];
      yfft[0]=xe*(yfft[0]+yfft[1]);
      yfft[1]=xe*(xa-yfft[1]);
      error=dofft(backward);
    }
  }

  return(error);
} //end of Curveone::realfft

int Curveone::dofft(int backward)
{ int i,j,k,ma,ncf;
  float xa,atheta,wpr,wpi,wtemp,tempr,tempi,wr,wi;

  ncf=nfft>>1;
  if (error==0)
  { j=0;
    for (i=0;i<ncf*2-1;i+=2)
    { if (i<j)
      { xa=yfft[i]; yfft[i]=yfft[j]; yfft[j]=xa;
        xa=yfft[i+1]; yfft[i+1]=yfft[j+1]; yfft[j+1]=xa;
      }
      k=ncf;
      while ((k>=2)&&(j>=k))
      { j-=k; k>>=1;
      }
      j+=k;
    }
    ma=2;
    while (ma<ncf*2)
    { atheta=2.0f*3.14159265f/(float)ma;
      if (backward!=0) atheta=-atheta;
      wtemp=(float)sin(0.5*(double)atheta);
      wpr=-2.0f*wtemp*wtemp;
      wpi=(float)sin((double)atheta);
      wr=1.0f; wi=0.0f;
      for (k=0;k<ma;k+=2)
      { for (i=k;i<ncf*2-ma-1;i+=ma+ma)
        { j=i+ma;
          tempr=wr*yfft[j]-wi*yfft[j+1];
          tempi=wr*yfft[j+1]+wi*yfft[j];
          yfft[j]=yfft[i]-tempr;
          yfft[j+1]=yfft[i+1]-tempi;
          yfft[i]+=tempr;
          yfft[i+1]+=tempi;
        }
        wtemp=wr;
        wr=wr*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
      }
      ma<<=1;
    }
    if (backward!=0) for (i=0;i<ncf*2;++i) yfft[i]/=float(ncf);
  }

  return(error);
} // end of Curveone::dofft

