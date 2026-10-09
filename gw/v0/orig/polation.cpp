#include <math.h>

#include "polation.h"

Polation::Polation()
{ error=107; init();
} // end of Polation::Polation

void Polation::init(int indata)
{ if ((error==107)||(error==103))
  { xdata=ydata=y2data=udata=NULL;
  }
  if (error==103)
  { if (indata<3) error=701;
    else if (indata>1024) error=702;
    else
    { error=101; ndata=indata; model=0;
      xdata=new float [ndata]; ydata=new float [ndata];
      y2data=new float [ndata]; udata=new float [ndata];
      if ((!xdata)||(!ydata)||(!y2data)||(!udata)) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} // end of Polation::init

Polation::~Polation()
{ if (xdata) delete [] xdata; if (ydata) delete [] ydata;
  if (y2data) delete [] y2data; if (udata) delete [] udata;
} // end of Polation::~Polation

float Polation::getmemory()
{ float xa;
  if (error==0) xa=(float)(sizeof(Polation)+ndata*4*sizeof(float));
  else xa=-1.0f;
  return(xa);
} //end of Polation::getmemory

int Polation::loadcurve(float *ixdata,float *iydata)
{ int i;

  yp1=ypn=1.0e30f;
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
    model=0;
  }

  return(error);
} //end of Polation::loadcurve

float Polation::fspline(float x)
{ int k,klow,khigh;
  float h,a,b,y;

  if ((error==0)&&(model!=1))
  { error=makespline();
    if (error==0) model=1;
  }
  if (error==0)
  { klow=0; khigh=ndata-1;
    while (khigh-klow>1)
    { k=(khigh+klow)/2;
      if (xdata[k] > x) khigh=k;
      else klow=k;
    }
    if (klow > ndata-2) klow=ndata-2;
    if (khigh < 1) khigh=1;
    if (khigh <= klow) khigh=klow+1;
    h=xdata[khigh]-xdata[klow];
    a=(xdata[khigh]-x)/h;
    b=(x-xdata[klow])/h;
    y=a*ydata[klow]+b*ydata[khigh]+
      ((a*a*a-a)*y2data[klow]+(b*b*b-b)*y2data[khigh])*h*h/6.0f;
  }
  else y=0.0f;
  return(y);
} //end of Polation::fspline

int Polation::makespline()
{ int i,k;
  float sig,p,qn,un;

  if (error==0)
  { if (yp1>=1.0e10) y2data[0]=udata[0]=0.0f;
    else
    { y2data[0]=-0.5f;
      udata[0]=(3.0f/(xdata[1]-xdata[0]))*
        ((ydata[1]-ydata[0])/(xdata[1]-xdata[0])-yp1);
    }
    for (i=1;i<ndata-1;++i)
    { sig=(xdata[i]-xdata[i-1])/(xdata[i+1]-xdata[i-1]);
      p=sig*y2data[i-1]+2.0f;
      y2data[i]=(sig-1.0f)/p;
      udata[i]=(ydata[i+1]-ydata[i])/(xdata[i+1]-xdata[i])-
        (ydata[i]-ydata[i-1])/(xdata[i]-xdata[i-1]);
      udata[i]=(6.0f*udata[i]/(xdata[i+1]-xdata[i-1])-
        sig*udata[i-1])/p;
    }
    if (ypn>=1.0e10) qn=un=0.0f;
    else
    { qn=0.5f;
      un=(3.0f/(xdata[ndata-1]-xdata[ndata-2]))*(ypn-(ydata[ndata-1]-
        ydata[ndata-2])/(xdata[ndata-1]-xdata[ndata-2]));
    }
    y2data[ndata-1]=(un-qn*udata[ndata-2])/(qn*y2data[ndata-2]+1.0f);
    for (k=ndata-2;k>=0;--k) y2data[k]=y2data[k]*y2data[k+1]+udata[k];
  }
  return(error);
} // end of Polation::makespline

float fitspline(int fitmath,int ndata,int nfit,int *yafit,
  float *afit,float *safit,float *xdata,float *ydata,float *ymod,
  float *dyda)
{ int i,j,k,error;
  float xc,xd,amin,amax,astep,reliable;
  float xafit[20];

  if (nfit>20) error=642;
  else
  { error=0;
    for (i=0;i<nfit;++i)
      xafit[i]=xdata[0]+i*(xdata[ndata-1]-xdata[0])/(nfit-1.0f);
    if (fitmath==1)
    { amin=amax=afit[0];
      for (i=0;i<nfit;++i)
      { if (amin>afit[i]) amin=afit[i];
        if (amax<afit[i]) amax=afit[i];
      }
      astep=(amax-amin)*0.01f; amin=(float)fabs(amax+amin)*0.001f;
      if (astep<amin) astep=amin;
      if (astep<1.0e-10) astep=1.0e-10f;

      for (k=0;k<nfit;++k)
      { if (yafit[k]>=0)
        { afit[k]+=astep;
          for (j=k+1;j<nfit;++j)
            if (-yafit[j]==k+1) afit[j]=afit[k]*safit[j];
          Polation aspline; aspline.init(nfit);
          if (error==0) error=aspline.loadcurve(xafit,afit);
          for (i=0;(error==0)&&(i<ndata);++i)
            dyda[k*ndata+i]=aspline.fspline(xdata[i]);
          afit[k]-=astep+astep;
          for (j=k+1;j<nfit;++j)
            if (-yafit[j]==k+1) afit[j]=afit[k]*safit[j];
          if (error==0) error=aspline.loadcurve(xafit,afit);
          for (i=0;(error==0)&&(i<ndata);++i)
          { ymod[i]=aspline.fspline(xdata[i]);
            dyda[k*ndata+i]=(dyda[k*ndata+i]-ymod[i])*0.5f/astep;
          }
          afit[k]+=astep;
          for (j=k+1;j<nfit;++j)
            if (-yafit[j]==k+1) afit[j]=afit[k]*safit[j];
        }
      }
    }
  }
  xc=xd=0.0f;
  if (error==0)
  { Polation aspline; aspline.init(nfit);
    if (error==0) error=aspline.loadcurve(xafit,afit);
    for (i=0;(error==0)&&(i<ndata);++i)
    { ymod[i]=aspline.fspline(xdata[i]);
      xc+=(ymod[i]-ydata[i])*(ymod[i]-ydata[i]);
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
} //end of fitspline

