//Some special functions in multiple scattering theory
#include <math.h>

#include "userutil.h"
#include "fcomplex.h"
#include "msfuncs.h"

Hankel::Hankel()
{ error=107; init();
}//end of Hankel::Hankel

Hankel::~Hankel()
{ if (hankmat) delete [] hankmat; if (hankarg) delete [] hankarg;
}//end of Hankel::~Hankel

void Hankel::init(int indata,int ilnum,int icmnum)
{ int k;
  if ((error==107)||(error==103))
  { hankmat=hankarg=NULL;
  }
  if (error==103)
  { if (indata<10) ndata=10;
    else if (indata>256) ndata=256;
    else ndata=indata;
    if (ilnum<0) lnum=0;
    else if (ilnum>60) lnum=60;
    else lnum=ilnum;
    if (icmnum<0) cmnum=0;
    else if (icmnum>5) cmnum=5;
    else cmnum=icmnum;
    if (error==0) k=sizeof(int); else k=sizeof(int);
    if ((k<4)&&(lnum>20)) lnum=20;
    hankmat=new Fcomplex [ndata*cmnum*lnum];
    hankarg=new Fcomplex [cmnum*lnum];
    if ((!hankmat)||(!hankarg)) error=102;
  }
  else if (error==106)
  { if (hankmat)
    { hankmat=new Fcomplex [ndata*cmnum*lnum];
      if (!hankmat) error=102;
    }
    if (hankarg)
    { hankarg=new Fcomplex [cmnum*lnum];
      if (!hankarg) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Hankel::init

int Hankel::getlength()
{ int ka;
  if (error==0) ka=sizeof(Hankel)+(ndata+1)*cmnum*lnum*
    sizeof(Fcomplex)+sizeof(int);
  else ka=0;
  return(ka);
} //end of Hankel::getlength

int Hankel::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Hankel);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=ndata*cmnum*lnum*sizeof(Fcomplex);
    if (ka>=0) ka=memorysend(dest,(char *)hankmat,ka,kb,length);
    kb=cmnum*lnum*sizeof(Fcomplex);
    if (ka>=0) ka=memorysend(dest,(char *)hankarg,ka,kb,length);
    if ((ka<0)||(ka!=length)) error=901;
  }
  return(error);
} //end of Hankel::paexport

int Hankel::paimport(char *source,int length)
{ int ka,kb;

  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Hankel);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=ndata*cmnum*lnum*sizeof(Fcomplex);
    if (ka>=0) ka=memoryrec((char *)hankmat,source,ka,kb,length);
    kb=cmnum*lnum*sizeof(Fcomplex);
    if (ka>=0) ka=memoryrec((char *)hankarg,source,ka,kb,length);
    if ((ka<0)||(ka!=length)) error=901;
  }
  return(error);
} //end of Hankel::paimport

int Hankel::makecurve()
{ int i,l,m,t;
  float invkx;
  if ((error==101)||(error==0))
  { error=0;
    for (i=0;i<ndata;++i)
    { invkx=(float)(i/(ndata-1.0)/4.0);
      for (l=0;l<lnum;++l)
      { for (m=0; m<cmnum;++m)
        { t=i*lnum*cmnum+l*cmnum+m;
          if ((l==0)&&(m==0)) hankmat[t]=1.0f;
          else if (l==0) hankmat[t]=0.0f;
          else if ((l==1)&&(m==0)) hankmat[t]=Fcomplex(1.0f,invkx);
          else if ((l==1)&&(m==1)) hankmat[t]=Fcomplex(0.0f,invkx);
          else if (l==1) hankmat[t]=0.0f;
          else if (m==0) hankmat[t]=hankmat[t-cmnum-cmnum]+
            Fcomplex(0.0f,(float)(l+l-1.0)*invkx)*hankmat[t-cmnum];
          else hankmat[t]=hankmat[t-cmnum-cmnum]+
            Fcomplex(0.0f,(float)(l+l-1.0)*invkx)*(hankmat[t-cmnum]+
            hankmat[t-cmnum-1]);
        }
      }
    }
    argument=0.0f;
    for (l=0;l<lnum;++l) for (m=0;m<cmnum;++m)
      hankarg[l*cmnum+m]=hankmat[l*cmnum+m];
  }
  return(error);
}// end of Hankel::makecurve

float Hankel::getmemory()
{ float xa;
  if (error==0)
    xa=(float)(sizeof(Hankel)+(ndata*cmnum*lnum+cmnum*lnum)*
      sizeof(Fcomplex));
  else xa=-1.0f;
  return(xa);
} //end of Hankel::getmemory

Fcomplex Hankel::fhankelfac(int al,int am,float invkx)
{ int i,l,m,t;
  float ya;
  Fcomplex cxa;
  if ((error==0)&&(fabs(invkx-argument)>1.0e-3))
  { argument=invkx;
    ya=(float)(4.0*(ndata-1.0)*invkx);
    i=(int)ya;
    if (i<0) i=0;
    else if (i>ndata-2) i=ndata-2;
    for (l=0;l<lnum;++l) for (m=0;m<cmnum;++m)
    { t=i*lnum*cmnum+l*cmnum+m;
      hankarg[l*cmnum+m]=hankmat[t]+
        (ya-i)*(hankmat[t+lnum*cmnum]-hankmat[t]);
    }
  }
  if ((error==0)&&(al<lnum)&&(am<cmnum)) cxa=hankarg[al*cmnum+am];
  else if (error==0) cxa=0.0f;
  else cxa=(float)error;

  return(cxa);
}//end of Hankel::fhankelfac

Fcomplex Hankel::fhankelfaca(int al,int am,float invkx)
{ int i,l,m,t;
  float ya;
  Fcomplex cxa;
  if (invkx!=argument)
  { argument=invkx;
    ya=(float)(4.0*(ndata-1.0)*invkx);
    i=(int)ya;
    if (i<0) i=0;
    else if (i>ndata-2) i=ndata-2;
    for (l=0;l<lnum;++l) for (m=0;m<cmnum;++m)
    { t=i*lnum*cmnum+l*cmnum+m;
      hankarg[l*cmnum+m]=hankmat[t]+
        (ya-i)*(hankmat[t+lnum*cmnum]-hankmat[t]);
    }
  }
  cxa=hankarg[al*cmnum+am];

  return(cxa);
}//end of Hankel::fhankelfaca

Expix::Expix()
{ error=107; init();
}//end of Expix::Expix

Expix::~Expix()
{ if (cexpix) delete [] cexpix; if (csinexp) delete [] csinexp;
}//end of Expix::~Expix

void Expix::init(int indata)
{ if ((error==107)||(error==103))
  { cexpix=csinexp=NULL;
  }
  if (error==103)
  { if (indata<128) ndata=128;
    else if (indata>8001) ndata=8001;
    else ndata=indata;
    if ((ndata&1)==0) ++ndata;
    mdata=(ndata>>1);
    cexpix=new Fcomplex [ndata];
    csinexp=new Fcomplex [ndata];
    if ((!cexpix)||(!csinexp)) error=102;
  }
  else if (error==106)
  { if (cexpix)
    { cexpix=new Fcomplex [ndata];
      if (!cexpix) error=102;
    }
    if (csinexp)
    { csinexp=new Fcomplex [ndata];
      if (!csinexp) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Expix::init

int Expix::getlength()
{ int ka;
  if (error==0) ka=sizeof(Expix)+ndata*2*sizeof(Fcomplex)+sizeof(int);
  else ka=0;
  return(ka);
} //end of Expix::getlength

int Expix::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Expix);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=ndata*sizeof(Fcomplex);
    if (ka>=0) ka=memorysend(dest,(char *)cexpix,ka,kb,length);
    kb=ndata*sizeof(Fcomplex);
    if (ka>=0) ka=memorysend(dest,(char *)csinexp,ka,kb,length);
    if ((ka<0)||(ka!=length)) error=901;
  }
  return(error);
} //end of Expix::paexport

int Expix::paimport(char *source,int length)
{ int ka,kb;

  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Expix);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=ndata*sizeof(Fcomplex);
    if (ka>=0) ka=memoryrec((char *)cexpix,source,ka,kb,length);
    kb=ndata*sizeof(Fcomplex);
    if (ka>=0) ka=memoryrec((char *)csinexp,source,ka,kb,length);
    if ((ka<0)||(ka!=length)) error=901;
  }
  return(error);
} //end of Expix::paimport

int Expix::makecurve()
{ int k;
  float xa;
  const float radian=(float)(3.14159265/180.0);
  if ((error==101)||(error==0))
  { error=0;
    for (k=0;k<ndata;++k)
    { xa=(float)((k-mdata)*360.0*radian/(ndata-1.0));
      cexpix[k]=Fcomplex((float)cos(xa),(float)sin(xa));
      csinexp[k]=(float)sin(xa)*cexpix[k];
    }
  }
  return(error);
}//end of Expix::makecurve();

float Expix::getmemory()
{ float xa;
  if (error==0) xa=(float)(sizeof(Expix)+ndata*2.0*sizeof(Fcomplex));
  else xa=-1.0f;
  return(xa);
} //end of Expix::getmemory

//xa use unit of degree
Fcomplex Expix::fexpix(float xa)
{ int k;
  while (xa>180.0) xa-=360.0f;
  while (xa<-180.0) xa+=360.0f;
  k=(int)(mdata+(ndata-1.0)*xa/360.0+0.5);
  if (k<0) k=0;
  else if (k>ndata-1) k=ndata-1;
  return(cexpix[k]);
}//end of Expix::fexpix

//xa use unit of degree
Fcomplex Expix::fsinexp(float xa)
{ int k;
  while (xa>180.0) xa-=360.0f;
  while (xa<-180.0) xa+=360.0f;
  k=(int)(mdata+(ndata-1.0)*xa/360.0+0.5);
  if (k<0) k=0;
  else if (k>ndata-1) k=ndata-1;
  return(csinexp[k]);
}//end of Expix::fsinexp

