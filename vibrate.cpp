#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "userutil.h"
#include "vibrate.h"

Vibration::Vibration()
{ error=107; init();
} //end of Vibration::Vibration

Vibration::~Vibration()
{ if (thermat) delete [] thermat;
} //end of Vibration::~Vibration

void Vibration::init()
{ if ((error==107)||(error==103))
  { thermat=NULL; npoint=0; datatype=951;
  }
  if (error==103)
  { thermat=new float [4000];
    if (!thermat) error=102;
  }
  else if (error==106)
  { if (thermat)
    { thermat=new float [4000];
      if (!thermat) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Vibration::init

int Vibration::getlength()
{ int ka;
  if (error==0) ka=sizeof(Vibration)+4000*sizeof(float)+sizeof(int);
  else ka=0;
  return(ka);
} //end of Vibration::getlength

int Vibration::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Vibration);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=4000*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)thermat,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Vibration::paexport

int Vibration::paimport(char *source,int length)
{ int ka,kb;
  
  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Vibration);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=4000*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)thermat,source,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Vibration::paimport

int Vibration::getnumpoint()
{ int k;
  if (error==0) k=npoint;
  else k=0;
  return(k);
} //end of Vibration::getnumpoint

float Vibration::getmemory()
{ float xa;
  if (error==0)
    xa=(float)(sizeof(Vibration)+4000*sizeof(float));
  else xa=-1.0f;
  return(xa);
} //end of Vibration::getmemory

int Vibration::loadparameter(float idensity,float imweight,
  float itdebye,float itsample)
{ int i,j;
  float x,y,z,qdebye,invalpha,invbeta,inva,localeps,pi;

  if ((error==101)||(error==0))
  { error=0;
    density=idensity; mweight=imweight;
    tdebye=itdebye; tsample=itsample;
    thernum=4000; thermax=40.0f;
    therstep=thermax/float(thernum);
    pi=3.14159265f; localeps=1.0e-5f;
    if ((density > 0.1) && (mweight > 0.1) && (tdebye > 0.1) &&
      (tsample > 0.1))
    { qdebye=(float)(3.29205*pow(density/mweight,1.0/3.0));
      z=72.7918f/tdebye/mweight;
    }
    else
    { qdebye=0.0f; z=0.0f;
    }
    for (i=0;i<thernum;++i)
    { inva=(float)i*therstep;
      y=invalpha=invbeta=0.0f;
      if (qdebye > localeps)
      { invalpha=tsample/tdebye; invbeta=inva/qdebye;
        y=(float)(1.0+4.0*invalpha*invalpha*pi*pi/6.0);
        for (j=1;j<20;++j)
        { x=(float)(exp(-j/invalpha)*(1.0/j/j+1.0/j/invalpha));
          y-=x*4.0f*invalpha*invalpha;
        }
      }
      if ((qdebye > localeps)&&(inva>=therstep))
      { y-=(float)(2.0*(1.0-cos(1.0/invbeta))*invbeta*invbeta);
        for (j=1;j<20;++j)
        { x=(float)(exp(-j/invalpha)*(j*sin(1.0/invbeta)+
            invalpha/invbeta*cos(1.0/invbeta))-invalpha/invbeta);
          x=x/(j*j+invalpha*invalpha/invbeta/invbeta);
          y+=x*4.0f*invalpha*invbeta;
        }
      }
      if (y > 0.0) thermat[i]=y*z;
      else thermat[i]=0.0f;
    }
  }
  return(error);
} //end of Vibration::loadparameter

float Vibration::fvibmsrd(float bondlength,float aweight)
{ int i;
  float y,value;

  if (error==0)
  { if (bondlength<1.0/therstep/thernum) y=thernum-2.0f;
    else y=1.0f/bondlength/therstep;
    if (y < 0.0) y=0.0f;
    else if (y > thernum-2.0) y=thernum-2.0f;
    i=(int)y;
    if ((mweight > 0.1) && (aweight > 0.1) && (tdebye > 0.1) &&
      (tsample > 0.1))
    { value=thermat[i]+(y-(float)i)*(thermat[i+1]-thermat[i]);
      value=value*(float)sqrt(mweight/aweight);
    }
    else value=0.0f;
  }
  else value=0.0f;
  return(value);
} //end of Vibration::fvibmsrd

int Vibration::savecurve(char *filename,char *usermessage,float aweight)
{ int i,ka,kb,begrow;
  float xa,xb,xmin,xmax,xstep;

  npoint=45;
  if (error==0)
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { i=0; begrow=10;
      while (usermessage[i]!='\0')
      { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
        ++i;
      }
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(npoint,6); fileout.space(5);
      fileout.string("datakind beginning-row linenumbers",0,1);
      fileout.string(usermessage,0,1);

      fileout.string(" correlated thermal vibrational mean square ");
      fileout.string("relative displacement",0,1);
      fileout.string("   versus interatomic distance",0,1);

      fileout.string("   density =");
      fileout.floating(density,7.2f,256);
      fileout.string(" g/cm3",0,1);
      fileout.string("   molecular weight =");
      fileout.floating(mweight,7.1f,256);
      fileout.string(" amu     atomic weight =");
      fileout.floating(aweight,7.1f,256);
      fileout.string(" amu",0,1);
      fileout.string("   debye and sample temperature =");
      fileout.floating(tdebye,7.1f,256);
      fileout.string(" and "); fileout.floating(tsample,7.1f,256);
      fileout.string(" K",0,1);
      fileout.string("idis",7,16); fileout.string("msrdis",14);
      fileout.string(" (units: angstrom angstrom2)",0,1);

      xmin=1.0f; xmax=xmin+40.0f; ka=npoint*4/9;
      xstep=(xmax-xmin)/8.0f/ka; xb=0.0f;
      for (i=0;i<ka;++i)
      { xb=xmin+(float)i*xstep;
        fileout.floating(xb,7.2f,256);
        fileout.floating(fvibmsrd(xb,aweight),14.5f,512+1);
      }
      xa=xb+xstep; kb=ka+npoint*2/9;
      xstep=(xmax-xmin)/8.0f/(kb-ka);
      for (i=ka;i<kb;++i)
      { xb=xa+(float)(i-ka)*xstep;
        fileout.floating(xb,7.2f,256);
        fileout.floating(fvibmsrd(xb,aweight),14.5f,512+1);
      }
      xa=xb+xstep; xstep=(xmax-xa)/(npoint-kb);
      for (i=kb;i<npoint;++i)
      { xb=xa+(float)(i-kb)*xstep;
        fileout.floating(xb,7.2f,256);
        fileout.floating(fvibmsrd(xb,aweight),14.5f,512+1);
      }
      if (error==0) error=fileout.geterror();
    }
  }
  return(error);
} //end of Vibration::savecurve

