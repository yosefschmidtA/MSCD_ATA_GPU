#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "userutil.h"
#include "meanpath.h"

Meanpath::Meanpath()
{ error=107; init();
} //end of Meanpath::Meanpath

Meanpath::~Meanpath()
{ if (imfpc) delete [] imfpc;
} //end of Meanpath::~Meanpath

void Meanpath::init()
{ if ((error==107)||(error==103))
  { datatype=941; npoint=0; imfpc=NULL;
  }
  if (error==103)
  { imfpc=new float [5];
    if (!imfpc) error=102;
  }
  else if (error==106)
  { if (imfpc)
    { imfpc=new float [5];
      if (!imfpc) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Meanpath::init

int Meanpath::getlength()
{ int ka;
  if (error==0) ka=sizeof(Meanpath)+5*sizeof(float)+sizeof(int);
  else ka=0;
  return(ka);
} //end of Meanpath::getlength

int Meanpath::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Meanpath);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=5*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)imfpc,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Meanpath::paexport

int Meanpath::paimport(char *source,int length)
{ int ka,kb;

  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Meanpath);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=5*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)imfpc,source,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Meanpath::paimport

int Meanpath::getnumpoint()
{ int k;
  if (error==0) k=npoint;
  else k=0;
  return(k);
} //end of Meanpath::getnumpoint

float Meanpath::getmemory()
{ float xa;
  if (error==0) xa=(float)(sizeof(Meanpath)+5*sizeof(float));
  else xa=-1.0f;
  return(xa);
} //end of Meanpath::getmemory

int Meanpath::loadparameter(float ivalence,float ibandgap,
  float idensity,float imweight)
{ float plasmon,localeps;

  localeps=1.0e-5f; kenergy=0.0f;
  if ((error==101)||(error==0))
  { error=0;
    if (ivalence<0.0) valence=0.0f;
    else if (ivalence>100.0) valence=100.0f;
    else if (ivalence>1.0) valence=(float)((int)ivalence);
    else valence=ivalence;
    if (ibandgap<0.0) bandgap=0.0f;
    else if (ibandgap>100.0) bandgap=100.0f;
    else bandgap=ibandgap;
    if (idensity<0.0) density=0.0f;
    else if (idensity>100.0) density=100.0f;
    else density=idensity;
    if (imweight<0.0) mweight=0.0f;
    else if (imweight>1000.0) mweight=1000.0f;
    else mweight=imweight;

    if ((valence>1.0-localeps)&&(density>0.1)&&(mweight>0.1))
    { formula=2;
      plasmon=(float)(28.821*sqrt(valence*density/mweight));
      imfpc[0]=plasmon;
      imfpc[1]=(float)(-0.0216+0.944/sqrt(plasmon*plasmon+
        bandgap*bandgap)+7.39e-4*density);
      imfpc[2]=(float)(0.191/sqrt(density));
      imfpc[3]=(float)(1.97-1.0955e-3*plasmon*plasmon);
      imfpc[4]=(float)(53.4-2.5041e-2*plasmon*plasmon);
    }
    else if ((valence>localeps)&&(valence<1.0-localeps)&&
      (bandgap>localeps))
    { formula=1;
      imfpc[0]=bandgap; imfpc[1]=valence;
    }
    else formula=0;
    kenergy=0.0f; invpath=0.0f;
  }
  return(error);
} //end of Meanpath::loadparameter

float Meanpath::finvpath(float wavevector)
{ float xa,xb,kcoef,localeps;
  localeps=1.0e-3f; kcoef=0.512331f;
  xa=wavevector*wavevector/kcoef/kcoef;
  if ((error==0)&&(wavevector<localeps)) error=717;
  else if ((error==0)&&(formula==2)&&
    ((wavevector<kenergy-localeps)||
    (wavevector>kenergy+localeps)))
  { kenergy=wavevector;
    xb=(float)(imfpc[1]*log(imfpc[2]*xa)-imfpc[3]/xa+
      imfpc[4]/xa/xa);
    invpath=xb*imfpc[0]*imfpc[0]/xa;
  }
  else if ((error==0)&&(formula==1)&&
    ((wavevector<kenergy-localeps)||
    (wavevector>kenergy+localeps)))
  { kenergy=wavevector;
    invpath=(float)(1.0/imfpc[0]/pow(xa,imfpc[1]));
  }
  else if ((error==0)&&(formula==0)) invpath=0.0f;
  else if (error!=0) invpath=0.0f;
  return(invpath);
}//end of Meanpath::finvpath

int Meanpath::savecurve(char *filename,char *usermessage)
{ int i,begrow,linenum;
  float xa,xb;
  filebuf foutbuf;

  linenum=npoint=61;
  if (error==0)
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { i=0; begrow=6;
      if (formula==2) ++begrow;
      while (usermessage[i]!='\0')
      { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
        ++i;
      }
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(linenum,6); fileout.space(5);
      fileout.string("datakind beginning-row linenumbers",0,1);
      fileout.string(usermessage,0,1);

      if (formula==2)
      { fileout.string("number of valence electrons =");
        fileout.integer((int)valence,4);
        fileout.space(5); fileout.string("bandgap energy =");
        fileout.floating(bandgap,7.2f,256);
        fileout.string(" eV",0,1);
      }
      else if (formula==1)
      { fileout.string("electron attenuation length =");
        fileout.floating(bandgap,8.4f,256);
        fileout.string("*energy^");
        fileout.floating(valence,5.3f,256);
        fileout.string(" angstroms",0,1);
      }
      else
      { fileout.string("electron wave attenuation due to ");
        fileout.string("inelastic process not considered",0,1);
      }
      if (formula==2)
      { fileout.string("density of bulk =");
        fileout.floating(density,7.2f,256);
        fileout.string(" g/cm3     molecular weight =");
        fileout.floating(mweight,7.1f,256);
        fileout.string(" amu",0,1);
      }
      fileout.string("k (1/angs)",10,16);
      fileout.string("lamda (angstrom)",17,1);

      for (i=0;i<linenum;++i)
      { xa=(float)(3.0+i*0.2); xb=finvpath(xa);
        if (xb>1.0e-5) xb=1.0f/xb;
        else xb=0.0f;
        fileout.floating(xa,7.2f,256);
        fileout.floating(xb,14.4f,1);
      }
      if (error==0) error=fileout.geterror();
    }
    else error+=(datatype<<16);
  }

  return(error);
} //end of Meanpath::savecurve

