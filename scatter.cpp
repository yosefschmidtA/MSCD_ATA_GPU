#include <iostream>
#include <fstream>

#include "userutil.h"
#include "phase.h"
#include "rotamat.h"
#include "msfuncs.h"
#include "scatter.h"

Scatter::Scatter()
{ error=107; init();
} //end of Scatter::Scatter

Scatter::~Scatter()
{ if (atomname) delete [] atomname;
  if (symbolname) delete [] symbolname;
  if (xbeta) delete [] xbeta; if (uevenelem) delete [] uevenelem;
  if (tevenelem) delete [] tevenelem;
  if (phaseshift) delete phaseshift; if (evenmat) delete evenmat;
  if (hanka) delete hanka; if (hankb) delete hankb;
} //end of Scatter::~Scatter

void Scatter::init(int alnum,int araorder)
{ int k;
  if ((error==107)||(error==103))
  { datatype=931; betanum=lnum=raorder=0; wavevec=alength=blength=0.0f;
    atomname=symbolname=NULL; xbeta=uevenelem=NULL; tevenelem=NULL;
    phaseshift=NULL; evenmat=NULL; hanka=hankb=NULL;
  }
  if (error==103)
  { if (alnum<0) lnum=0;
    else if (alnum>60) lnum=60;
    else lnum=alnum;
    if (araorder<0) raorder=0;
    else if (araorder>4) raorder=4;
    else raorder=araorder;
    if (raorder==4) eledim=15;
    else if (raorder==3) eledim=10;
    else if (raorder==2) eledim=6;
    else if (raorder==1) eledim=3;
    else eledim=1;
    if (error==0) k=sizeof(int); else k=sizeof(int);
    if ((k<4)&&(lnum>20)) lnum=20;
    if (k<4) betanum=61; else betanum=61;
    atomname=new char [80]; symbolname=new char [80];
    xbeta=new float [betanum];
    uevenelem=new float [betanum*eledim*eledim*2];
    tevenelem=new Fcomplex [betanum*eledim*eledim];
    if ((!atomname)||(!symbolname)||(!xbeta)||
      (!uevenelem)||(!tevenelem)) error=102;
  }
  else if (error==106)
  { if (atomname)
    { atomname=new char [80];
      if (!atomname) error=102;
    }
    if (symbolname)
    { symbolname=new char [80];
      if (!symbolname) error=102;
    }
    if (xbeta)
    { xbeta=new float [betanum];
      if (!xbeta) error=102;
    }
    if (uevenelem)
    { uevenelem=new float [betanum*eledim*eledim*2];
      if (!uevenelem) error=102;
    }
    if (tevenelem)
    { tevenelem=new Fcomplex [betanum*eledim*eledim];
      if (!tevenelem) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Scatter::init

int Scatter::getlength()
{ int ka;
  if (error==0)
  { ka=sizeof(int)+sizeof(Scatter);
    if (atomname) ka+=80;
    if (symbolname) ka+=80;
    if (xbeta) ka+=betanum*sizeof(float);
    if (uevenelem) ka+=betanum*eledim*eledim*2*sizeof(float);
    if (tevenelem) ka+=betanum*eledim*eledim*sizeof(Fcomplex);
  }
  else ka=0;
  return(ka);
} //end of Scatter::getlength

int Scatter::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Scatter);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=80;
    if (ka>=0) ka=memorysend(dest,atomname,ka,kb,length);
    kb=80;
    if (ka>=0) ka=memorysend(dest,symbolname,ka,kb,length);
    kb=betanum*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)xbeta,ka,kb,length);
    kb=betanum*eledim*eledim*2*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)uevenelem,ka,kb,length);
    kb=betanum*eledim*eledim*sizeof(Fcomplex);
    if (ka>=0) ka=memorysend(dest,(char *)tevenelem,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Scatter::paexport

int Scatter::paimport(char *source,int length)
{ int ka,kb;

  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Scatter);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=80;
    if (ka>=0) ka=memoryrec(atomname,source,ka,kb,length);
    kb=80;
    if (ka>=0) ka=memoryrec(symbolname,source,ka,kb,length);
    kb=betanum*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)xbeta,source,ka,kb,length);
    kb=betanum*eledim*eledim*2*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)uevenelem,source,ka,kb,length);
    kb=betanum*eledim*eledim*sizeof(Fcomplex);
    if (ka>=0) ka=memoryrec((char *)tevenelem,source,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Scatter::paimport

int Scatter::getnumpoint()
{ return(betanum*eledim*eledim);
} //end of Scatter::getnumpoint

float Scatter::getmemory()
{ float xa;
  xa=(float)getlength();
  return(xa);
} //end of Scatter::getmemory

int Scatter::loadcurve(char *filename,char *iatomname,
  char *isymbolname)
{
  if ((error==101)||(error==0))
  { error=0;
    stringcopy(atomname,iatomname,25);
    stringcopy(symbolname,isymbolname,10);
  }

  if (error==0)
  { phaseshift=new Phaseshift;
    phaseshift->init();
    if (error==0)
      error=phaseshift->loadcurve(filename,atomname,symbolname,lnum);
    if (error==0) lnum=phaseshift->getlnum();
  }

  return(error);
} //end of Scatter::loadcurve

Fcomplex Scatter::sevenelem(int ma,int na,int mb,int nb,float akin,
  float vka,float vkb,float beta)
{ int al,ka,kb,kelem,kharm,alnum;
  float xa;
  Fcomplex cxa,cxb,cxc,cxd;

  if (error==0)
  { ka=iabs(na); kb=iabs(mb)+iabs(nb);
    kelem=evenmat->getkelem(ma,mb);
    kharm=evenmat->getkharm(ma,mb);
    alnum=phaseshift->getalnum(akin);
    cxa=0.0f;
    for (al=0;al<alnum;++al)
    { xa=evenmat->rotharma(al,kelem,kharm,beta);
      cxb=phaseshift->fsinexpa(akin,al);
      cxc=hankb->fhankelfaca(al,kb,vkb);
      cxd=hanka->fhankelfaca(al,ka,vka);
      cxa+=xa*cxb*cxc*cxd;
    }
  }

  return(cxa);
} //end of Scatter::sevenelem

int Scatter::makefactor(float akin,float lenga,float lengb)
{ int i,j,k,p,q,t,ma,mb,na,nb,elenum;
  float vka,vkb,beta;
  Fcomplex cxa;
  int lamda[64];
  const float radian=(float)(3.1415926535/180.0);

  if (akin<3.0) akin=3.0f;
  else if (akin>25.0) akin=25.0f;
  if (lenga<1.0) lenga=1.0f;
  if (lengb<1.0) lengb=1.0f;
  vka=1.0f/(akin*lenga); vkb=1.0f/(akin*lengb);
  wavevec=akin; alength=lenga; blength=lengb;
  elenum=eledim*eledim;

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
  
  if (error==0)
  { evenmat=new Rotamat;
    evenmat->init(lnum,raorder);
    if (error==0) error=evenmat->makecurve();
  }

  if (error==0)
  { hanka=new Hankel; hankb=new Hankel;
    if (error==0) hanka->init(101,lnum,raorder+1);
    if (error==0) hankb->init(101,lnum,raorder+1);
    if (error==0) error=hanka->makecurve();
    if (error==0) error=hankb->makecurve();
  }

  for (i=0;(error==0)&&(i<betanum);++i)
  { if (betanum>1) beta=float(i*180.0/(betanum-1));
    else beta=0.0f;
    xbeta[i]=beta;
    for (p=0;p<eledim;++p)
    { for (q=0;q<eledim;++q)
      { t=i*elenum+p*eledim+q;
        ma=lamda[p]; na=lamda[32+p];
        mb=lamda[q]; nb=lamda[32+q];
        k=(ma+mb)&1;
        if ((ma==0)&&(mb<0)&&(k==0))
          tevenelem[t]=tevenelem[t-1];
        else if ((ma==0)&&(mb<0))
          tevenelem[t]=-tevenelem[t-1];
        else if ((ma<0)&&(mb==0)&&(k==0))
          tevenelem[t]=tevenelem[t-eledim];
        else if ((ma<0)&&(mb==0))
          tevenelem[t]=-tevenelem[t-eledim];
        else if ((ma<0)&&(mb>0)&&(k==0))
          tevenelem[t]=tevenelem[t-eledim+1];
        else if ((ma<0)&&(mb>0))
          tevenelem[t]=-tevenelem[t-eledim+1];
        else if ((ma<0)&&(mb<0)&&(k==0))
          tevenelem[t]=tevenelem[t-eledim-1];
        else if ((ma<0)&&(mb<0))
          tevenelem[t]=-tevenelem[t-eledim-1];
        else
        { cxa=sevenelem(ma,na,mb,nb,akin,vka,vkb,beta);
          tevenelem[t]=vka*cxa;
        }
        cxa=tevenelem[t];
        uevenelem[t+t]=cabs(cxa);
        uevenelem[t+t+1]=arg(cxa)/radian;
      }
    }
  }

  return(error);
} //end of Scatter::makefactor

int Scatter::savecurve(char *filename,char *usermessage)
{ int i,p,q,t,ma,na,begrow,elenum;
  float xa,xb;
  int lamda[64];

  elenum=eledim*eledim;

  for (p=0;(error==0)&&(p<32);++p)
  { if ((p==0)||(p==3)||(p==10)) ma=0;
    else if (p<3) ma=3-p*2;
    else if (p<6) ma=18-p*4;
    else if (p<8) ma=13-p*2;
    else if (p<10) ma=51-p*6;
    else if (p<13) ma=46-p*4;
    else if (p<15) ma=108-p*8;
    else ma=0;

    if (p==10) na=2;
    else if ((p==3)||(p==6)||(p==7)||(p==11)||(p==12)) na=1;
    else na=0;

    lamda[p]=ma; lamda[32+p]=na;
  }

  if (error==0)
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { i=0;
      if (eledim<2) begrow=9; else begrow=12;
      while (usermessage[i]!='\0')
      { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
        ++i;
      }
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(0,6); fileout.space(5);
      fileout.string("datakind beginning-row linenumbers",0,1);
      fileout.string(usermessage,0,1);
      fileout.space(3); fileout.string(symbolname);
      fileout.space(3); fileout.string(atomname);
      fileout.space(1); fileout.string(" scattering factors");
      if (eledim>1)
      { fileout.string(" in ("); fileout.integer(eledim);
        fileout.string("x"); fileout.integer(eledim);
        fileout.string(") Rehr-Albers matrix",0,1);
        fileout.string("no",5);
        for (i=0;i<eledim;++i) fileout.integer(i,4);
        fileout.string("ma",5,16);
        for (i=0;i<eledim;++i) fileout.integer(lamda[i],4);
        fileout.string("na",5,16);
        for (i=0;i<eledim;++i) fileout.integer(lamda[32+i],4);
      }
      fileout.string("   scattering factors versus scattering angles",
        0,32+1);
      fileout.string("     parameters: number of angles and");
      fileout.string(" quantum momenta, wavevector",0,1);
      fileout.string("       (1/angstrom), first and second bonding");
      fileout.string(" length (angstroms)",0,1);
      fileout.string("     columns:  angle (degree) factors");
      if (eledim<2) fileout.string(" phase (degree) (");
      else fileout.string(" (");
      fileout.integer(eledim); fileout.string("x");
      fileout.integer(eledim); fileout.string(")",0,1);
      fileout.integer(betanum,9); fileout.integer(lnum,9);
      fileout.floating(wavevec,9.2f,256);
      fileout.floating(alength,9.2f,256);
      fileout.floating(blength,9.2f,256+1);
      for (i=0;i<betanum;++i)
      { xa=xbeta[i];
        fileout.floating(xa,6.1f,256);
        for (p=0;p<eledim;++p)
        { for (q=0;q<eledim;++q)
          { t=i*elenum+q*eledim+p;
            xa=uevenelem[t+t]; xb=uevenelem[t+t+1];
            if (q==6) fileout.space(26,16);
            else if (q==10) fileout.space(16,16);
            fileout.floating(xa,10.4f,256);
            if (eledim<2) fileout.floating(xb,10.1f,256);
          }
          if (p<eledim-1) fileout.space(6,16);
          else fileout.newline();
        }
      }
      fileout.newline();
      if (error==0) error=fileout.geterror();
    }
    else error+=(datatype<<16);
  }

  return(error);
} //end of Scatter::savecurve

