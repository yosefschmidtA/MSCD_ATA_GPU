#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "userutil.h"
#include "phase.h"

Phaseshift::Phaseshift()
{ error=107; init();
} //end of Phaseshift::Phaseshift

Phaseshift::~Phaseshift()
{ if (atomname) delete [] atomname;
  if (symbolname) delete [] symbolname;
  if (phasea) delete [] phasea; if (phaseb) delete [] phaseb;
  if (phasec) delete [] phasec;
} //end of Phaseshift::~Phaseshift

void Phaseshift::init()
{ if ((error==107)||(error==103))
  { datatype=711; ndata=0; palnum=0;
    atomname=symbolname=NULL; phasea=phaseb=NULL; phasec=NULL;
  }
  if (error==103)
  { atomname=new char [25]; symbolname=new char [10];
    phasea=new float [256*61]; phaseb=new float [61];
    phasec=new Fcomplex [61];
    if ((!atomname)||(!symbolname)||(!phasea)||(!phaseb)) error=102;
  }
  else if (error==106)
  { if (atomname)
    { atomname=new char [25];
      if (!atomname) error=102;
    }
    if (symbolname)
    { symbolname=new char [10];
      if (!symbolname) error=102;
    }
    if (phasea)
    { phasea=new float [256*61];
      if (!phasea) error=102;
    }
    if (phaseb)
    { phaseb=new float [61];
      if (!phaseb) error=102;
    }
    if (phasec)
    { phasec=new Fcomplex [61];
      if (!phasec) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Phaseshift::init

int Phaseshift::getlength()
{ int ka;
  if (error==0) ka=sizeof(Phaseshift)+(25+10)*sizeof(char)+
    (256*61+61)*sizeof(float)+61*sizeof(Fcomplex)+
    sizeof(int);
  else ka=0;
  return(ka);
} //end of Phaseshift::getlength

int Phaseshift::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Phaseshift);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=25;
    if (ka>=0) ka=memorysend(dest,atomname,ka,kb,length);
    kb=10;
    if (ka>=0) ka=memorysend(dest,symbolname,ka,kb,length);
    kb=256*61*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)phasea,ka,kb,length);
    kb=61*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)phaseb,ka,kb,length);
    kb=61*sizeof(Fcomplex);
    if (ka>=0) ka=memorysend(dest,(char *)phasec,ka,kb,length);
    if ((ka<0)||(ka!=length)) error=901;
  }
  return(error);
} //end of Phaseshift::paexport

int Phaseshift::paimport(char *source,int length)
{ int ka,kb;
  
  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Phaseshift);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=25;
    if (ka>=0) ka=memoryrec(atomname,source,ka,kb,length);
    kb=10;
    if (ka>=0) ka=memoryrec(symbolname,source,ka,kb,length);
    kb=256*61*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)phasea,source,ka,kb,length);
    kb=61*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)phaseb,source,ka,kb,length);
    kb=61*sizeof(Fcomplex);
    if (ka>=0) ka=memoryrec((char *)phasec,source,ka,kb,length);
    if ((ka<0)||(ka!=length)) error=901;
  }
  return(error);
} //end of Phaseshift::paimport

int Phaseshift::getlnum()
{ int k;
  if (error==0) k=lnum;
  else k=0;
  return(k);
} //end of Phaseshift::getlnum

int Phaseshift::getnumpoint()
{ int k;
  if (error==0) k=ndata;
  else k=0;
  return(k);
} //end of Phaseshift::getnumpoint

float Phaseshift::getmemory()
{ float xa;
  if (error==0)
    xa=(float)(sizeof(Phaseshift)+35*sizeof(char)+
      (256*61+61)*sizeof(float)+61*sizeof(Fcomplex));
  else xa=-1.0f;
  return(xa);
} //end of Phaseshift::getmemory

int Phaseshift::getalnum(float wavevec)
{ int alnum;
  alnum=makephase(wavevec);
  return(alnum);
} //end of Phaseshift::getalnum

int Phaseshift::loadcurve(char *filename,char *iatomname,
  char *isymbolname,int alnum)
{ int i,j,k,begrow;
  float xa,xb,xc;
  const float radian=(float)(3.14159265/180.0);

  if ((error==101)||(error==0))
  { error=0; lnum=0;
    for (i=0;i<256*61;++i) phasea[i]=0.0f;
    std::ifstream fin(filename, std::ios::in);
    if (!fin) error=201+(datatype<<16);
    else
    { fin >> xa; k=(int)xa;
      if (fin && (k==datatype))
      { fin >> xb >> xc; skipendline(fin);
        begrow=(int)xb; ndata=(int)xc;
        for (i=1;i<begrow-1;++i) skipendline(fin);
        if (fin)
        { fin >> xb; k=lnum=(int)xb;
          if ((ndata<=0)||(ndata-1==k))
          { ndata=k; fin >> xb; lnum=(int)xb;
          }
          else --ndata;
          skipendline(fin);
        }
        if (ndata<0) ndata=0;
        else if (ndata>256) ndata=256;
        if (lnum<0) lnum=0;
        else if (lnum>60) lnum=60;
        if ((alnum>0)&&(lnum>alnum)) lnum=alnum;
        for(i=0;fin&&(i<ndata);++i)
        { for (j=0;j<lnum+1;++j)
          { if ((i==0)||(j>0)) fin >> xb;
            else
            { xb=0.0f; xc=phasea[(i-1)*61];
              while (xb<xc) fin >> xb;
            }
            if (j==0) phasea[i*61+j]=xb;
            else phasea[i*61+j]=xb/radian;
          }
          skipendline(fin);
        }
        if (i<ndata) ndata=i;
      }
      else if (fin && (xa > 1.0e5)) error=209+(datatype<<16);
      else if (!fin)
      { fin.clear(); skipendline(fin);
        i=0;
        while (fin && (i<256))
        { fin >> phasea[i*61] >> xb;
          k=(int)xb;
          if (k<0) k=0;
          else if (k+1>60) lnum=60;
          else if (k+1 > lnum) lnum=k+1;
          for (j=0;fin&&(j<k+1);++j)
          { if (j<lnum)
            { fin >> xa >> xb >> xa; phasea[i*61+j+1]=xb/radian;
            }
            else fin >> xa >> xa >> xa;
          }
          if (fin) ++i;
        }
        ndata=i;
      }
      else
      { i=0;
        while (fin && (i<256))
        { if (i==0) phasea[i*61]=xa;
          else fin >> phasea[i*61];
          k=40;
          lnum=60;
          for (j=0;fin&&(j<k);++j)
          { if (j<lnum)
            { fin >> xb; phasea[i*61+j+1]=xb/radian;
            }
            else fin >> xa;
          }
          if (fin&&(phasea[i*61]>1.0)) ++i;
          else break;
        }
        ndata=i;
        for (i=0;i<ndata;++i)
          phasea[i*61]=(float)(0.512331*sqrt(phasea[i*61]));
      }
    }
  }
  if (error==0)
  { i=0; while (iatomname[i]==' ') ++i;
    j=0;
    while ((j<16)&&(iatomname[j+i]!='\0'))
    { atomname[j]=iatomname[j+i];
      ++j;
    }
    atomname[j]='\0';
    i=0; while (isymbolname[i]==' ') ++i;
    j=0;
    while ((j<8)&&(isymbolname[j+i]!='\0'))
    { symbolname[j]=isymbolname[j+i];
      ++j;
    }
    symbolname[j]='\0';
  }
  for (i=0;(i<ndata-1)&&(error==0);++i)
    if (phasea[i*61+61]-phasea[i*61]<1.0e-5) error=707+(datatype<<16);
  if ((error==0)&&(ndata<0)) error=204+(datatype<<16);
  else if (error==0)
  { if ((alnum>0)&&(lnum>alnum)) lnum=alnum;
    if (lnum>60) lnum=60;
    for (j=lnum;j>0;--j)
    { for (i=0;i<ndata;++i)
      { xb=phasea[i*61+j];
        while (xb>180.0) xb-=360.0f;
        while (xb<-180.0) xb+=360.0f;
        if ((xb>1.0e-5)||(xb<-1.0e-5)) i=ndata+1;
      }
      if (i>ndata) break;
    }
    lnum=j;

    pwavea=pwaveb=phasea[0]-100.0f;
  }
  else if (error==0) error=204+(datatype<<16);
  return(error);
} //end of Phaseshift::loadcurve

float Phaseshift::fphase(float wavevec,int lquantum)
{ float xa;
  if ((error==0)&&((lquantum<0)||(lquantum>lnum-1))) xa=0.0f;
  else if ((error==0)&&(fabs(wavevec-pwavea)<1.0e-3))
    xa=phaseb[lquantum+1];
  else if (error==0)
  { makephase(wavevec); xa=phaseb[lquantum+1];
  }
  else xa=0.0f;

  return(xa);
}//end of Phaseshift::fphase

Fcomplex Phaseshift::fsinexp(float wavevec,int lquantum)
{ Fcomplex cxa;
  if ((error==0)&&((lquantum<0)||(lquantum>lnum-1))) cxa=0.0f;
  else if ((error==0)&&(fabs(wavevec-pwaveb)<1.0e-3))
    cxa=phasec[lquantum+1];
  else if (error==0)
  { makesinexp(wavevec); cxa=phasec[lquantum+1];
  }
  else cxa=0.0f;
  return(cxa);
} //end of Phaseshift::fsinexp

int Phaseshift::makephase(float wavevec)
{ int i,j,alnum;
  float xa,xb,xc,ya;
  xa=(float)fabs(wavevec-pwavea);
  if ((error==0)&&(xa>1.0e-3))
  { pwavea=wavevec; alnum=0;
    for (i=0;i<ndata;++i) if (wavevec<phasea[i*61]) break;
    if (i<0) i=0; else if (i>ndata-1) i=ndata-1;
    phaseb[0]=wavevec;
    for (j=lnum-1;j>=0;--j)
    { xb=phasea[(i-1)*61+j+1];
      xc=phasea[i*61+j+1];
      if (xb==xc) ya=xb;
      else
      { while (xc-xb>90.0) xc-=180.0f;
        while (xc-xb<-90.0) xc+=180.0f;
        ya=xb+(xc-xb)*(wavevec-phasea[(i-1)*61])/
          (phasea[i*61]-phasea[(i-1)*61]);
        while (ya>180.0) ya-=360.0f;
        while (ya<-180.0) ya+=360.0f;
      }
      if ((ya>-1.0e-5)&&(ya<1.0e-5)) ya=0.0f;
      else if (alnum<=0) alnum=j+1;
      phaseb[j+1]=ya;
    }
    palnum=alnum;
  }
  else alnum=palnum;
  return(alnum);
} //end of Phaseshift::makephase

int Phaseshift::makesinexp(float wavevec)
{ int j,alnum;
  float xa,xb,xc,xd;
  const float radian=(float)(3.14159265/180.0);

  xb=(float)fabs(wavevec-pwaveb);
  if ((error==0)&&(xb>1.0e-3))
  { pwaveb=wavevec;
    xa=(float)fabs(wavevec-pwavea);
    if (xa>1.0e-3) alnum=makephase(wavevec); else alnum=palnum;
    phasec[0]=0.0f;
    for (j=0;j<lnum;++j)
    { if (j<alnum)
      { xc=(float)sin(phaseb[j+1]*radian);
        xd=(float)cos(phaseb[j+1]*radian);
        phasec[j+1]=Fcomplex(xc*xd,xc*xc);
      }
      else phasec[j+1]=0.0f;
    }
  }
  else alnum=palnum;
  return(alnum);
} //end of Phaseshift::makesinexp

int Phaseshift::savecurve(char *filename,char *usermessage)
{ int i,j,begrow;
  float xa,xb;
  const float radian=(float)(3.14159265/180.0);
  if (error==0)
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { i=0; begrow=7;
      while (usermessage[i]!='\0')
      { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
        ++i;
      }
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(0,6); fileout.space(5);
      fileout.string("datakind beginning-row linenumbers",0,1);
      fileout.string(usermessage,0,1);
      fileout.string(symbolname,10); fileout.string(atomname,20);
      fileout.string(" phase shift data",0,2);
      fileout.string("   parameters: number of wave vectors and");
      fileout.string(" quantum momenta",0,1);
      fileout.string("   columns: k (1/angstrom)   ");
      if (lnum<1) fileout.string("no phase shift data",0,1);
      else
      { fileout.string("phase (l = 0 - ");
        fileout.integer(lnum-1); fileout.string(") (radian)",0,1);
      }
      fileout.integer(ndata,9); fileout.integer(lnum,9,1);
      xa=0.0f;
      for (i=0;i<ndata;++i)
      { for (j=0;j<lnum+1;++j)
        { if (j==0) xb=xa=phasea[i*61];
          else xb=fphase(xa,j-1)*radian;
          if (j==0) fileout.floating(xb,9.4f,256);
          else if ((j>1)&&((j%4)==1))
          { fileout.space(9,16); fileout.floating(xb,14.4f,512);
          }
          else fileout.floating(xb,14.4f,512);
        }
        fileout.newline();
      }
      if (error==0) error=fileout.geterror();
    }
    else error+=(datatype<<16);
  }

  return(error);
} //end of Phaseshift::savecurve

int Phaseshift::kconfine(float *kmin,float *kmax,int *alnum)
{ int k;
  float xa,xb;

  if (error==0)
  { xa=phasea[0];
    xb=phasea[(ndata-1)*61]*phasea[(ndata-1)*61]-
      0.512331f*0.512331f*25.0f;
    if (xb>0.0) xb=(float)sqrt(xb);
    if (*kmin<xa) *kmin=xa; if (*kmin>xb) *kmin=xb;
    if (*kmax<*kmin) *kmax=*kmin; if (*kmax>xb) *kmax=xb;
    if ((*alnum==0)||(*alnum<lnum)) *alnum=lnum;
    if (error==0) k=sizeof(int); else k=sizeof(int);
    if ((k<4)&&(*alnum>20)) *alnum=20;
  }

  return(error);
} //end of Phaseshift::kconfine

Fcomplex Phaseshift::fsinexpa(float wavevec,int lquantum)
{ Fcomplex cxa;
  if (wavevec!=pwaveb) makesinexp(wavevec);
  cxa=phasec[lquantum+1];
  return(cxa);
} //end of Phaseshift::fsinexpa

