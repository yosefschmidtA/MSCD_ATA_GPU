#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "userutil.h"
#include "radmat.h"

Radialmatrix::Radialmatrix()
{ error=107; init();
}//end of Radialmatrix::Radialmatrix

Radialmatrix::~Radialmatrix()
{ if (atomname) delete [] atomname;
  if (symbolname) delete [] symbolname;
  if (subshell) delete [] subshell; if (radmat) delete [] radmat;
}//end of Radialmatrix::~Radialmatrix

void Radialmatrix::init()
{ if ((error==107)||(error==103))
  { datatype=721;
    atomname=symbolname=subshell=NULL; radmat=NULL;
  }
  if (error==103)
  { atomname=new char [25]; symbolname=new char [10];
    subshell=new char [10]; radmat=new float [256*5];
    if ((!atomname)||(!symbolname)||(!subshell)||(!radmat)) error=102;
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
    if (subshell)
    { subshell=new char [10];
      if (!subshell) error=102;
    }
    if (radmat)
    { radmat=new float [256*5];
      if (!radmat) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Radialmatrix::init

int Radialmatrix::getlength()
{ int ka;
  if (error==0) ka=sizeof(Radialmatrix)+(25+10+10)*sizeof(char)+
    256*5*sizeof(float)+sizeof(int);
  else ka=0;
  return(ka);
} //end of Radialmatrix::getlength

int Radialmatrix::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Radialmatrix);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=25;
    if (ka>=0) ka=memorysend(dest,atomname,ka,kb,length);
    kb=10;
    if (ka>=0) ka=memorysend(dest,symbolname,ka,kb,length);
    kb=10;
    if (ka>=0) ka=memorysend(dest,subshell,ka,kb,length);
    kb=256*5*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)radmat,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Radialmatrix::paexport

int Radialmatrix::paimport(char *source,int length)
{ int ka,kb;

  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Radialmatrix);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  if (error==0)
  { kb=25;
    if (ka>=0) ka=memoryrec(atomname,source,ka,kb,length);
    kb=10;
    if (ka>=0) ka=memoryrec(symbolname,source,ka,kb,length);
    kb=10;
    if (ka>=0) ka=memoryrec(subshell,source,ka,kb,length);
    kb=256*5*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)radmat,source,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Radialmatrix::paimport

int Radialmatrix::getnumpoint()
{ int k;
  if (error==0) k=ndata;
  else k=0;
  return(k);
} //end of Radialmatrix::getnumpoint

int Radialmatrix::getnumfinal()
{ int k;
  if ((error==0)&&((subshell[1]=='s')||(subshell[1]=='S'))) k=2;
  else if (error==0) k=4;
  else k=0;
  return(k);
} //end of Radialmatrix::getnumfinal

float Radialmatrix::getmemory()
{ float xa;
  if (error==0) xa=(float)(sizeof(Radialmatrix)+
    (25+10+10)*sizeof(char)+256*5*sizeof(float));
  else xa=-1.0f;
  return(xa);
} //end of Radialmatrix::getmemory

int Radialmatrix::loadcurve(char *filename,char *iatomname,
  char *isymbolname,char *isubshell)
{ int i,j,k,begrow,npoint;
  float xa,xb,xc;
  char cha,chb;
  const float radian=(float)(3.14159265/180.0);

  i=begrow=npoint=0;
  if ((error==101)||(error==0))
  { error=0;
    i=0; while (iatomname[i]==' ') ++i;
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
    i=0; while (isubshell[i]==' ') ++i;
    cha=isubshell[i];
    ++i; while (isubshell[i]==' ') ++i;
    chb=isubshell[i];
    subshell[0]=cha;
    if ((chb=='s')||(chb=='S'))
    { i=0; subshell[1]='s';
    }
    else if ((chb=='p')||(chb=='P'))
    { i=1; subshell[1]='p';
    }
    else if ((chb=='d')||(chb=='D'))
    { i=2; subshell[1]='d';
    }
    else if ((chb=='f')||(chb=='F'))
    { i=3; subshell[1]='f';
    }
    else if ((chb=='g')||(chb=='G'))
    { i=4; subshell[1]='g';
    }
    else error=623;
    j=(int)(cha-'0');
    if ((error==0) && ((j<=i)||(j>6))) error=623;
    subshell[2]='\0';
  }
  if (error==0)
  { std::ifstream fin(filename,std::ios::in);
    if (!fin) error=201+(datatype<<16);
    else
    { fin >> xa; k=(int)xa;
      if (fin && (k==datatype))
      { fin >> xb >> xc; skipendline(fin);
        begrow=(int)xb; npoint=(int)xc;
        if (npoint <0) npoint=0;
        else if (npoint>256) npoint=256;
        for (i=1;i<begrow-1;++i) skipendline(fin);
        for(i=0;fin&&(i<npoint);++i)
        { fin >> radmat[i*5+0] >> radmat[i*5+1] >> radmat[i*5+2];
          if ((subshell[1] != 's') && (subshell[1] != 'S'))
            fin >> radmat[i*5+3] >> radmat[i*5+4];
          skipendline(fin);
        }
        if (!fin) error=204+(datatype<<16);
      }
      else if (fin && (xa > 1.0e5)) error=209+(datatype<<16);
      else
      { i=0;
        while (fin && (i<256))
        { skipendline(fin);
          if ((subshell[1] != 's') && (subshell[1] != 'S'))
            fin >> radmat[i*5+0] >> radmat[i*5+1] >> radmat[i*5+3]
              >> radmat[i*5+2] >> radmat[i*5+4];
          else
          fin >> radmat[i*5+0] >> radmat[i*5+1] >> radmat[i*5+2];
          if (fin&&(radmat[i*5+0]>1.0)) ++i;
          else break;
        }
        npoint=i;
        for (i=0;(i<npoint)&&(error==0);++i)
        { if (radmat[i*5+0]>xa)
            radmat[i*5+0]=0.512331f*(float)sqrt(radmat[i*5+0]-xa);
          else error=717+(datatype<<16);
        }
      }
    }
  }
  for (i=0;(i<npoint-2)&&(error==0);++i)
  { if (radmat[i*5+5]-radmat[i*5]<1.0e-5) error=707+(datatype<<16);
    radmat[i*5+2]/=radian;
    if ((subshell[1]!='s')&&(subshell[1]!='S')) radmat[i*5+4]/=radian;
  }
  if ((error==0)&&(npoint>0)) ndata=npoint;
  else if (error==0) error=204+(datatype<<16);
  return(error);
} //end of Radialmatrix::loadcurve

float Radialmatrix::famphase(float wavevec,int code)
{ int i;
  float xa,xb,ya;

  if ((error==0)&&((code<0)||(code>4))) ya=0.0f;
  else if (error==0)
  { for (i=0;i<ndata;++i) if (wavevec<radmat[i*5]) break;
    if (i<=0) ya=radmat[code];
    else if (i>ndata-1) ya=radmat[(ndata-1)*5+code];
    else
    { xa=radmat[(i-1)*5+code];
      xb=radmat[i*5+code];
      if ((code==2)||(code==4))
      { while (xb-xa>180.0) xb-=360.0f;
        while (xb-xa<-180.0) xb+=360.0f;
      }
      ya=xa+(xb-xa)*(wavevec-radmat[(i-1)*5])/
        (radmat[i*5]-radmat[(i-1)*5]);
      if ((code==2)||(code==4))
      { while (ya>180.0) ya-=360.0f;
        while (ya<-180.0) ya+=360.0f;
      }
    }
  }
  else ya=0.0f;

  return(ya);
}//end of Radialmatrix::famphase

int Radialmatrix::savecurve(char *filename,char *usermessage)
{ int i,begrow;
  float xa;
  const float radian=(float)(3.14159265/180.0);

  if (error==0)
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { i=0; begrow=6;
      while (usermessage[i]!='\0')
      { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
        ++i;
      }
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(ndata,6); fileout.space(5);
      fileout.string("datakind beginning-row linenumbers",0,1);
      fileout.string(usermessage,0,1); fileout.string(symbolname,10);
      fileout.space(1); fileout.charfill(subshell[0],1);
      fileout.charfill(subshell[1],1);
      fileout.string(atomname,18);
      fileout.string(" radial matrix data",0,2);
      fileout.string("k",12); fileout.string("R(li+1)",15);
      fileout.string("phase(li+1)",12);
      if ((subshell[1]!='s') && (subshell[1]!='S'))
      { fileout.string("R(li-1)",15);
        fileout.string("phase(li-1)",12,1);
      }
      else fileout.newline();
      for (i=0;i<ndata;++i)
      { xa=radmat[i*5];
        fileout.floating(xa,12.4f,256);
        fileout.floating(famphase(xa,1),15.4f,512);
        fileout.floating(famphase(xa,2)*radian,15.4f,512);
        if ((subshell[1]!='s')&&(subshell[1]!='S'))
        { fileout.floating(famphase(xa,3),15.4f,512);
          fileout.floating(famphase(xa,4)*radian,15.4f,512);
        }
        fileout.newline();
      }
      if (error==0) error=fileout.geterror();
    }
    else error+=(datatype<<16);
  }
  return(error);
} //end of Radialmatrix::savecurve

int Radialmatrix::kconfine(float *kmin,float *kmax)
{ float xa,xb;
  if (error==0)
  { xa=radmat[0];
    xb=radmat[(ndata-1)*5]*radmat[(ndata-1)*5]-
      0.512331f*0.512331f*25.0f;
    if (xb>0.0) xb=(float)sqrt(xb);
    if (*kmin<xa) *kmin=xa; if (*kmin>xb) *kmin=xb;
    if (*kmax<*kmin) *kmax=*kmin; if (*kmax>xb) *kmax=xb;
  }

  return(error);
} //end of Radialmatrix::kconfine

