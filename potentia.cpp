#include <iostream>
#include <fstream>
#include <iomanip>

#include "userutil.h"
#include "potentia.h"

Potential::Potential()
{ error=107; init();
} //end of Potential::Potential

Potential::~Potential()
{ if (atomname) delete [] atomname;
  if (symbolname) delete [] symbolname;
  if (radist) delete [] radist; if (radpot) delete [] radpot;
} //end of Potential::~Potential

void Potential::init()
{ if ((error==107)||(error==103))
  { datatype=811; atomname=symbolname=NULL; radist=radpot=NULL;
  }
  if (error==103)
  { atomname=new char [25]; symbolname=new char [10];
    radist=new float [256]; radpot=new float [256];
    if ((!atomname)||(!symbolname)||(!radist)||(!radpot)) error=102;
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
    if (radist)
    { radist=new float [256];
      if (!radist) error=102;
    }
    if (radpot)
    { radpot=new float [256];
      if (!radpot) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Potential::init

int Potential::getlength()
{ int ka;
  if (error==0) ka=sizeof(int)+sizeof(Potential)+
    (25+10)*sizeof(char)+2*256*sizeof(float);
  else ka=0;
  return(ka);
} //end of Potential::getlength

int Potential::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Potential);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=25;
    if (ka>=0) ka=memorysend(dest,(char *)atomname,ka,kb,length);
    kb=10;
    if (ka>=0) ka=memorysend(dest,(char *)symbolname,ka,kb,length);
    kb=256*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)radist,ka,kb,length);
    kb=256*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)radpot,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Potential::paexport

int Potential::paimport(char *source,int length)
{ int ka,kb;

  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Potential);
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
    if (ka>=0) ka=memoryrec((char *)atomname,source,ka,kb,length);
    kb=10;
    if (ka>=0) ka=memoryrec((char *)symbolname,source,ka,kb,length);
    kb=256*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)radist,source,ka,kb,length);
    kb=256*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)radpot,source,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Potential::paimport

int Potential::getnumpoint()
{ int k;
  if (error==0) k=ndata;
  else k=0;
  return(k);
} //end of Potential::getnumpoint

float Potential::getmemory()
{ float xa;
  if (error==0) xa=(float)(sizeof(Potential)+291*sizeof(char)+
    512*sizeof(int));
  else xa=-1.0f;
  return(xa);
} //end of Potential::getmemory

int Potential::loadcurve(char *filename,char *iatomname,
  char *isymbolname)
{ int i,j,k,begrow,npoint;
  float xa;

  if ((error==101)||(error==0))
  { error=0; i=begrow=npoint=0;
    std::ifstream fin(filename,std::ios::in);
    if (!fin) error=201+(datatype<<16);
    else
    { fin >> radist[0] >> radpot[0];
      k=(int)radist[0];
      if (fin && (k==datatype))
      { begrow=(int)radpot[0];
        fin >> xa; skipendline(fin);
        npoint=(int)xa;
        if (npoint<0) npoint=0;
        else if (npoint>256) npoint=256;
        for (j=1;j<begrow-1;++j) skipendline(fin);
        if (fin)
        { fin >> radist[0] >> radpot[0];
          skipendline(fin);
        }
        if (!fin) error=204+(datatype<<16);
      }
      else if (fin && (k > 100)) error=209+(datatype<<16);
      else if (!fin) error=204+(datatype<<16);
      else skipendline(fin);

      if ((error==0) && fin) i=1;
      while (fin && (error==0) && ((i<256) || ((npoint>0)&&(i<npoint))))
      { fin >> radist[i] >> radpot[i];
        skipendline(fin);
        if (fin) ++i;
      }
      if ((error==0)&&(npoint==0))
      { for (j=0;j<i;++j)
        { radist[j]*=0.529167f; radpot[j]*=13.605f*0.529167f;
        }
        if (radist[0]<1.0e-10) radist[0]=0.0f;
      }
      if (error==0) npoint=i;
    }
    if (error==0)
    { for (i=0;i<npoint-2;++i)
        if (radist[i+1]-radist[i]<1.0e-10) error=707+(datatype<<16);
    }
    if ((error==0)&&(npoint>0)) ndata=npoint;
    else if (error==0) error=204+(datatype<<16);
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
  return(error);
} //end of Potential::loadcurve

float Potential::fpoten(float distance)
{ int i;
  float poten;
  if (error==0)
  { for (i=0; i<ndata; ++i) if (distance<radist[i]) break;
    if (i<=0) poten=radpot[0];
    else if (i>ndata-1) poten=radpot[ndata-1];
    else poten=radpot[i-1]+(radpot[i]-radpot[i-1])*
      (distance-radist[i-1])/(radist[i]-radist[i-1]);
  }
  else poten=0.0f;
  return(poten);
} //end of Potential::fpoten

int Potential::savecurve(char *filename,char *usermessage)
{ int i,begrow;
  float xa;
  
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
      fileout.string(usermessage,0,1);
      fileout.string(symbolname,10); fileout.string(atomname,20);
      fileout.string(" potential data",0,2);
      fileout.string("r (angstrom)",15);
      fileout.string("rV (eV-angs)",15,1);
      for (i=0;i<ndata;++i)
      { xa=radist[i];
        fileout.floating(xa,15.6f,512);
        fileout.floating(fpoten(xa),15.6f,512+1);
      }
      if (error==0) error=fileout.geterror();
    }
    else error+=(datatype<<16);
  }
  return(error);
} //end of Potential::savecurve

