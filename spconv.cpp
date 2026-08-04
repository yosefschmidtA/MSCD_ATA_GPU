//program to convert the potential file
#include <iostream>

#include "userinfo.h"
#include "userutil.h"
#include "xpspec.h"

int parainput(int argc,char **argv,char *finname,char *foutname,
  int *filewrite)
{ int i,k,datatype,write,error;
  char membuf[80];

  error=0; k=0;
  Textout conout;
  ++k;
  conout.string(
    "Enter name of input traditional format xps spectrum file",0,1);
  if (argc<=k) cin >> finname;
  else
  { stringcopy(finname,argv[k]);
    conout.string(finname,0,1);
  }
  ++k;
  conout.string("Enter name of output mscd format spectrum file",0,1);
  if (argc<=k) cin >> foutname;
  else
  { stringcopy(foutname,argv[k]);
    conout.string(foutname,0,1);
  }
  ++k;

  std::ifstream fin(foutname,std::ios::in);
  if (stringcomp(finname,foutname)==0) datatype=0;
  else if (fin) fin >> datatype;
  if ((fin)&&(datatype>=840)&&(datatype<850))
  { write=0;
    while (write<=0)
    { conout.string("Overwrite (O), Replace (R) or Append (A) ? ");
      if ((write<0)||(argc<=k)) cin >> membuf;
      else
      { stringcopy(membuf,argv[k]);
        conout.string(membuf,0,1);
      }
      i=0;
      while ((membuf[i]==' ')||(membuf[i]=='\t')) ++i;
      if ((membuf[i]=='O')||(membuf[i]=='o')) write=1;
      else if ((membuf[i]=='R')||(membuf[i]=='r')) write=2;
      else if ((membuf[i]=='A')||(membuf[i]=='a')) write=3;
      else write=-1;
    }
  }
  else write=1;
  *filewrite=write;

  return(error);
} //end of parainput

int main(int argc,char **argv)
{ int filewrite,error;
  char finname[80],foutname[80];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.20f,"SPCONV",
    "Yufeng Chen and Michel A Van Hove",
    "Conversion of XPS spectrum data to mscd format");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0)
    error=parainput(argc,argv,finname,foutname,&filewrite);

  if (error==0)
  { Spectrum spectrum;
    spectrum.init();
    if ((error==0)&&((filewrite==2)||(filewrite==3)))
      error=spectrum.loadcurve(foutname);
    if (error==0) error=spectrum.loadcurve(finname);
    if ((error==0)&&(filewrite==2)) error=spectrum.clean(1);
    else if (error==0) error=spectrum.clean(0);
    if (error==0)
    { conout.string("converted ",0,16);
      conout.integer(spectrum.getnumpoint());
      conout.string(" points in ");
      conout.integer(spectrum.getnumcurve());
      conout.string(" curves",0,1);
    }
    if (error==0)
      error=spectrum.savecurve(foutname,userinfo.usermsg);
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

