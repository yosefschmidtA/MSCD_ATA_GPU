//program to convert the potential file
#include <iostream>

#include "userinfo.h"
#include "userutil.h"
#include "xpspec.h"

int parainput(int argc,char **argv,char *finname,char *foutname,
  int *npeak,float *epeak)
{ int j,k,m,error;
  float xa;

  error=0; k=0; xa=0.0f;
  Textout conout;
  ++k;
  conout.string(
    "Enter name of input mscd format xps spectrum file",0,1);
  if (argc<=k) cin >> finname;
  else
  { stringcopy(finname,argv[k]);
    conout.string(finname,0,1);
  }
  ++k;
  conout.string("Enter name of output xps peak file",0,1);
  if (argc<=k) cin >> foutname;
  else
  { stringcopy(foutname,argv[k]);
    conout.string(foutname,0,1);
  }
  ++k;
  conout.string("Enter nnumber of peaks (1-3) : ");
  if (argc<=k) cin >> xa;
  else
  { xa=stringtofloat(argv[k]);
    conout.integer((int)xa,4,1);
  }
  m=(int)xa;
  if (m<0) m=0;
  else if (m>3) m=3;
  ++k;
  if (m>1)
  { for (j=0;j<m;++j)
    { conout.string("Enter binding energy of peak ");
      conout.integer(j+1);
      conout.string(" : ");
      if (argc<=k) cin >> xa;
      else
      { xa=stringtofloat(argv[k+j]);
        conout.floating(xa,8.2f,1);
      }
      if (xa<0.0) m=j;
      else epeak[j]=xa;
    }
  }
  else if (m>0) epeak[0]=0.0f;
  *npeak=m;

  return(error);
} //end of parainput

int main(int argc,char **argv)
{ int npeak,error;
  char finname[80],foutname[80];
  float epeak[3];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.25f,"XPSPEAK",
    "Yufeng Chen and Michel A Van Hove",
    "XPS spectrum peakfit positions and intensities");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0)
    error=parainput(argc,argv,finname,foutname,&npeak,epeak);

  if (error==0)
  { Spectrum spectrum;
    spectrum.init();
    if (error==0) error=spectrum.loadcurve(finname);
    if (error==0) error=spectrum.clean(1);
    if (error==0) error=spectrum.peakfit(npeak,epeak,1,1);
    if (error==0)
    { conout.string("processed ",0,16);
      conout.integer(spectrum.getnumpoint());
      conout.string(" points in ");
      conout.integer(spectrum.getnumcurve());
      conout.string(" curves",0,1);
    }
    if (error==0) error=spectrum.savepeak(foutname,userinfo.usermsg);
    if (error==0)
    { conout.string("writing fitted curves into file xpslist.txt ...",
        0,1);
      error=spectrum.savelist("xpslist.txt",userinfo.usermsg);
    }
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

