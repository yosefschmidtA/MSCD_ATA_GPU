//program to convert the radial matrix file
#include <iostream>

#include "userinfo.h"
#include "userutil.h"
#include "meanpath.h"

int parainput(int argc,char **argv,char *foutname,float *valence,
  float *bandgap,float *density,float *mweight)
{ int k,error;
  float xa,localeps;

  error=0; k=0; localeps=1.0e-5f;
  Textout conout;
  ++k;
  conout.string("Enter name of output mean free path file",0,1);
  if (argc<=k) cin >> foutname;
  else
  { stringcopy(foutname,argv[k]);
    conout.string(foutname,0,1);
  }

  ++k;
  conout.string("Enter number of valence electrons (>=1) or ");
  conout.string("exponent of energy (0.5-0.9): ");
  if (argc<=k) cin >> xa;
  else
  { xa=stringtofloat(argv[k]);
    if (xa>=1.0) conout.integer((int)xa,8,1);
    else conout.floating(xa,8.3f,256+1);
  }
  if (xa<0.0) *valence=0.0f;
  else if (xa>100.0) *valence=100.0f;
  else *valence=xa;

  ++k;
  if (*valence>1.0-localeps)
    conout.string("Enter bandgap energy (in eV): ");
  else
    conout.string("Enter coefficient of energy (0.02-0.3): ");
  if (argc<=k) cin >> xa;
  else
  { xa=stringtofloat(argv[k]);
    conout.floating(xa,8.4f,256+1);
  }
  if (xa<0.0) *bandgap=0.0f;
  else if (xa>100.0) *bandgap=100.0f;
  else *bandgap=xa;

  if (*valence > 1.0-localeps)
  { ++k;
    conout.string("Enter density of the bulk (in g/cm3): ");
    if (argc<=k) cin >> xa;
    else
    { xa=stringtofloat(argv[k]);
      conout.floating(xa,8.2f,256+1);
    }
    if (xa<1.0) *density=1.0f;
    else if (xa>100.0) *density=100.0f;
    else *density=xa;

    ++k;
    conout.string("Enter molecular or atomic weight (in amu): ");
    if (argc<=k) cin >> xa;
    else
    { xa=stringtofloat(argv[k]);
      conout.floating(xa,8.1f,1);
    }
    if (xa<1.0) *mweight=1.0f;
    else if (xa>1000.0) *mweight=100.0f;
    else *mweight=xa;
  }
  else *density=*mweight=0.0f;
  return(error);
} //end of parainput

int main(int argc,char **argv)
{ int ndata,error;
  float valence,bandgap,density,mweight;
  char foutname[80];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.20f,"CALMFP","Yufeng Chen and Michel A Van Hove",
    "Calculation of correlated thermal vibrational displacement");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0) parainput(argc,argv,foutname,&valence,&bandgap,
    &density,&mweight);

  if (error==0)
  { Meanpath meanpath;
    meanpath.init();
    if (error==0)
      error=meanpath.loadparameter(valence,bandgap,density,mweight);
    if (error==0) error=meanpath.savecurve(foutname,userinfo.usermsg);
    if (error==0)
    { ndata=meanpath.getnumpoint();
      conout.string("calculated "); conout.integer(ndata);
      conout.string(" data points",0,1);
    }
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

