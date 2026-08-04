//program to calculate scattering factors
#include <iostream>

#include "userinfo.h"
#include "userutil.h"
#include "scatter.h"

int parainput(int argc,char **argv,char *finname,char *foutname,
  char *atomname,char *symbolname,float *akin,float *lenga,
  float *lengb,int *lnum,int *raorder)
{ int k,error;
  float xa;

  error=0; k=0;
  Textout conout;
  if (error==0)
  { conout.string("Enter name of input phase shift file",0,1);
    ++k;
    if (argc<=k) cin >> finname;
    else
    { stringcopy(finname,argv[k]);
      conout.string(finname,0,1);
    }
  }

  if (error==0)
  { conout.string("Enter name of output scattering factor file",0,1);
    ++k;
    if (argc<=k) cin >> foutname;
    else
    { stringcopy(foutname,argv[k]);
      conout.string(foutname,0,1);
    }
  }

  if (error==0)
  { conout.string("Enter name of the atom (e.g. Copper): ");
    ++k;
    if (argc<=k) cin >> atomname;
    else
    { stringcopy(atomname,argv[k]);
      conout.string(atomname,0,1);
    }
  }

  if (error==0)
  { conout.string("Enter symbol of the atom (e.g. Cu): ");
    ++k;
    if (argc<=4) cin >> symbolname;
    else
    { stringcopy(symbolname,argv[k]);
      conout.string(symbolname,0,1);
    }
  }

  if (error==0)
  { conout.string("Enter wave vector in inverse angstrom ");
    conout.string("(3.0-25.0): ");
    ++k;
    if (argc<=k) cin >> xa;
    else
    { xa=stringtofloat(argv[k]);
      conout.floating(xa,8.1f,1);
    }
    if (xa<3.0) *akin=3.0f;
    else if (xa>25.0) *akin=25.0f;
    else *akin=xa;
  }

  if (error==0)
  { conout.string("Enter first bonding length in angstrom ");
    conout.string("(1.0-100.0): ");
    ++k;
    if (argc<=k) cin >> xa;
    else
    { xa=stringtofloat(argv[k]);
      conout.floating(xa,8.3f,1);
    }
    if (xa<1.0) *lenga=1.0f;
    else if (xa>100.0) *lenga=100.0f;
    else *lenga=xa;
  }

  if (error==0)
  { conout.string("Enter second bonding length in angstrom ");
    conout.string("(1.0-100.0): ");
    ++k;
    if (argc<=k) cin >> xa;
    else
    { xa=stringtofloat(argv[k]);
      conout.floating(xa,8.3f,1);
    }
    if (xa<1.0) *lengb=1.0f;
    else if (xa>100.0) *lengb=100.0f;
    else *lengb=xa;
  }

  if (error==0)
  { conout.string("Enter maximum number of quantum momentum lnum ");
    conout.string("(0-60): ");
    ++k;
    if (argc<=k) cin >> xa;
    else
    { xa=stringtofloat(argv[k]);
      conout.integer((int)xa,4,1);
    }
    if (xa<0) *lnum=0;
    else if (xa>60) *lnum=60;
    else *lnum=(int)xa;
  }

  if (error==0)
  { conout.string("Enter Rehr-Albers approximation order ");
    conout.string("(0-4): ");
    ++k;
    if (argc<=k) cin >> xa;
    else
    { xa=stringtofloat(argv[k]);
      conout.integer((int)xa,4,1);
    }
    if (xa<0) *raorder=0;
    else if (xa>4) *raorder=4;
    else *raorder=(int)xa;
  }

  return(error);
} //end of parainput

int main(int argc,char **argv)
{ int ndata,lnum,raorder,error;
  float akin,lenga,lengb;
  char finname[80],foutname[80],atomname[80],symbolname[80];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.20f,"CALFAC",
    "Yufeng Chen and Michel A Van Hove",
    "Calculation of photoelectron scattering factors");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0) parainput(argc,argv,finname,foutname,atomname,
    symbolname,&akin,&lenga,&lengb,&lnum,&raorder);

  if (error==0)
  { Scatter scatter;
    scatter.init(lnum,raorder);
    if (error==0)
      error=scatter.loadcurve(finname,atomname,symbolname);
    if (error==0)
    { ndata=scatter.getnumpoint();
      conout.string("calculating ",0,16); conout.integer(ndata,8);
      conout.string(" data points ...",0,1);
    }
    if (error==0) error=scatter.makefactor(akin,lenga,lengb);
    if (error==0) error=scatter.savecurve(foutname,userinfo.usermsg);
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

