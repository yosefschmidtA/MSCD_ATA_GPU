//program to convert the radial matrix file
#include <iostream>

#include "userinfo.h"
#include "userutil.h"
#include "radmat.h"

int parainput(int argc,char **argv,char *finname,char *foutname,
  char *atomname,char *symbolname,char *subshell)
{ int error;

  error=0;
  Textout conout;
  conout.string(
    "Enter name of input traditional format radial matrix file",0,1);
  if (argc<=1) cin >> finname;
  else
  { stringcopy(finname,argv[1]);
    conout.string(finname,0,1);
  }
  conout.string(
    "Enter name of output mscd format radial matrix file",0,1);
  if (argc<=2) cin >> foutname;
  else
  { stringcopy(foutname,argv[2]);
    conout.string(foutname,0,1);
  }
  conout.string("Enter name of the atom (e.g. Copper): ");
  if (argc<=3) cin >> atomname;
  else
  { stringcopy(atomname,argv[3]);
    conout.string(atomname,0,1);
  }
  conout.string("Enter symbol of the atom (e.g. Cu): ");
  if (argc<=4) cin >> symbolname;
  else
  { stringcopy(symbolname,argv[4]);
    conout.string(symbolname,0,1);
  }
  conout.string("Enter name of the subshell (e.g. 3p): ");
  if (argc<=5) cin >> subshell;
  else
  { stringcopy(subshell,argv[5]);
    conout.string(subshell,0,1);
  }

  return(error);
} //end of parainput

int main(int argc,char **argv)
{ int i,j,k,error;
  char finname[80],foutname[80],atomname[25],symbolname[10],
    subshell[10];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.20f,"RMCONV",
    "Yufeng Chen and Michel A Van Hove",
    "Conversion of traditional radial matrix file to mscd format");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0) error=
    parainput(argc,argv,finname,foutname,atomname,symbolname,subshell);

  if (error==0)
  { Radialmatrix radmatrix; radmatrix.init();
    error=radmatrix.loadcurve(finname,atomname,symbolname,subshell);
    if (error==0)
    { i=radmatrix.getnumfinal(); j=radmatrix.getnumpoint();
      k=i*j;
      conout.string("calculated "); conout.integer(k);
      conout.string(" elements for "); conout.integer(j);
      conout.string(" energies",0,1);
    }
    if (error==0)
      error=radmatrix.savecurve(foutname,userinfo.usermsg);
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

