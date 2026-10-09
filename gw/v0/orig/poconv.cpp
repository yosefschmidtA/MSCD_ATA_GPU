//program to convert the potential file
#include <iostream>

#include "userinfo.h"
#include "userutil.h"
#include "potentia.h"

int parainput(int argc,char **argv,char *finname,char *foutname,
  char *atomname,char *symbolname)
{ int error;

  error=0;
  Textout conout;
  conout.string(
    "Enter name of input traditional format potential file",0,1);
  if (argc<=1) cin >> finname;
  else
  { stringcopy(finname,argv[1]);
    conout.string(finname,0,1);
  }
  conout.string("Enter name of output mscd format potential file",0,1);
  if (argc<=2) cin >> foutname;
  else
  { stringcopy(foutname,argv[2]);
    conout.string(foutname,0,1);
  }
  cout << "Enter name of the atom (e.g. Copper): ";
  if (argc<=3) cin >> atomname;
  else
  { stringcopy(atomname,argv[3]);
    conout.string(atomname,0,1);
  }
  cout << "Enter symbol of the atom (e.g. Cu): ";
  if (argc<=4) cin >> symbolname;
  else
  { stringcopy(symbolname,argv[4]);
    conout.string(symbolname,0,1);
  }

  return(error);
} //end of parainput

int main(int argc,char **argv)
{ int error;
  char finname[80],foutname[80],atomname[25],symbolname[10];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.20f,"POCONV",
    "Yufeng Chen and Michel A Van Hove",
    "Conversion of traditional potential data to mscd format");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0)
    error=parainput(argc,argv,finname,foutname,atomname,symbolname);

  if (error==0)
  { Potential potential;
    potential.init();
    if (error==0) 
      error=potential.loadcurve(finname,atomname,symbolname);
    if (error==0)
    { conout.string("calculated ");
      conout.integer(potential.getnumpoint());
      conout.string(" points",0,1);
    }
    if (error==0)
      error=potential.savecurve(foutname,userinfo.usermsg);
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

