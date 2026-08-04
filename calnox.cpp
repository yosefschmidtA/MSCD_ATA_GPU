//program to calculate normalized chi from photoelectron diffraction chi
#include <iostream>

#include "userinfo.h"
#include "userutil.h"
#include "pdinten.h"

int parainput(int argc,char **argv,char *finname,char *foutname)
{ int error;

  error=0;
  Textout conout;
  conout.string(
    "Enter name of input photoelectron diffraction chi file",0,1);
  if (argc<=1) cin >> finname;
  else
  { stringcopy(finname,argv[1]);
    conout.string(finname,0,1);
  }

  conout.string(
    "Enter name of output photoelectron diffraction chi file",0,1);
  if (argc<=2) cin >> foutname;
  else
  { stringcopy(foutname,argv[2]);
    conout.string(foutname,0,1);
  }
  return(error);
} //end of parainput

int main(int argc,char **argv)
{ int ndata,ncurve,error;
  float unknown;
  char finname[80],foutname[80];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.20f,"CALNOX",
    "Yufeng Chen and Michel A Van Hove",
    "Calculation of normalized chi from photoelectron diffraction chi");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0) parainput(argc,argv,finname,foutname);

  if (error==0)
  { Pdintensity pdintensity;
    if (error==0) pdintensity.init();
    if (error==0) error=pdintensity.loadintensity(finname,finname,
      2,&unknown,&unknown,&unknown,&unknown,&unknown);
    if (error==0) error=pdintensity.chicalc(12);
    if (error==0)
    { ndata=pdintensity.getnumpoint();
      ncurve=pdintensity.getnumcurve();
      conout.string("calculated "); conout.integer(ndata);
      conout.string(" data points in "); conout.integer(ncurve);
      if (ncurve<2) conout.string(" curve",0,1);
      else conout.string(" curves",0,1);
    }
    if (error==0)
      error=pdintensity.savecurve(foutname,userinfo.usermsg);
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

