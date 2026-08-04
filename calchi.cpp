//program to calculate chi from photoelectron diffraction intensity
#include <iostream>

#include "userinfo.h"
#include "userutil.h"
#include "pdinten.h"

int parainput(int argc,char **argv,char *finname,char *foutname)
{ int error;

  error=0;
  Textout conout;
  conout.string("Enter name of input photoelectron");
  conout.string(" diffraction intensity file",0,1);
  if (argc<=1) cin >> finname;
  else
  { stringcopy(finname,argv[1],40);
    conout.string(finname,0,1);
  }

  conout.string("Enter name of output photoelectron diffraction ");
  conout.string("intensity file",0,1);
  if (argc<=2) cin >> foutname;
  else
  { stringcopy(foutname,argv[2],40);
    conout.string(foutname,0,1);
  }
  return(error);
} //end of parainput

int main(int argc,char **argv)
{ int k,ndata,ncurve,error;
  float unknown;
  char finname[80],foutname[80];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.20f,"CALCHI",
    "Yufeng Chen and Michel A Van Hove",
    "Calculation of chi from photoelectron diffraction intensity");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0) parainput(argc,argv,finname,foutname);

  if (error==0)
  { Pdintensity pdintensity;
    pdintensity.init();
    k=stringcomp(foutname,"transfer",8);
    if (k==0)
    { error=pdintensity.loadintensity(finname,finname,
        2,&unknown,&unknown,&unknown,&unknown,&unknown);
    }
    else if (error==0)
    { error=pdintensity.loadintensity(finname,finname,
        1,&unknown,&unknown,&unknown,&unknown,&unknown);
      if (error==0) error=pdintensity.chicalc(11);
    }
    if (error==0)
    { ndata=pdintensity.getnumpoint(); ncurve=pdintensity.getnumcurve();
      if (k!=0) conout.string("calculated ");
      else conout.string("transferred ");
      conout.integer(ndata);
      conout.string(" data points in "); conout.integer(ncurve);
      if (ncurve<2) conout.string(" curve",0,1);
      else conout.string(" curves",0,1);
    }
    if ((error==0)&&(k!=0))
      error=pdintensity.savecurve(foutname,userinfo.usermsg);
    else if (error==0)
      error=pdintensity.savecurve(finname,userinfo.usermsg);
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

