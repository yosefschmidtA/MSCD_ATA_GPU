#include <iostream>

#include "userinfo.h"
#include "userutil.h"
#include "vibrate.h"

int parainput(int argc,char **argv,char *foutname,float *density,
  float *mweight,float *aweight,float *tdebye,float *tsample)
{ int error;
  float xa;

  error=0;
  Textout conout;
  conout.string(
    "Enter name of output vibrational displacement file",0,1);
  if (argc<=1) cin >> foutname;
  else
  { stringcopy(foutname,argv[1]);
    conout.string(foutname,0,1);
  }
  conout.string("Enter density of the bulk (in g/cm3): ");
  if (argc<=2) cin >> xa;
  else
  { xa=stringtofloat(argv[2]);
    conout.floating(xa,8.2f,1);
  }
  if (xa<1.0) *density=1.0f;
  else if (xa>100.0) *density=100.0f;
  else *density=xa;

  conout.string("Enter molecular or atomic weight (in amu): ");
  if (argc<=3) cin >> xa;
  else
  { xa=stringtofloat(argv[3]);
    conout.floating(xa,8.1f,1);
  }
  if (xa<1.0) *mweight=1.0f;
  else if (xa>1000.0) *mweight=100.0f;
  else *mweight=xa;

  conout.string("Enter atomic weight for specified atom (in amu): ");
  if (argc<=4) cin >> xa;
  else
  { xa=stringtofloat(argv[4]);
    conout.floating(xa,8.1f,1);
  }
  if (xa<1.0) *aweight=1.0f;
  else if (xa>1000.0) *aweight=1000.0f;
  else *aweight=xa;

  conout.string("Enter debye temperature (in K): ");
  if (argc<=5) cin >> xa;
  else
  { xa=stringtofloat(argv[5]);
    conout.floating(xa,8.1f,1);
  }
  if (xa<0.0) *tdebye=0.0f;
  else if (xa>1000.0) *tdebye=1000.0f;
  else *tdebye=xa;

  conout.string("Enter sample temperature (in K): ");
  if (argc<=6) cin >> xa;
  else
  { xa=stringtofloat(argv[6]);
    conout.floating(xa,8.1f,1);
  }
  if (xa<0.0) *tsample=0.0f;
  else if (xa>1000.0) *tsample=1000.0f;
  else *tsample=xa;

  return(error);
} //end of parainput

int main(int argc,char **argv)
{ int ndata,error;
  float density,mweight,aweight,tdebye,tsample;
  char foutname[80];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.20f,"CALVIB",
    "Yufeng Chen and Michel A Van Hove",
    "Calculation of correlated thermal vibrational displacement");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0) parainput(argc,argv,foutname,&density,&mweight,&aweight,
    &tdebye,&tsample);
  if (error==0)
  { Vibration vibration;
    vibration.init();
    if (error==0)
      error=vibration.loadparameter(density,mweight,tdebye,tsample);
    if (error==0)
      error=vibration.savecurve(foutname,userinfo.usermsg,aweight);
    if (error==0)
    { ndata=vibration.getnumpoint();
      conout.string("calculated "); conout.integer(ndata);
      conout.string(" data points",0,1);
    }
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

