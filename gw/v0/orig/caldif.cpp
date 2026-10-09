//program to calculate reliability between two photoelectron diffraction
//  chi data curves
#include <iostream>
#include <iomanip>

#include "userinfo.h"
#include "userutil.h"
#include "pdinten.h"

int parainput(int argc,char **argv,char *fina,char *finb)
{ int error;

  error=0;
  Textout conout;
  conout.string(
    "Enter name of first photoelectron diffraction chi file",0,1);
  if (argc<=1) cin >> fina;
  else
  { stringcopy(fina,argv[1]);
    conout.string(fina,0,1);
  }

  conout.string(
    "Enter name of second photoelectron diffraction chi file",0,1);
  if (argc<=2) cin >> finb;
  else
  { stringcopy(finb,argv[2]);
    conout.string(finb,0,1);
  }
  return(error);
} //end of parainput

float dispreliable(float areliable,float breliable,
  float creliable,float dreliable,
  float ereliable,float freliable,
  float greliable,float hreliable)
{ Textout conout;
  conout.string("reliability factors of two photoelectron diffraction",
    0,16);
  conout.string(" chi data",0,1);
  conout.string("   afac = ");
  conout.floating(areliable,9.4f,256);
  conout.string("      bfac = ");
  conout.floating(breliable,9.4f,256+2);
  conout.string("   afac=sum((chi1-chi2)*(chi1-chi2))/");
  conout.string("sum(chi1*chi1+chi2*chi2)",0,1);
  conout.string("   bfac=sum(chi1*chi1-chi2*chi2)/");
  conout.string("sum(chi1*chi1+chi2*chi2)",0,1);

  conout.string("Van Hove and Fadley group's reliability factors",
    0,16+1);
  conout.string("   rfa1 = ");
  conout.floating(creliable,9.4f,256);
  conout.string("      rfa2 = ");
  conout.floating(dreliable,9.4f,256+1);
  conout.string("   rfa4 = ");
  conout.floating(ereliable,9.4f,256);
  conout.string("      rfa5 = ");
  conout.floating(freliable,9.4f,256+1);
  conout.string("   rfa3 = ");
  conout.floating(greliable,9.4f,256);
  conout.string("      rfa6 = ");
  conout.floating(hreliable,9.4f,256+2);
  conout.string("   rfa1=sum(abs(chical-chiexp))/");
  conout.string("sum(abs(chiexp))",0,1);
  conout.string("   rfa2=sum((chical-chiexp)*(chical-chiexp))/");
  conout.string("sum(chiexp*chiexp)",0,1);
  conout.string("   rfa4=sum(chical'-chiexp')/");
  conout.string("sum(chiexp')",0,1);
  conout.string("   rfa5=sum((chical'-chiexp')*(chical'-chiexp'))/");
  conout.string("sum(chiexp'*chiexp')",0,1);
  conout.string("   rfa3=percentage of angle range over which chical");
  conout.string("and calexp",0,1);
  conout.string("        have slopes of different sign(+/-)",0,1);
  conout.string("   rfa6=(rfac1+rfac2+rfac4+rfac5+rfac3)/");
  conout.string("5",0,1);

  return(0.0f);
} //end of dispreliable

int main(int argc,char **argv)
{ int ndata,ncurve,error;
  float amin,amax,astep,areliable,breliable,creliable,dreliable,
    ereliable,freliable,greliable,hreliable,unknown;
  char fina[80],finb[80];

  error=0;
  Userinfo userinfo; userinfo.init();
  userinfo.loadparameter(1.21f,"CALDIF",
    "Yufeng Chen and Michel A Van Hove",
    "Calculation of reliability of photoelectron diffraction chi data");
  Textout conout;
  conout.string(userinfo.usermsg,0,1);
  if (userinfo.version==0.0) error=102;
  if (error==0) parainput(argc,argv,fina,finb);

  amin=0.0f; amax=1000.0f; astep=0.01f;
  if (error==0)
  { Pdintensity pdintensity;
    pdintensity.init();
    if (error==0) error=pdintensity.loadintensity(fina,finb,
      14,&amin,&amax,&astep,&unknown,&unknown);
    if (error==0)
      error=pdintensity.reliability(&areliable,&breliable,
        &creliable,&dreliable,&ereliable,&freliable,&greliable,
        &hreliable);
    if (error==0)
    { ndata=pdintensity.getnumpoint();
      ncurve=pdintensity.getnumcurve();
      conout.string("calculated "); conout.integer(ndata);
      conout.string(" data points in "); conout.integer(ncurve);
      if (ncurve<2) conout.string(" curve",0,1);
      else conout.string(" curves",0,1);
    }
    if (error==0) dispreliable(areliable,breliable,creliable,
      dreliable,ereliable,freliable,greliable,hreliable);
  }

  Errorinfo errorinfo(error);
  conout.string(errorinfo.errmsg,0,16+1);

  return(error);
} //end of main

