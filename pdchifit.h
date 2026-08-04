/*
----------------------------------------------------------------------
  C++ classes designed for photoelectron diffraction software package
  Yufeng Chen, Michel A. Van Hove
  Lawrence Berkeley National Laboratory (LBNL), Berkeley, CA 94720
  Copyright (c) Van Hove Group 1997-1998. All rights reserved.
----------------------------------------------------------------------

  Yufeng Chen, LBNL MS 2-100, 1 Cyclotron Road, Berkeley, CA 94720
  Phone (510) 486-4581 Fax (510) 486-5530 Email. ychen@LBL.gov

  M. A. Van Hove, LBNL MS 66-100, 1 Cyclotron Road, Berkeley, CA 94720
  Phone (510) 486-6160 Fax (510) 486-4995 Email. MAVanhove@LBL.gov
----------------------------------------------------------------------

List of classes

1. Pdchifit
Description: curve fitting routines using netsearch, downhill simplex,
  and non-linear least-square Marquadt's method, designed for
  photoelectron diffraction calculation.
Constructors and Destructors:
  Pdchifit(int ndata,int nfit) --- ndata is the number of data of the
    curve to be fitted. nfit is the number of parameters to be fitted
    or held.
Member Functions:
  float getmemory() --- get memory it took in unit of bytes
  float loadcurve(int fitmath,int trymax,float tolerance,
    float *xdata,float *ydata,float *afit,float *safit,
    float (*fitfunc)(int fitmath,int ndata,int nfit,int *yafit,
    float *afit,float *safit,float *xdata,float *ydata,float *ymod,
    float *dyda,Mscdrun *mscdrun)) --- load curve to be fitted and
      initial input parameters.

    fitmath --- control of fitting routines. fitmath=1 using Marquadt's
      method only. fitmath=2 using Downhill simplex and then Marquadt's
      method. fitmath=3 using Netsearch first, then Downhill simplex,
      then Marquadt's method. fitmath=4 using Netsearch only.
    trymax --- the maximum number of trials
    tolerance --- convergence tolerance to be achieved.
    xdata and ydata --- arguments and values to be fitted.
    afit --- fitting parameters and their initial values
    safit --- control the fitting for each parameter, zero entry is
      given for the component that should be held fixed at its input
      value, positive entry is given for the component the trial step
      size in a netsearch or downhill simplex method, negative integer
      (-n) entry is given for the component that tie its value to the
      nth component such that the ratio of their values is kept fixed
      as determined by their input values.
    fitfunc(int fitmath,int ndata,int nfit,int *yafit,
      float *afit,float *safit,float *xdata,float *ydata,float *ymod,
      float *dyda,Mscdrun *mscdrun) --- the function to be fitted.

  float dofit(float *ybest,float *abest,float *reliability,
      int *lastrynum) --- enforce the fitting procedure.
    ybest --- the best fitted values corresponding to each argument
    abest --- the best fitted parameters
    reliability --- reliability factors between the best-fit value
      and their original input data
    lastrynum --- number of trials took by the fitting procedure.
----------------------------------------------------------------------
*/

#ifndef __PDCHIFIT_H
#define __PDCHIFIT_H

#include "mscdrun.h"

class Pdchifit
{ private:
    int ndata,nfit,mfit,fitmath,trynum,trymax,error;
    float tolerance,alamda,ochisq,chisq,reliable;
    int *yafit;
    float *xdata,*ydata,*afit,*safit,*invsig,*atry,*dafit,*beta,*covar,
      *alpha,*ymod,*dyda;
    Mscdrun *mscdrun;
  private:
    float (*fitfunc) (int fitmath,int ndata,int nfit,int *yafit,
      float *afit,float *safit,float *xdata,float *ydata,float *ymod,
      float *dyda,Mscdrun *mscdrun);
  public:
    Pdchifit();
    ~Pdchifit();
    void init(int ndata=0,int nfit=0);
    float getmemory();
    int loadcurve(int fitmath,int trymax,float tolerance,
      float *xdata,float *ydata,float *afit,float *safit,
      Mscdrun *imscdrun,
      float (*fitfunc)(int fitmath,int ndata,int nfit,int *yafit,
        float *afit,float *safit,float *xdata,float *ydata,float *ymod,
        float *dyda,Mscdrun *mscdrun));
    int dofit(float *ybest,float *abest,float *reliability,
      int *lastrynum);
  private:
    int pdangsearch();
    int netsearch();
    int downsimplex();
    int levenmarq();
    float downhill(int high,float factor);
    int marqmin();
    int marqcof();
    int gausjord(float *amat,float *bmat,int na,int nb);
  // xdata[ndata],ydata[ndata],invsig[ndata],afit[nfit],safit[nfit],
  // yafit[nfit],dafit[nfit],atry[nfit],beta[nfit],covar[nfit+1][nfit],
  // alpha[nfit][nfit],ymod[ndata],dyda[nfit][ndata]
};

#endif //__PDCHIFIT_H

