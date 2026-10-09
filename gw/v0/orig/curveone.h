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

1. Polation
Description: interpolation procedure
Constructors and Destructor:
  Polation(int ndata) --- ndata is the number of points of the curve.
Member Functions:
  float getmemory() --- get memory it took in unit of bytes
  int loadcurve(float *xdata,float *ydata) --- load curve to be
    interpolated.
  float fspline(float x) --- return the spline interpolated value of
    the input x.


List of functions

1.  fitspline --- spline interpolation function for fitting, used by
      curvefit object

----------------------------------------------------------------------
*/

#ifndef __CURVEONE_H
#define __CURVEONE_H

#ifndef NULL
#define NULL 0
#endif

class Curveone
{ private:
    int ndata,nfft,nsmooth,nderiv,error;
    float *xdata,*ydata,*ysmooth,*yfderiv,*ysderiv,*xfft,*yfft,*zfft;
  public:
    Curveone();
    ~Curveone();
    void init(int ndata=0);
    float getmemory();
    int loadcurve(float *xdata,float *ydata);
    float fsmooth(float xa);
    float ffderiv(float xa);
    float fsderiv(float xa);
    int smoothcurve(int insmooth=3);
    int derivative(int insmooth=3);
    int realfft(int backward=0);
    int getsdmini(float xa,float *xm,float *ym);
private:
    int dofft(int backward=0);
    int fftsmooth(int insmooth=3);
  // xdata[ndata],ydata[ndata],ysmooth[ndata],xfft[nfft],yfft[nfft];
};

#endif //__CURVEONE_H

