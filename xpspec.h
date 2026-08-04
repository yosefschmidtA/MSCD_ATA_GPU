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
1. Spectrum
Description: convert xps spectrum file from other format to mscd format
Constructors and Destructors:
  Spectrum() --- constructor
Member Functions:
  int getnumpoint() --- get number of data points
  int getnumcurve() --- get number of curves
  float getmemory() --- get memory it took in unit of bytes
  int loadcurve(char *filename) --- load curve to be converted from
    a file with name filename
  int loadcomment(char *commentbuf,int bufsize) --- load comment
  int clean(int replace=0) --- clean the spectrum curves, sort the
    curves in order with photon energy. If replace set to 1, then
    remove curves with photon energy same as the next curve.
  int savecurve(char *filename,char *usermessage) --- save xps
    spectrum into file with name filename
  int savepeak(char *filename,char *usermessage) --- save xps
    spectrum peak intensities and their kinetic energies
  int savelist(char *filename,char *usermessage) --- save xps
    spectrum peak fit curve list
  int peakfit(int npeak,float *epeak,int backassociate=0,
    int peakassociate=0)
    --- fit Gaussian peaks of peak number npeak and peak positions
      in epeak.
      The Shirley background will be fitted and then fixed
      if backassociate is set to non-zero.
      The peak width and peak separation will be fitted
      and then fixed if associate is set to non-zero.


----------------------------------------------------------------------
*/

#ifndef __XPSPEC_H
#define __XPSPEC_H

class Spectrum
{ private:
    int datatype,ncurve,npoint,mpoint,npeak,nfit,numint,commentsize,
      xpsmode,error;
    char *comment;
    float *ebind,*xpscount,*curvepar,*emission,*bestfit;
  public:
    Spectrum();
    ~Spectrum();
    void init();
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int getnumpoint();
    int getnumcurve();
    float getmemory();
    int loadcurve(char *filename);
    int clean(int replace=0);
    int loadcomment(char *commentbuf,int bufsize);
    int peakfit(int npeak,float *epeak,int backassociate=0,
      int peakassociate=0);
    int savecurve(char *filename,char *usermessage);
    int savepeak(char *filename,char *usermessage);
    int savelist(char *filename,char *usermessage);
};

float fitxpspeak(int fitmath,int ndata,int nfit,int *yafit,
  float *afit,float *safit,float *xdata,float *ydata,float *ymod,
  float *dyda);

#endif //__XPSPEC_H
/*
----------------------------------------------------------------------
dimmensions
  ebind[mpoint] xpscount[mpoint] curvepar[ncurve*8]
  comment[commentsize] emission[ncurve*numint]
  bestfit[ncurve*nfit]


curvepar  0        1       2      3        4
          ephoton  theta   phi    memadd   memend+1

          5              6-7
          reliability    reserved

emission  0        1       2     3
          ephoton  theta   phi   reliability

          4                    5               6            7
          peak-binding-energy  peak-intensity  peak-height  peak-width

          8                    9               10           11
          peak-binding-energy  peak-intensity  peak-height  peak-width

          12                   13              14           15
          peak-binding-energy  peak-intensity  peak-height  peak-width

----------------------------------------------------------------------
*/

