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

1. Pdintensity
Description: load intensity data from file, calculate chi or reliability
Constructors and Destructor:
  Pdintensity() --- constructor.
Member Functions:
  int getnumcurve() --- get number of curves
  int getnumpoint() --- get number of points
  float getmemory() --- get memory it took in unit of bytes
  int loadparameter(int idatatype,float kmin,float kmax,float kstep,
    float dtmin,float dtmax,float dtstep,float dpmin,float dpmax,
    float dpstep)
    --- calculate number of curves and data points, and allocate memory
  int loadpoint(int pnumber,float kout,float dtheta,
    float dphi,float ltheta,float lphi,float intensity,float chical,
    float chiexp) --- load one data point
  int setbeamangle(float ltheta,float lphi) --- set polarization angle
    for each data point
  int loadcomment(char *commentbuf,int commentsize) --- load comment
  int loadintensity(char *calfile,char *expfile,int loadmode,
    float *amin,float *amax,float *astep,float *bmin,float *cmin)
    --- load intensity data from file, amin amax and astep are minimum,
    maxmum and step of the scanning parameter, bmin and cmin are other
    two parameters
  int chicalc(int chicommand) --- calculate chi of the intensity
  int chicurve(int ndata,int nfit,float *xdata,float *ydata,
    float *ychi) --- calculate for one curve
  int normcurve(int ndata,float *ydata,float *ychi) --- calculate
    normalized chi for one curve
  int reliability(float *areliable,float *breliable,
      float *creliable=NULL, float *dreliable=NULL,
      float *ereliable=NULL, float *freliable=NULL,
      float *greliable=NULL, float *hreliable=NULL);
    --- calculate reliabilities between calculated and experimental
    chi data
  int savecurve(char *filename,char *usermessage) --- save intensity
    data
  int makemission(float adtheta,float adphi,float devstep)
    --- make a deviation rotation of emission direction
  int getfitdata(float *xdata,float *ydata) --- set xdata in ascending
    order of wave vector, ydata its experimental chi value
  int getcurvepar(int icurve,int *ndata,float *akout,
    float *adtheta,float *adphi,float *weightc,float *weightk)
    --- output parameters of the curve icurve
  int getpoint(int pnumber,float *kout,float *dtheta,
    float *dphi,float *ltheta,float *lphi,float *intensity,
    float *chical,float *chiexp) --- get parameters of pnumber point
----------------------------------------------------------------------
*/

#ifndef __PDINTEN_H
#define __PDINTEN_H

class Pdintensity
{ private:
    int datatype,npoint,nkout,ndtheta,ndphi,ndangle,ncurve,commentsize,
      error;
    float basemem;
    char *comment;
    float *pdata,*curvepar;
  public:
    Pdintensity();
    ~Pdintensity();
    void init();
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int getnumcurve();
    int getnumpoint();
    float getmemory();
    int loadparameter(int idatatype,float kmin,float kmax,float kstep,
      float dtmin,float dtmax,float dtstep,float dpmin,float dpmax,
      float dpstep);
    int loadpoint(int pnumber,float kout,float dtheta,
      float dphi,float ltheta,float lphi,float intena,float intenb,
      float chical,float chiexp);
    int setbeamangle(float ltheta,float lphi);
    int loadintensity(char *calfile,char *expfile,int loadmode,
      float *amin,float *amax,float *astep,float *bmin,float *cmin);
    int chicalc(int chicommand);
    int reliability(float *areliable,float *breliable,
      float *creliable=NULL, float *dreliable=NULL,
      float *ereliable=NULL, float *freliable=NULL,
      float *greliable=NULL, float *hreliable=NULL);
    int savecurve(char *filename,char *usermessage);
    int makemission(float adtheta,float adphi,
      float devstep);
    int getfitdata(float *xdata,float *ydata);
    int getcurvepar(int icurve,int *ndata,float *akout,
      float *adtheta,float *adphi,float *weightc,float *weightk);
    int getpoint(int pnumber,float *kout,float *dtheta,
      float *dphi,float *ltheta,float *lphi,float *intena,float *intenb,
      float *chical,float *chiexp);
  private:
    int loadcomment(char *commentbuf,int commentsize);
    int chicurve(int ndata,int nfit,float *xdata,float *ydata,
      float *ychi);
    int normcurve(int ndata,float *ydata,float *ychi);
};

#endif //__PDINTEN_H

/*
----------------------------------------------------------------------
pdata[npoint*12]  curvepar[ncurve*12]


pdata   0   1      2    3       4     5          6         7    8
        k   theta  phi  ptheta  pphi  intensity  reference chi  chiexp

        9-11
        reserved

curvepar     0     1     2     3     4     5       6       7      8
scanuni=1,5  kmin  kmax  kstep theta phi   weightc weightk numbeg numend
scanuni=2,6  tmin  tmax  tstep k     phi   weightc weightk numbeg numend
scanuni=3,7  pmin  pmax  pstep k     theta weightc weightk numbeg numend
scanuni=4,8  0.0   0.0   0.0   k     0.0   weightc weightk numbeg numend

             9-11
             reserved for later use

loadmode in loadintensity
     1   load calfile without chi
     2   load calfile
     3   load expfile
     4   load calfile and expfile with identical parameters
    11   load calfile without chi, using defined amin,amax,astep
    12   load calfile, using defined amin,amax,astep
    13   load expfile, using defined amin,amax,astep
    14   load calfile and expfile, using defined amin,amax,astep
 111-438 load expfile, with datatype=loadmode, and using defined amin,
           amax,astep
         (amin,amax,astep) is one of (kmin,kmax,kstep), (dtmin,dtmax,
           dtstep) and (dpmin,dpmax,dpstep)

chicommand
  11   calculate chi only for calculated intensity
  12   calculate normalization only for calculated chi data
  13   calculate chi and normalization for calculated intensity
  22   calculate normalization only for experimental chi data
  23   same as 22
  31   same as 11
  32   calculate normalization only for both calculated and experimental
         chi data
  33   calculate chi and normalization for calculated intensity, and
         normalization for experimental chi data
----------------------------------------------------------------------
*/
