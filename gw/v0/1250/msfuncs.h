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

1. Hankel
Description: calculate polynomial factor of spherical Hankel function
  times z^m/m! = clm(z), where z=-i/kx
Constructors and Destructor:
  Hankel(int ndata,int lnum,int cmnum) --- constructor, where ndata is
    the number of arguments, lnum is the number of quantum indecies,
    cmnum is the number of Rehr-Albers expansion indecies. Recommended
    number ndata=101, lnum=20, cmnum=4
Member Functions:
  int init() --- initialization, return unzero if error occured
  float getmemory() --- get memory it took in unit of bytes
  Fcomplex fhankelfac(int al,int am,float invkx) --- return clm(z) at
    z=-i*invkx for quantum number al and Rehr-Albers expansion number
    am.
  Fcomplex fhankelfaca(int al,int am,float invkx) --- return clm(z) at
    z=-i*invkx for quantum number al and Rehr-Albers expansion number
    am, without check.

2. Expix
Description: calculate exp(ix)
Constructors and Destructor:
  Expix(int ndata) --- constructor, where ndata is the number of
    internal arguments, 4001 is recommended
Member Functions:
  int init() --- initialization, return unzero if error occured
  float getmemory() --- get memory it took
  Fcomplex fexpix(float x) --- return exp(ix) for argument x in unit of
    degree
  Fcomplex fsinexp(float x) --- return sin(x)exp(ix) for argument x in
    unit of degree
----------------------------------------------------------------------
*/

#ifndef __MSFUNCS_H
#define __MSFUNCS_H

#include "fcomplex.h"

class Hankel
{ private:
    int ndata,cmnum,lnum,error;
    float argument;
    Fcomplex *hankmat,*hankarg;
  public:
    Hankel();
    ~Hankel();
    void init(int indata=0,int lnum=0,int cmnum=0);
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int makecurve();
    float getmemory();
    Fcomplex fhankelfac(int al,int am,float invkx);
    Fcomplex fhankelfaca(int al,int am,float invkx);
};

class Expix
{ private:
    int ndata,mdata,error;
    Fcomplex *cexpix, *csinexp;
  public:
    Expix();
    ~Expix();
    void init(int ndata=0);
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int makecurve();
    float getmemory();
    Fcomplex fexpix(float x);
    Fcomplex fsinexp(float x);
};

#endif //__MSFUNCS_H

