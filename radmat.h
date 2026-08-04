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

1. Radialmatrix
Description: conversion of radial matrix element
Constructors and Destructor:
  Radialmatrix() --- constructor.
Member Functions:
  int getnumpoint() --- get number of data points
  int getnumfinal() --- get number of final states
  float getmemory() --- get memory it took in unit of bytes
  int loadcurve(char *filename,char *atomname,char *symbolname,
    char *subshell)
  float famphase(float wavevec, int code) --- return the radial matrix
    element amplitude for li+1 channel if code=1, phase for li+1 channel
    if code=2, amplitude for li-1 channel if code=3, and phase for li-1
    channel if code=4 at wavevec
  int savecurve(char *filename,char *usermessage) --- save the
    radial matrix element data versus wave vector into a file
    with filename
  int kconfine(float *kmin,float *kmax) --- confine kmin and kmax within
    the wave vector range of the matrix data
----------------------------------------------------------------------
*/

#ifndef __RADMAT_H
#define __RADMAT_H

class Radialmatrix
{ private:
    int datatype,ndata,error;
    char *atomname,*symbolname,*subshell;
    float *radmat;
  public:
    Radialmatrix();
    ~Radialmatrix();
    void init();
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int getnumpoint();
    int getnumfinal();
    float getmemory();
    int loadcurve(char *filename,char *atomname,char *symbolname,
      char *subshell);
    float famphase(float wavevec,int code);
    int savecurve(char *filename,char *usermessage);
    int kconfine(float *kmin,float *kmax);
};

#endif //__RADMAT_H

