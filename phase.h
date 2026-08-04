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

1. Phaseshift
Description: conversion of scattering phase shift data
Constructors and Destructor:
  Phaseshift() --- constructor.
Member Functions:
  int getlnum() --- get lnum
  int getnumpoint() --- get number of data points
  float getmemory() --- get memory it took in unit of bytes
  int getalnum(float wavevec) --- return number of non-zero phase shift
    data
  int loadcurve(char *filename,char *atomname,char *symbolname,
    char *subshell)
  float fphase(float wavevec, int lquantum) --- return the scattering
    phase shift at wave vector = wavevec for quantum number lquantum
  Fcomplex fsinexp(float wavevec,int lquantum) --- return the
    scattering phase shift t matrix element sin(phase)*exp(i*phase)
  int savecurve(char *filename,char *usermessage) --- save the
    phase shift data versus wave vector into a file
    with filename
  int kconfine(float *kmin,float *kmax,int *alnum) --- confine kmin and
    kmax within the wave vector range of the phase shift data, and
    have alnum equals to lnum if alnum is greater
  int makephase(float wavevec) --- make phase shift data for all
    quantum numbers at a specific wave vector wavevec, and return the
    number of non-zero phase shift data
  Fcomplex fsinexpa(float wavevec,int lquantum) --- return the
    scattering phase shift t matrix element sin(phase)*exp(i*phase)
    without check
----------------------------------------------------------------------
*/

#ifndef __PHASE_H
#define __PHASE_H

#include "fcomplex.h"

class Phaseshift
{ private:
    int datatype,lnum,palnum,ndata,error;
    float pwavea,pwaveb;
    char *atomname,*symbolname;
    float *phasea,*phaseb;
    Fcomplex *phasec;
  public:
    Phaseshift();
    ~Phaseshift();
    void init();
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int getlnum();
    int getnumpoint();
    float getmemory();
    int getalnum(float wavevec);
    int loadcurve(char *filename,char *atomname,char *symbolname,
      int lnum=0);
    float fphase(float wavevec,int lquantum);
    Fcomplex fsinexp(float wavevec,int lquantum);
    int savecurve(char *filename,char *usermessage);
    int kconfine(float *kmin,float *kmax,int *alnum);
    Fcomplex fsinexpa(float wavevec,int lquantum);
private:
    int makephase(float wavevec);
    int makesinexp(float wavevec);
};

#endif //__PHASE_H

