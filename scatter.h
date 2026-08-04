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

1. Scatter
Description: execute a job request
Constructors and Destructor:
  Scatter()  --- constructor
Member Functions:
  void init(int lnum,int raorder) --- initialization, lnum is the
    number of quantum momenta, raorder is Rehr-Albers approximation
    order. raorder=2 corresponding to 6x6 scattering factor matrix,
    raorder=1 to 3x3 matrix, raorder=0 to 1x1
  int loadcurve(char *filename,char *atomname,char *symbolname)
    --- load curve, filename is the name of phase shift file for
      atom atomname and symbol symbolname
  int getlength() --- get length of the object
  int paexport(char *dest,int length=0) --- paexport object to buffer
  int paimport(char *source,int length=0) --- paimport an object from buffer
  int getnumpoint() --- get number of points
  float getmemory() --- get memory of the object
  Fcomplex sevenelem(int ma,int na,int mb,int nb,
    float akin,float vka,float vkb,float beta) --- calculate one
      scattering factor element
  int makefactor(float akin,float lenga,float lengb) --- make scattering
    factor matrix for wave vector akin, and bonding lengths lenga and
    lengb, akin uses unit of inverse angstrom, lenga and lengb angstroms
  int savecurve(char *filename,char *usermessage) --- save scattering
    factor matrix into file filename



----------------------------------------------------------------------
*/

#ifndef __SCATTER_H
#define __SCATTER_H

#include "userutil.h"
#include "fcomplex.h"
#include "phase.h"
#include "rotamat.h"
#include "msfuncs.h"

class Scatter
{ private:
    int datatype,lnum,raorder,betanum,eledim,error;
    float wavevec,alength,blength;
    char *atomname,*symbolname;
    float *xbeta,*uevenelem;
    Fcomplex *tevenelem;
    Phaseshift *phaseshift;
    Rotamat *evenmat;
    Hankel *hanka,*hankb;
  public:
    Scatter();
    ~Scatter();
    void init(int lnum=0,int raorder=0);
    int loadcurve(char *filename,char *atomname,char *symbolname);
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int getnumpoint();
    float getmemory();
    Fcomplex sevenelem(int ma,int na,int mb,int nb,
      float akin,float vka,float vkb,float beta);
    int makefactor(float akin,float lenga,float lengb);
    int savecurve(char *filename,char *usermessage);
};

#endif //__SCATTER_H

/*
----------------------------------------------------------------------
dimmensions
  tevenelem[ndata*(eledim*eledim+1)]
    --- matrix element of triple vertex event

----------------------------------------------------------------------
*/
