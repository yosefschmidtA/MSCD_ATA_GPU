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

1. Meanpath
Description: calculation of inelastic mean free path
Constructors and Destructor:
  Meanpath() --- constructor.
Member Functions:
  int getnumpoint() --- get number of points
  float getmemory() --- get memory it took in unit of bytes
  int loadparameter(float valence,float bandgap,
    float density,float mweight)
    valence --- number of valence electrons (>=1) or exponent of energy
      (0.5-0.9)
    bandgap --- bandgap energy (in eV) or coefficient of energy
      (0.02-0.3)
    density --- density of the bulk (in g/cm3)
    mweight --- molecular or atomic weight (in amu)
  float finvpath(float wavevector) --- return the inverse mean free path
    for specified wave vector
  int savecurve(char *filename,char *usermessage) --- save the
    mean free path versus wave vector into a file with filename
----------------------------------------------------------------------
*/

#ifndef __MEANPATH_H
#define __MEANPATH_H

class Meanpath
{ private:
    int formula,datatype,npoint,error;
    float valence,bandgap,density,mweight,kenergy,invpath;
    float *imfpc;
  public:
    Meanpath();
    ~Meanpath();
    void init();
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int getnumpoint();
    float getmemory();
    int loadparameter(float ivalence,float ibandgap,
      float idensity,float imweight);
    float finvpath(float wavevector);
    int savecurve(char *filename,char *usermessage);
};

#endif //__MEANPATH_H
