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

1. Vibration
Description: calculate correlated thermal vibrational mean square
  relative displacement versus interatomic distances
Constructors and Destructor:
  Vibration --- no arguments
Member Functions:
  int getnumpoint() --- get number of data points
  float getmemory() --- get memory it took in unit of bytes
  int loadparameter(float density,float mweight,float tdebye,
    float tsample) --- load parameters
    density --- density of the material in unit of g/cm3
    mweight --- molecular weight of the material in atomic unit
    tdebye --- debye temperature of the materials
    tsample --- environmental temperature at the sample
  float fvibmsrd(float length,float aweight) --- return the mean
    square displacement for the input bonding length and atomic weight
    aweight --- atomic weight of the specified atom in atomic unit
  int savecurve(char *filename,char *usermessage,float aweight)
    --- save the displacement versus interatomic distance into a file
    with filename
----------------------------------------------------------------------
*/

#ifndef __VIBRATE_H
#define __VIBRATE_H

class Vibration
{ private:
    int thernum,npoint,datatype,error;
    float density,mweight,tdebye,tsample,therstep,thermax;
    float *thermat;
  public:
    Vibration();
    ~Vibration();
    void init();
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int getnumpoint();
    float getmemory();
    int loadparameter(float density,float mweight,float tdebye,
      float tsample);
    float fvibmsrd(float length,float aweight);
    int savecurve(char *filename,char *usermessage,float aweight);
  private:
    int makecurve();
};

#endif //__VIBRATE_H

