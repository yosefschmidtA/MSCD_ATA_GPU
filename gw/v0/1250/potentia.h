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

1. Potential
Description: calculation of electronic potential
Constructors and Destructor:
  Potential() --- constructor.
Member Functions:
  int getnumpoint() --- get number of points
  float getmemory() --- get memory it took in unit of bytes
  int loadcurve(char *filename,char *atomname,char *symbolname)
  float fpoten(float distance) --- return the potential times radial
    distance at specified radial distance
  int savecurve(char *filename,char *usermessage) --- save the
    potential times radial distance versus radial distance into a file
    with filename
----------------------------------------------------------------------
*/

#ifndef __POTENTIA_H
#define __POTENTIA_H

class Potential
{ private:
    int datatype,ndata,error;
    char *atomname,*symbolname;
    float *radist,*radpot;
  public:
    Potential();
    ~Potential();
    void init();
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int getnumpoint();
    float getmemory();
    int loadcurve(char *filename,char *atomname,char *symbolname);
    float fpoten(float distance);
    int savecurve(char *filename,char *usermessage);
};

#endif //__POTENTIA_H

