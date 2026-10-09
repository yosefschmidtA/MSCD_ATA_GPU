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

1. Rotamat
Description: rotation matrix, Albert Messiah, Quantum Mechanics, vol II,
  John Wiley & Sons, Inc. Appendix IV. eqn (C.72)
Constructors and Destructor:
  Rotmat(int lnum,int maxmag) --- constructor, where lnum is
    the number of angular momenta, maxmag is the maximum value of
    magnetic momentum
Member Functions:
  int init() --- initialization, return unzero if error occured
  float getmemory() --- get memory it took in unit of bytes
  float rotelem(int al,int ma,int mb,float beta) --- return rotation
    matrix element at angle beta in degree for specified indecies
    (al,ma,mb)
  float rotharm(int al,int ma,int mb,float beta) --- return rotation
    matrix element at angle beta in degree for specified indecies
    (al,ma,mb), times a combined spherical harmonic normalization factor
    (2l+1)*Nlmb/Nlma*(-1)mb, J.J. Rehr and R.C. Albers, Phys. Rev. B 41,
    8139 (1990) eqn. (10) and (12)
  float termination(int al,int ma,int mb,float beta) --- return rotation
    matrix element at angle beta in degree for specified indecies
    (al,ma,mb), times a spherical harmonic normalization factor
    Nlm*(-1)m
  int getkelem(int ma,int mb) --- get address of rotation matrix
    element
  int getkharm(int ma,int mb) --- get address of spherical harmonic
    coefficient
  float rotharma(int al,int kelem,int kharm,float beta) --- return
    rotation matrix element at angle beta in degree for specified
    indecies of kelem and kharm
  int makerotation(float beta) --- calculate all the matrix elements
    at the specified angle beta in degree.
----------------------------------------------------------------------
*/

#ifndef __ROTAMAT_H
#define __ROTAMAT_H

class Rotamat
{ private:
    int error,lnum,maxmag,magnum,lamdum,betanum;
    float pbeta;
    float *rotmata,*rotmatb,*rotmatc;
  public:
    Rotamat();
    ~Rotamat();
    void init(int alnum=0,int amaxmag=0);
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int makecurve();
    float getmemory();
    float rotelem(int al,int ma,int mb,float beta);
    float rotharm(int al,int ma,int mb,float beta);
    float termination(int al,int ma,int mb,float beta);
    int getkelem(int ma,int mb);
    int getkharm(int ma,int mb);
    float rotharma(int al,int kelem,int kharm,float beta);
  private:
    int makerotation(float beta);
};

#endif //__ROTAMAT_H

/*
----------------------------------------------------------------------
rotmata[61*lnum*lamdum] rotmatb[lnum*lamdum] rotmatc[lnum*lamdum]

----------------------------------------------------------------------
*/

