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

1. Cartesian coordinate system
Description: operations of the Cartesian coordinate system
Constructors and Destructor:
  Cartesia() --- constructor.
  Cartesia(float x, float y, float z) --- (x,y,z) is the initial
    coordinates
Member Functions:
  void loadcoordinates(float x_val,float y_val,float z_val) ---
    load coordinates for x, y, z
  void loadthetaphi(float atheta,float aphi) --- load theta and phi
    values in unit of degree, x=sin(atheta)*cos(aphi),
    y=sin(atheta)*sin(aphi), z=cos(atheta)
  void getcoordinates(float *x_val,float *y_val,float *z_val) ---
    get x y and z into x_val, y_val and z_val
  float length() --- return the length of the vector from (000) to (xyz)
  float xylength() --- return the length of the vector of projection in
    xy plane, from (000) to (xy0)
  float theta() --- return the polar angle of the vector in unit of
    degree
  float phi() --- return the azimuthal angle of the vector in unit of
    degree
  float costheta(Cartesia& __p2) --- return the cosine of the angle
    between the vector and another vector __p2.
  Cartesia euler(float alpha,float beta,float gamma) --- return the new
    vector after an Euler rotation (alpha, beta, gamma) in unit of
    degrees. The rotation is the result of the following sequence of
    three rotations:
    (i) a rotation of the vector of angle gamma about Oz axis, the frame
        or the cartesian coordinate system is fixed.
    (ii) a rotation of the vector of angle beta about Oy axis, the frame
        or the cartesian coordinate system is fixed.
    (iii) a rotation of the vector of angle alpha about Oz axis, the
        frame or the cartesian coordinate system is fixed.
    Note: This series of rotations is equivalent to the following
      rotations:
    (i) a rotation of angle alpha about Oz axis for the vector and the
        frame together, fram S goes into S', Oy axis goes into Oy'
    (ii) a rotation of angle beta about new Oy' axis for the vector and
        the frame together, frame S' goes into S'', Oz axis goes into
        Oz''
    (iii) a rotation of angle gamma about Oz'' axis for the vector and
        the frame together, frame S'' goes into S''', Oy' goes into
        Oy'''.

----------------------------------------------------------------------
*/

#ifndef __CARTESIA_H
#define __CARTESIA_H

#include <iostream>

class Cartesia
{ private:
    float x,y,z;
  public:
    Cartesia(float x_val=0.0,float y_val=0.0,float z_val=0.0);
    void loadcoordinates(float x_val=0.0,float y_val=0.0,
      float z_val=0.0);
    void loadthetaphi(float atheta=0.0,float aphi=0.0);
    void getcoordinates(float *x_val,float *y_val,float *z_val);
    float length();
    float xylength();
    float theta();
    float phi();
    float costheta(Cartesia& __p2);
    Cartesia euler(float alpha,float beta,float gamma);

    friend Cartesia operator+(Cartesia&, Cartesia&);
    friend Cartesia operator-(Cartesia&, Cartesia&);
    friend Cartesia operator*(Cartesia&, Cartesia&);
    friend float operator/(Cartesia&, Cartesia&); // dot multiplication
    friend int operator==(Cartesia&, Cartesia&);
    friend int operator!=(Cartesia&, Cartesia&);
    Cartesia& operator+=(Cartesia&);
    Cartesia& operator-=(Cartesia&);
    Cartesia& operator*=(Cartesia&);
    Cartesia operator+();
    Cartesia operator-();
    friend std::ostream& operator<<(std::ostream& os, Cartesia& __p);
    friend std::istream& operator>>(std::istream& is, Cartesia& __p);
};

#endif //__CARTESIA_H

