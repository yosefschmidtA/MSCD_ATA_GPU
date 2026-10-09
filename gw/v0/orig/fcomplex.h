/* fcomplex.h
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

    Float complex Number Library - Include File
    class Fcomplex:  declarations for float complex numbers.

*/

#ifndef __FCOMPLEX_H
#define __FCOMPLEX_H

#include <math.h>
#include <iostream>

class Fcomplex
{
  // Implementation
  private:
    float re, im;

  public:
    // constructors
    Fcomplex(float __re_val, float __im_val=0.0);
    Fcomplex();

  // Fcomplex manipulations
    friend float real(Fcomplex);   // the real part
    friend float imag(Fcomplex);   // the imaginary part
    friend Fcomplex conj(Fcomplex);  // the Fcomplex conjugate
    friend float norm(Fcomplex);   // the square of the magnitude
    friend float arg(Fcomplex);    // the angle in the plane

    // Create a Fcomplex object given polar coordinates
    friend Fcomplex polar(float __mag, float __angle=0);

    // Overloaded ANSI C math functions
    friend float cabs(Fcomplex);
    friend Fcomplex acos(Fcomplex);
    friend Fcomplex asin(Fcomplex);
    friend Fcomplex atan(Fcomplex);
    friend Fcomplex cos(Fcomplex);
    friend Fcomplex cosh(Fcomplex);
    friend Fcomplex exp(Fcomplex);
    friend Fcomplex log(Fcomplex);
    friend Fcomplex log10(Fcomplex);
    friend Fcomplex pow(Fcomplex __base, float __expon);
    friend Fcomplex pow(float __base, Fcomplex __expon);
    friend Fcomplex pow(Fcomplex __base, Fcomplex __expon);
    friend Fcomplex sin(Fcomplex);
    friend Fcomplex sinh(Fcomplex);
    friend Fcomplex sqrt(Fcomplex);
    friend Fcomplex tan(Fcomplex);
    friend Fcomplex tanh(Fcomplex);

    // Binary Operator Functions
    friend Fcomplex operator+(Fcomplex, Fcomplex);
    friend Fcomplex operator+(float, Fcomplex);
    friend Fcomplex operator+(Fcomplex, float);
    friend Fcomplex operator-(Fcomplex, Fcomplex);
    friend Fcomplex operator-(float, Fcomplex);
    friend Fcomplex operator-(Fcomplex, float);
    friend Fcomplex operator*(Fcomplex, Fcomplex);
    friend Fcomplex operator*(Fcomplex, float);
    friend Fcomplex operator*(float, Fcomplex);
    friend Fcomplex operator/(Fcomplex, Fcomplex);
    friend Fcomplex operator/(Fcomplex, float);
    friend Fcomplex operator/(float, Fcomplex);
    friend int operator==(Fcomplex, Fcomplex);
    friend int operator!=(Fcomplex, Fcomplex);
    Fcomplex operator+=(Fcomplex);
    Fcomplex operator+=(float);
    Fcomplex operator-=(Fcomplex);
    Fcomplex operator-=(float);
    Fcomplex operator*=(Fcomplex);
    Fcomplex operator*=(float);
    Fcomplex operator/=(Fcomplex);
    Fcomplex operator/=(float);
    Fcomplex operator+();
    Fcomplex operator-();

    friend std::ostream& operator<<(std::ostream& os, Fcomplex __z);
    friend std::istream& operator>>(std::istream& is, Fcomplex __z);

};

#endif // __FCOMPLEX_H

