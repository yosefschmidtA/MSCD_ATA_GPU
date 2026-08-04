#include <iostream>
#include <math.h>

#include "fcomplex.h"
#include "userutil.h"

// Fcomplex functions

Fcomplex::Fcomplex(float __re_val, float __im_val)
{ re = __re_val; im = __im_val;
} //end of Fcomplex::Fcomplex

Fcomplex::Fcomplex()
{ re = im = 0.0f;
} //end of Fcomplex::Fcomplex

Fcomplex Fcomplex::operator+()
{ return *this;
} //end of Fcomplex::operator+

Fcomplex Fcomplex::operator-()
{ return Fcomplex(-re, -im);
} //end of Fcomplex::operator-

// Definitions of compound-assignment operator member functions

Fcomplex Fcomplex::operator+=(Fcomplex __z2)
{ re += __z2.re; im += __z2.im;
  return *this;
} //end of Fcomplex::operator+=

Fcomplex Fcomplex::operator+=(float __re_val2)
{ re += __re_val2;
  return *this;
} //end of Fcomplex::operator+=

Fcomplex Fcomplex::operator-=(Fcomplex __z2)
{ re -= __z2.re; im -= __z2.im;
  return *this;
} //end of Fcomplex::operator-=

Fcomplex Fcomplex::operator-=(float __re_val2)
{ re -= __re_val2;
  return *this;
} //end of Fcomplex::operator-=

Fcomplex Fcomplex::operator*=(float __re_val2)
{ re *= __re_val2; im *= __re_val2;
  return *this;
} //end of Fcomplex::operator*=

Fcomplex Fcomplex::operator/=(float __re_val2)
{ re /= __re_val2; im /= __re_val2;
  return *this;
} //end of Fcomplex::operator/=

Fcomplex Fcomplex::operator*=(Fcomplex __z2)
{ float temp;
  temp=re*__z2.re-im*__z2.im;
  im=re*__z2.im+im*__z2.re; re=temp;
  return *this;
} //end of Fcomplex::operator*=

Fcomplex Fcomplex::operator/=(Fcomplex __z2)
{ float tempa,tempb;
  tempa=(__z2.re*__z2.re+__z2.im*__z2.im);
  tempb=(re*__z2.re+im*__z2.im)/tempa;
  im=(im*__z2.re-re*__z2.im)/tempa; re=tempb;
  return *this;
} //end of Fcomplex::operator/=


// Definitions of non-member Fcomplex functions

float real(Fcomplex __z)
{ return __z.re;
} //end of real

float imag(Fcomplex __z)
{ return __z.im;
} //end of imag

Fcomplex conj(Fcomplex __z)
{ return Fcomplex(__z.re, -__z.im);
} //end of conj

// Definitions of non-member binary operator functions

Fcomplex operator+(Fcomplex __z1, Fcomplex __z2)
{ return Fcomplex(__z1.re + __z2.re, __z1.im + __z2.im);
} //end of operator+

Fcomplex operator+(float __re_val1, Fcomplex __z2)
{ return Fcomplex(__re_val1 + __z2.re, __z2.im);
} //end of operator+

Fcomplex operator+(Fcomplex __z1, float __re_val2)
{ return Fcomplex(__z1.re + __re_val2, __z1.im);
} //end of operator+

Fcomplex operator-(Fcomplex __z1, Fcomplex __z2)
{ return Fcomplex(__z1.re - __z2.re, __z1.im - __z2.im);
} //end of operator-

Fcomplex operator-(float __re_val1, Fcomplex __z2)
{ return Fcomplex(__re_val1 - __z2.re, -__z2.im);
} //end of operator-

Fcomplex operator-(Fcomplex __z1, float __re_val2)
{ return Fcomplex(__z1.re - __re_val2, __z1.im);
} //end of operator-

Fcomplex operator*(Fcomplex __z1, float __re_val2)
{ return Fcomplex(__z1.re*__re_val2, __z1.im*__re_val2);
} //end of operator*

Fcomplex operator*(float __re_val1, Fcomplex __z2)
{ return Fcomplex(__z2.re*__re_val1, __z2.im*__re_val1);
} //end of operator*

Fcomplex operator/(Fcomplex __z1, float __re_val2)
{ return Fcomplex(__z1.re/__re_val2, __z1.im/__re_val2);
} //end of operator/

int operator==(Fcomplex __z1, Fcomplex __z2)
{ return __z1.re == __z2.re && __z1.im == __z2.im;
} //end of operator==

int operator!=(Fcomplex __z1, Fcomplex __z2)
{ return __z1.re != __z2.re || __z1.im != __z2.im;
} //end of operator!=

Fcomplex operator*(Fcomplex __z1, Fcomplex __z2)
{ return Fcomplex(__z1.re*__z2.re-__z1.im*__z2.im,
    __z1.re*__z2.im+__z1.im*__z2.re);
} //end of Fcomplex operator*

Fcomplex operator/(Fcomplex __z1, Fcomplex __z2)
{ float tempa,tempb,tempc;
  tempa=(__z2.re*__z2.re+__z2.im*__z2.im);
  tempb=(__z1.re*__z2.re+__z1.im*__z2.im)/tempa;
  tempc=(__z1.im*__z2.re-__z1.re*__z2.im)/tempa;
  return Fcomplex(tempb, tempc);
} //end of Fcomplex operator/

Fcomplex operator/(float __re_val1, Fcomplex __z2)
{ float temp;
  temp=__re_val1/(__z2.re*__z2.re+__z2.im*__z2.im);
  return Fcomplex(__z2.re*temp, -__z2.im*temp);
} //end of Fcomplex operator/

float cabs(Fcomplex __z)
{ return((float)sqrt(__z.re*__z.re+__z.im*__z.im));
} //end of cabs

float norm(Fcomplex __z)
{ return (__z.re*__z.re+__z.im*__z.im);
} //end of norm

float arg(Fcomplex __z)
{ float temp;
  const float localeps=1.0e-30f;
  if ((__z.im>localeps)||(__z.im<-localeps)||(__z.re>localeps)||
    (__z.re<-localeps))
    temp=(float)atan2((double)__z.im,(double)__z.re);
  else temp=0.0f;
  return (temp);
} //end of arg

Fcomplex polar(float __mag, float __angle)
{ return Fcomplex(__mag*(float)cos((double)__angle),
    __mag*(float)sin((double)__angle));
} //end of polar

ostream& operator<<(ostream& os, Fcomplex __z)
{ return os << '(' << __z.re << ", " << __z.im << ')';
} //end of operator<<

istream& operator>>(istream& is, Fcomplex __z)
{ skipdelim(is); is>>__z.re; skipdelim(is); is>>__z.im;
  return is;
} //end of operator>>

