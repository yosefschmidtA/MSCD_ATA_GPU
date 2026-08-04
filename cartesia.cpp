#include <iostream>
#include <math.h>

#include "userutil.h"
#include "cartesia.h"

Cartesia::Cartesia(float x_val,float y_val,float z_val)
{ x=x_val; y=y_val; z=z_val;
} //end of Cartesia::Cartesia

Cartesia Cartesia::operator+()
{ return *this;
} //end of Cartesia::operator+

Cartesia Cartesia::operator-()
{ return Cartesia(-x,-y,-z);
} //end of Cartesia::operator-

Cartesia& Cartesia::operator+=(Cartesia& __p2)
{ x+=__p2.x; y+=__p2.y; z+=__p2.z;
  return *this;
} //end of Cartesia::operator+=

Cartesia& Cartesia::operator-=(Cartesia& __p2)
{ x-=__p2.x; y-=__p2.y; z-=__p2.z;
  return *this;
} //end of Cartesia::operator-=

Cartesia operator+(Cartesia& __p1, Cartesia& __p2)
{ return Cartesia(__p1.x+__p2.x,__p1.y+__p2.y,__p1.z+__p2.z);
} //end of operator+

Cartesia operator-(Cartesia& __p1, Cartesia& __p2)
{ return Cartesia(__p1.x-__p2.x,__p1.y-__p2.y,__p1.z-__p2.z);
} //end of operator-

int operator==(Cartesia& __p1, Cartesia& __p2)
{ return __p1.x == __p2.x && __p1.y == __p2.y && __p1.z == __p2.z;
} //end of operator==

int operator!=(Cartesia& __p1, Cartesia& __p2)
{ return ((__p1.x!=__p2.x)||(__p1.y!=__p2.y)||(__p1.z!=__p2.z));
} //end of operator!=

void Cartesia::loadcoordinates(float x_val,float y_val,float z_val)
{ x=x_val; y=y_val; z=z_val;
  return;
} //end of Cartesia::loadcoordinates

void Cartesia::loadthetaphi(float atheta,float aphi)
{ const float radian=(float)(3.14159265/180.0);
  float temp;
  temp=(float)sin(atheta*radian); z=(float)cos(atheta*radian);
  x=temp*(float)cos(aphi*radian); y=temp*(float)sin(aphi*radian);
  return;
} //end of Cartesia::loadthetaphi

void Cartesia::getcoordinates(float *x_val,float *y_val,float *z_val)
{ *x_val=x; *y_val=y; *z_val=z;
  return;
} //end of Cartesia::getcoordinates

float Cartesia::length()
{ return ((float)sqrt(x*x+y*y+z*z));
} //end of Cartesia::length

float Cartesia::xylength()
{ return ((float)sqrt(x*x+y*y));
} //end of Cartesia::xylength

float Cartesia::theta()
{ float xa,ya;
  const float radian=(float)(3.14159265/180.0);
  xa=(float)sqrt(x*x+y*y);
  if ((xa==0.0)&&(z==0.0)) ya=0.0f;
  else ya=(float)atan2(sqrt(x*x+y*y),z)/radian;
  return(ya);
} //end of Cartesia::theta

float Cartesia::phi()
{ float ya;
  const float radian=(float)(3.14159265/180.0);
  if ((x==0.0)&&(y==0.0)) ya=0.0f;
  else ya=(float)atan2(y,x)/radian;
  return(ya);
} //end of Cartesia::phi

float Cartesia::costheta(Cartesia& __p2)
{ float xa,xb,xc,ya;
  Cartesia __p1(x,y,z);
  xa=length();
  xb=__p2.length();
  if ((fabs(xa)<1.0e-30)&&(fabs(xb)>1.0e-30))
    ya=__p2.xylength()/xb;
  else if ((fabs(xa)>1.0e-30)&&(fabs(xb)<1.0e-30))
    ya=xylength()/xa;
  else if ((fabs(xa)>1.0e-30)&&(fabs(xb)>1.0e-30))
  { __p1-=__p2; xc=__p1.length();
    ya=(xa*xa+xb*xb-xc*xc)/(xa*xb*2.0f);
  }
  else ya=1.0f;
  return(ya);
} //end of Cartesia::costheta

/*
----------------------------------------------------------------------
Cartesia Cartesia::euler(float alpha,float beta,float gamma)
  --- return the new vector after an Euler rotation (alpha,beta,gamma).
    The rotation is the result of the following sequence of three
      rotations:
    (i) a rotation of the vector of angle gamma about Oz axis, the frame
        or the cartesian coordinate system is fixed.
    (ii) a rotation of the vector of angle beta about Oy axis, the frame
        or the cartesian coordinate system is fixed.
    (iii) a rotation of the vector of angle alpha about Oz axis, the
        frame or the cartesian coordinate system is fixed.
----------------------------------------------------------------------
*/
Cartesia Cartesia::euler(float alpha,float beta,float gamma)
{ float xa,ya,za,sina,sinb,sing,cosa,cosb,cosg;
  const float radian=(float)(3.14159265/180.0);

  sina=(float)sin(alpha*radian); sinb=(float)sin(beta*radian);
  sing=(float)sin(gamma*radian);
  cosa=(float)cos(alpha*radian); cosb=(float)cos(beta*radian);
  cosg=(float)cos(gamma*radian);

  xa=(cosa*cosb*cosg-sina*sing)*x-(cosa*cosb*sing+sina*cosg)*y+
     cosa*sinb*z;
  ya=(sina*cosb*cosg+cosa*sing)*x-(sina*cosb*sing-cosa*cosg)*y+
     sina*sinb*z;
  za=-sinb*cosg*x+sinb*sing*y+cosb*z;
  return Cartesia(xa,ya,za);
} //end of Cartesia::euler

Cartesia& Cartesia::operator*=(Cartesia& __p2)
{ float tempx,tempy,tempz;
  tempx=y*__p2.z-z*__p2.y;
  tempy=z*__p2.x-x*__p2.z;
  tempz=x*__p2.y-y*__p2.x;
  x=tempx; y=tempy; z=tempz;
  return *this;
} //end of Cartesia::operator*=

Cartesia operator*(Cartesia& __p1, Cartesia& __p2)
{ float tempx,tempy,tempz;
  tempx=__p1.y*__p2.z-__p1.z*__p2.y;
  tempy=__p1.z*__p2.x-__p1.x*__p2.z;
  tempz=__p1.x*__p2.y-__p1.y*__p2.x;
  return Cartesia(tempx,tempy,tempz);
} //end of operator*

float operator/(Cartesia& __p1, Cartesia& __p2)
{ return (__p1.x*__p2.x+__p1.y*__p2.y+__p1.z*__p2.z);
} //end of operator/

ostream& operator<<(ostream& os, Cartesia& __p)
{ return os << '(' << __p.x << ", " << __p.y << ", " << __p.z << ')';
} //end of operator<<

istream& operator>>(istream& is, Cartesia& __p)
{ skipdelim(is); is>>__p.x; skipdelim(is); is>>__p.y;
  skipdelim(is); is>>__p.z;
  return is;
} //end of operator>>

