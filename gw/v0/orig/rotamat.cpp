#include <iostream>
#include <math.h>

#include "userutil.h"
#include "rotamat.h"

Rotamat::Rotamat()
{ error=107; init();
} //end of Rotamat::Rotamat

Rotamat::~Rotamat()
{ if (rotmata) delete [] rotmata; if (rotmatb) delete [] rotmatb;
  if (rotmatc) delete [] rotmatc;
} //end of Rotamat::~Rotamat

void Rotamat::init(int alnum,int amaxmag)
{ int k;
  if ((error==107)||(error==103))
  { rotmata=rotmatb=rotmatc=NULL;
  }
  if (error==103)
  { if (alnum<0) lnum=0;
    else if (alnum>60) lnum=60;
    else lnum=alnum;
    if (error==0) k=sizeof(int); else k=sizeof(int);
    if ((k<4)&&(lnum>20)) lnum=20;
    if (k<4) betanum=61; else betanum=361;
    if (amaxmag<0) maxmag=magnum=lamdum=0;
    else if ((k<4)&&(amaxmag>2))
    { maxmag=2; magnum=maxmag+maxmag+1; lamdum=(maxmag+1)*(maxmag+1);
    }
    else if (amaxmag>4)
    { maxmag=4; magnum=maxmag+maxmag+1; lamdum=(maxmag+1)*(maxmag+1);
    }
    else
    { maxmag=amaxmag; magnum=maxmag+maxmag+1;
      lamdum=(maxmag+1)*(maxmag+1);
    }
    rotmata=new float [betanum*lnum*lamdum];
    rotmatb=new float [lnum*lamdum];
    rotmatc=new float [lnum*lamdum];
    if ((!rotmata)||(!rotmatb)||(!rotmatc)) error=102;
  }
  else if (error==106)
  { if (rotmata)
    { rotmata=new float [betanum*lnum*lamdum];
      if (!rotmata) error=102;
    }
    if (rotmatb)
    { rotmatb=new float [lnum*lamdum];
      if (!rotmatb) error=102;
    }
    if (rotmatc)
    { rotmatc=new float [lnum*lamdum];
      if (!rotmatc) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Rotamat::init

int Rotamat::getlength()
{ int ka;
  if (error==0) ka=sizeof(Rotamat)+(betanum+2)*lnum*lamdum*
    sizeof(float)+sizeof(int);
  else ka=0;
  return(ka);
} //end of Rotamat::getlength

int Rotamat::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Rotamat);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=betanum*lnum*lamdum*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)rotmata,ka,kb,length);
    kb=lnum*lamdum*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)rotmatb,ka,kb,length);
    kb=lnum*lamdum*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)rotmatc,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Rotamat::paexport

int Rotamat::paimport(char *source,int length)
{ int ka,kb;
  
  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Rotamat);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=betanum*lnum*lamdum*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)rotmata,source,ka,kb,length);
    kb=lnum*lamdum*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)rotmatb,source,ka,kb,length);
    kb=lnum*lamdum*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)rotmatc,source,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Rotamat::paimport

/*
----------------------------------------------------------------------
index of rotmata and rotmatb, the rotation matrix element
k     0     1     2     3     4     5     6     7     8
ma    0     1     1     1     2     2     2     2     2
mb    0    -1     0     1    -2    -1     0     1     2

k     9    10    11    12    13    14    15
ma    3     3     3     3     3     3     3
mb   -3    -2    -1     0     1     2     3

k    16    17    18    19    20    21    22    23    24
ma    4     4     4     4     4     4     4     4     4
mb   -4    -3    -2    -1     0     1     2     3     4


ma=  -4    -3    -2    -1     0     1     2     3     4
mb=
-4                                                   24
-3                                             15    23
-2                                        8    14    22
-1                                  3     7    13    21
 0                  k =       0     2     6    12    20
 1                                  1     5    11    19
 2                                        4    10    18
 3                                              9    17
 4                                                   16

index of rotmatc, the spherical harmonic normalization factor
k     0     1     2     3     4     5     6     7     8
ma    0     0     1     1     0     1     2     2     2
mb    0     1     0     1     2     2     0     1     2

k     9    10    11    12    13    14    15
ma    0     1     2     3     3     3     3
mb    3     3     3     0     1     2     3

k    16    17    18    19    20    21    22    23    24
ma    0     1     2     3     4     4     4     4     4
mb    4     4     4     4     0     1     2     3     4


ma=   0     1     2     3     4
mb=
0     0     2     6    12    20
1     1     3     7    13    21
2     4     5     8    14    22
3     9    10    11    15    23
4    16    17    18    19    24

----------------------------------------------------------------------
*/
int Rotamat::makecurve()
{ int i,j,k,al,ma,mb,ta,tb,p;
  int itemp[8];
  double xa,xb,xc,rotsum;
  double rtemp[8];
  double *rotfac,*rotcoef;
  const float radian=(float)(3.14159265/180.0);

  rotfac=rotcoef=NULL;
  if ((error==101)||(error==0))
  { rotfac=new double [lnum+maxmag+1];
    rotcoef=new double [lnum];
    if ((!rotfac)||(!rotcoef)) error=102;
    else error=0;
  }
  if (error==0)
  { for (i=0;i<lnum+maxmag+1;++i)
    { if (i<2) rotfac[i]=1.0;
      else rotfac[i]=i*rotfac[i-1];
    }
  }

  for (al=0;(error==0)&&(al<lnum);++al)
  { for (k=0;k<lamdum;++k)
    { if (k<1) ma=mb=0;
      else if (k<4)
      { ma=1; mb=k-2;
      }
      else if (k<9)
      { ma=2; mb=k-6;
      }
      else if (k<16)
      { ma=3; mb=k-12;
      }
      else
      { ma=4; mb=k-20;
      }
      for (j=0;j<=al;++j)
      { itemp[0]=al+mb; itemp[1]=al-mb; itemp[2]=al+ma; itemp[3]=al-ma;
        itemp[4]=al+mb-j; itemp[5]=al-ma-j; itemp[6]=j-mb+ma;
        itemp[7]=j;
        for (p=0;p<8;++p)
        { if (itemp[p]<0) p=10;
          else if (p<4) rtemp[p]=sqrt(rotfac[itemp[p]]);
          else rtemp[p]=rotfac[itemp[p]];
        }
        if (p<10)
        { xb=rtemp[0]/rtemp[4]*rtemp[1]/rtemp[5];
          xc=rtemp[2]/rtemp[6]*rtemp[3]/rtemp[7];
          if ((j&1)!=0) xb=-xb;
          rotcoef[j]=xb*xc;
        }
        else rotcoef[j]=0.0;
      }
      for (i=0;i<betanum;++i)
      { xa=(double)(i*180.0*radian/(betanum-1)); xb=cos(xa/2.0);
        xc=sin(xa/2.0);
        rotsum=0.0;
        for (j=0;j<=al;++j)
        { xa=rotcoef[j];
          ta=al+al+mb-ma-j-j; tb=j+j-mb+ma;
          if ((fabs(xa)>1.0e-10)&&((ta==0)||(fabs(xb)>1.0e-10))&&
            ((tb==0)||(fabs(xc)>1.0e-10)))
          { if (ta!=0) xa*=pow(xb,ta);
            if (tb!=0) xa*=pow(xc,tb);
            rotsum+=xa;
          }
        }
        rotmata[i*lnum*lamdum+al*lamdum+k]=(float)rotsum;
      }

      if (k<2)
      { ma=0; mb=k;
      }
      else if (k<4)
      { ma=1; mb=k-2;
      }
      else if (k<6)
      { ma=k-4; mb=2;
      }
      else if (k<9)
      { ma=2; mb=k-6;
      }
      else if (k<12)
      { ma=k-9; mb=3;
      }
      else if (k<16)
      { ma=3; mb=k-12;
      }
      else if (k<20)
      { ma=k-16; mb=4;
      }
      else
      { ma=4; mb=k-20;
      }
      if ((al-ma>=0)&&(al-mb>=0)&&(al+ma>=0)&&(al+mb>=0))
      { xb=rotfac[al+ma]/rotfac[al+mb]*rotfac[al-mb]/rotfac[al-ma];
        xb=(double)(al+al+1.0)*sqrt(xb);
        if ((mb&1)!=0) xb=-xb;
      }
      else xb=0.0;
      rotmatc[al*lamdum+k]=(float)xb;
    }
  }

  pbeta=0.0f;
  for (al=0;(error==0)&&(al<lnum);++al)
  { for (k=0;k<lamdum;++k) rotmatb[al*lamdum+k]=rotmata[al*lamdum+k];
  }

  if (rotfac) delete [] rotfac; if (rotcoef) delete [] rotcoef;

  return(error);
} //end of Rotamat::makecurve

float Rotamat::getmemory()
{ float xa;
  if (error==0) xa=(float)(sizeof(Rotamat)+
    (betanum+2)*lnum*lamdum*sizeof(float));
  else xa=-1.0f;
  return(xa);
} //end of Rotamat::getmemory

//beta unit: degree
int Rotamat::makerotation(float beta)
{ int i,k,al;
  float abeta,xa;

  xa=(float)fabs(beta-pbeta); abeta=(float)fabs(beta);
  if ((error==0)&&(xa>0.1)&&(abeta<181.0))
  { pbeta=beta;
    xa=(float)(abeta*(betanum-1.0)/180.0); i=(int)xa;
    if (i<0) i=0;
    else if (i>betanum-2) i=betanum-2;
    for (al=0;al<lnum;++al)
    { for (k=0;k<lamdum;++k)
      { rotmatb[al*lamdum+k]=rotmata[i*lnum*lamdum+al*lamdum+k]+
          (xa-i)*(rotmata[(i+1)*lnum*lamdum+al*lamdum+k]-
          rotmata[i*lnum*lamdum+al*lamdum+k]);
        if ((beta<0.0)&&((k==2)||(k==5)||(k==7)||(k==10)||(k==12)||
          (k==14)||(k==17)||(k==19)||(k==21)||(k==23)))
          rotmatb[al*lamdum+k]=-rotmatb[al*lamdum+k];
      }
    }
  }
  return(error);
} //end of Rotamat::makerotation

//beta unit: degree
float Rotamat::rotelem(int al,int ma,int mb,float beta)
{ int k,ta,tb,sign;
  float element;

  if (error!=0) element=0.0f;
  else if ((al>=lnum)||(iabs(ma)>maxmag)||(iabs(mb)>maxmag)||
    (iabs(ma)>al)||(iabs(mb)>al))
    element=0.0f;
  else
  { if (fabs(beta-pbeta)>0.1) error=makerotation(beta);
    if (ma>=iabs(mb))
    { ta=ma; tb=mb; sign=1;
    }
    else if ((ma<0)&&(-ma>=iabs(mb)))
    { ta=-ma; tb=-mb;
      if ((ma&1)==(mb&1)) sign=1; else sign=-1;
    }
    else if (mb>0)
    { ta=mb; tb=ma;
      if ((ma&1)==(mb&1)) sign=1; else sign=-1;
    }
    else
    { ta=-mb; tb=-ma; sign=1;
    }
    k=ta*(ta+1)+tb;
    if (sign>=0) element=rotmatb[al*lamdum+k];
    else element=-rotmatb[al*lamdum+k];
  }
  return (element);
} //end of Rotamat::rotelem

//beta unit: degree
float Rotamat::rotharm(int al,int ma,int mb,float beta)
{ int k,ta,tb;
  float element;

  if (error!=0) element=0.0f;
  else if ((al>=lnum)||(iabs(ma)>maxmag)||(iabs(mb)>maxmag)||
    (iabs(ma)>al)||(iabs(mb)>al))
    element=0.0f;
  else
  { element=rotelem(al,ma,mb,beta);
    ta=iabs(ma); tb=iabs(mb);
    if ((ta==0)&&(tb<2)) k=tb;
    else if ((ta==1)&&(tb<2)) k=tb+2;
    else if ((ta<2)&&(tb==2)) k=ta+4;
    else if ((ta==2)&&(tb<3)) k=tb+6;
    else if ((ta<3)&&(tb==3)) k=ta+9;
    else if ((ta==3)&&(tb<4)) k=tb+12;
    else if ((ta<4)&&(tb==4)) k=ta+16;
    else k=tb+20;
    element*=rotmatc[al*lamdum+k];
  }
  return (element);
} //end of Rotamat::rotharm

//beta unit: degree
float Rotamat::termination(int al,int ma,int mb,float beta)
{ int k,ta,tb;
  float element;

  if (error!=0) element=0.0f;
  else if ((al>=lnum)||(iabs(ma)>maxmag)||(iabs(mb)>maxmag)||
    (iabs(ma)>al)||(iabs(mb)>al))
    element=0.0f;
  else
  { element=rotelem(al,ma,mb,beta);
    ta=0; tb=iabs(mb); if (error==0) ta=0;
    if ((ta==0)&&(tb<2)) k=tb;
    else if ((ta==1)&&(tb<2)) k=tb+2;
    else if ((ta<2)&&(tb==2)) k=ta+4;
    else if ((ta==2)&&(tb<3)) k=tb+6;
    else if ((ta<3)&&(tb==3)) k=ta+9;
    else if ((ta==3)&&(tb<4)) k=tb+12;
    else if ((ta<4)&&(tb==4)) k=ta+16;
    else k=tb+20;
    element*=rotmatc[al*lamdum+k]/(float)sqrt(al+al+1.0);
  }
  return (element);
} //end of Rotamat::termination

int Rotamat::getkelem(int ma,int mb)
{ int k,ta,tb,sign;
  if (ma>=iabs(mb))
  { ta=ma; tb=mb; sign=1;
  }
  else if ((ma<0)&&(-ma>=iabs(mb)))
  { ta=-ma; tb=-mb;
    if ((ma&1)==(mb&1)) sign=1; else sign=-1;
  }
  else if (mb>0)
  { ta=mb; tb=ma;
    if ((ma&1)==(mb&1)) sign=1; else sign=-1;
  }
  else
  { ta=-mb; tb=-ma; sign=1;
  }
  if (sign>0) k=ta*(ta+1)+tb;
  else k=-ta*(ta+1)-tb;
  return(k);
} //end of Rotamat::getkelem

int Rotamat::getkharm(int ma,int mb)
{ int k,ta,tb;
  ta=iabs(ma); tb=iabs(mb);
  if ((ta==0)&&(tb<2)) k=tb;
  else if ((ta==1)&&(tb<2)) k=tb+2;
  else if ((ta<2)&&(tb==2)) k=ta+4;
  else if ((ta==2)&&(tb<3)) k=tb+6;
  else if ((ta<3)&&(tb==3)) k=ta+9;
  else if ((ta==3)&&(tb<4)) k=tb+12;
  else if ((ta<4)&&(tb==4)) k=ta+16;
  else k=tb+20;
  return(k);
} //end of Rotamat::getkharm

//beta unit: degree
float Rotamat::rotharma(int al,int kelem,int kharm,float beta)
{ float element;

  if (al<lnum)
  { if (beta!=pbeta) error=makerotation(beta);
    if (kelem>=0) element=rotmatb[al*lamdum+kelem];
    else element=-rotmatb[al*lamdum-kelem];
    element*=rotmatc[al*lamdum+kharm];
  }
  else element=0.0f;

  return (element);
} //end of Rotamat::rotharma

