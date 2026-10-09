#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

#include "polation.h"
#include "curvefit.h"
#include "cartesia.h"
#include "userutil.h"
#include "pdinten.h"

Pdintensity::Pdintensity()
{ error=107; init();
} //end of Pdintensity::Pdintensity

Pdintensity::~Pdintensity()
{ if (pdata) delete [] pdata; if (comment) delete [] comment;
  if (curvepar) delete [] curvepar;
} //end of Pdintensity::~Pdintensity

void Pdintensity::init()
{ if ((error==107)||(error==103))
  { datatype=111; npoint=ncurve=commentsize=0; basemem=0.0f;
    pdata=NULL; comment=NULL; curvepar=NULL;
  }
  if (error==106)
  { if (comment)
    { comment=new char [commentsize];
      if (!comment) error=102;
    }
    if (pdata)
    { pdata=new float [npoint*12];
      if (!pdata) error=102;
    }
    if (curvepar)
    { curvepar=new float [ncurve*12];
      if (!curvepar) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Pdintensity::init

int Pdintensity::getlength()
{ int ka;
  if (error==0)
  { ka=sizeof(int)+sizeof(Pdintensity);
    if (comment) ka+=commentsize*sizeof(char);
    if (pdata) ka+=npoint*12*sizeof(float);
    if (curvepar) ka+=ncurve*12*sizeof(float);
  }
  else ka=0;
  return(ka);
} //end of Pdintensity::getlength

int Pdintensity::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Pdintensity);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=commentsize;
    if (ka>=0) ka=memorysend(dest,(char *)comment,ka,kb,length);
    kb=npoint*12*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)pdata,ka,kb,length);
    kb=ncurve*12*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)curvepar,ka,kb,length);
    if ((ka<0)||(ka!=length)) error=901;
  }
  return(error);
} //end of Pdintensity::paexport

int Pdintensity::paimport(char *source,int length)
{ int ka,kb;

  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Pdintensity);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=commentsize;
    if (ka>=0) ka=memoryrec((char *)comment,source,ka,kb,length);
    kb=npoint*12*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)pdata,source,ka,kb,length);
    kb=ncurve*12*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)curvepar,source,ka,kb,length);
    if ((ka<0)||(ka!=length)) error=901;
  }
  return(error);
} //end of Pdintensity::paimport

int Pdintensity::loadparameter(int idatatype,float kmin,float kmax,
  float kstep,float dtmin,float dtmax,float dtstep,float dpmin,
  float dpmax,float dpstep)
{ int i,j,k,ma,mb,scanmode,scanhun,scanten,scanuni;
  float akout,atheta,aphi,ptstep;
  const float radian=(float)(3.14159265/180.0);

  if ((idatatype<=100)||(idatatype>=500)) error=207+(datatype<<16);
  datatype=scanmode=idatatype;
  if ((error==101)||(error==0))
  { error=0; scanhun=(scanmode%1000)/100; scanten=(scanmode%100)/10;
    scanuni=scanmode%10;

    if (kstep>0.0) nkout=(int)((kmax-kmin)/kstep+1.001); else nkout=1;
    nkout=confine(nkout,1,256);
    kmax=kmin+kstep*(nkout-1);
    if (dtstep>0.0) ndtheta=(int)((dtmax-dtmin)/dtstep+1.001);
    else ndtheta=1;
    ndtheta=confine(ndtheta,1,256);
    dtmax=(float)(dtmin+dtstep*(ndtheta-1));
    if ((dpstep>0.0)&&(dpmax>dpmin)&&(dtmax>0.1))
      ndphi=(int)((dpmax-dpmin)/dpstep+1.001);
    else ndphi=1;
    ndphi=confine(ndphi,1,256);
    dpmax=dpmin+dpstep*(ndphi-1);

    if ((scanuni==2)||(scanuni==3)||(scanuni==6)||(scanuni==7)||
      (ndtheta==1)||(ndphi==1)) ndangle=ndtheta*ndphi;
    else
    { ndangle=0;
      for (j=0;j<ndtheta;++j)
      { atheta=dtmin+j*dtstep;
        if (atheta<0.1) ++ndangle;
        else
        { ptstep=dpstep/(float)sin(atheta*radian);
          k=(int)((dpmax-dpmin)/ptstep+1.5);
          ndangle+=k;
        }
      }
    }
    npoint=nkout*ndangle;
    if ((scanuni==1)||(scanuni==5)) ncurve=ndangle;
    else if ((scanuni==2)||(scanuni==6)) ncurve=nkout*ndphi;
    else if ((scanuni==3)||(scanuni==7)) ncurve=nkout*ndtheta;
    else if ((scanuni==4)||(scanuni==8)) ncurve=nkout;
    else ncurve=0;
    if (ncurve<=0) error=703+(datatype<<16);
    else if (ncurve>100) error=704+(datatype<<16);
    else if (npoint<=0) error=701+(datatype<<16);
    else if (npoint>30000) error=702+(datatype<<16);

    if ((scanuni==0)&&(ndtheta>1)&&(ndphi>1)) scanuni=4;
    else if ((scanuni==0)&&(nkout>=10)&&(nkout>=ndtheta)&&
      (nkout>=ndphi)) scanuni=1;
    else if ((scanuni==0)&&(ndtheta>=10)&&(ndtheta>=nkout)&&
      (ndtheta>=ndphi)) scanuni=2;
    else if ((scanuni==0)&&(ndphi>=10)&&(ndphi>=nkout)&&
      (ndphi>=ndtheta)) scanuni=3;
    else if ((scanuni==4)&&(ndtheta==1)&&(ndphi==1)) scanuni=1;
    else if ((scanuni==4)&&(ndphi==1)) scanuni=2;
    else if ((scanuni==4)&&(ndtheta==1)) scanuni=3;
    else if ((scanuni==8)&&(ndtheta==1)&&(ndphi==1)) scanuni=5;
    else if ((scanuni==8)&&(ndphi==1)) scanuni=6;
    else if ((scanuni==8)&&(ndtheta==1)) scanuni=7;
    else if ((scanuni<1)||(scanuni>8)) scanuni=1;

    if ((scanhun>1)&&((scanuni==1)||(scanuni==5))&&(nkout<10))
      scanhun=1;
    else if ((scanhun>1)&&((scanuni==2)||(scanuni==6))&&(ndtheta<10))
      scanhun=1;
    else if ((scanhun>1)&&((scanuni==3)||(scanuni==7))&&(ndphi<10))
      scanhun=1;
    else if ((scanhun>1)&&((scanuni==4)||(scanuni==8))&&(ndtheta<10)&&
      (ndphi<8)) scanhun=1;
    else if ((scanhun!=2)&&(scanhun!=3)) scanhun=1;
    if ((scanten<1)||(scanten>3)) scanten=1;
    datatype=scanmode=scanhun*100+scanten*10+scanuni;

    if (pdata) delete [] pdata;
    if (curvepar) delete [] curvepar;
    pdata=new float [npoint*12]; curvepar=new float [ncurve*12];
    if ((!pdata)||(!curvepar)) error=102;
    else
    { error=0;
      for (i=0;i<npoint*12;++i) pdata[i]=0.0f;
      for (i=0;i<ncurve;++i)
      { if ((scanuni==1)||(scanuni==5))
        { if (ndtheta<=1)
          { atheta=dtmin; aphi=dpmin+i*dpstep;
          }
          else if (ndphi<=1)
          { atheta=dtmin+i*dtstep; aphi=dpmin;
          }
          else
          { ma=0; atheta=aphi=0.0f;
            for (j=0;j<ndtheta;++j)
            { atheta=dtmin+j*dtstep;
              if (atheta<0.1)
              { aphi=dpmin; ++ma;
              }
              else
              { ptstep=dpstep/(float)sin(atheta*radian);
                k=(int)((dpmax-dpmin)/ptstep+1.5);
                if (k<2) ptstep=dpmax-dpmin;
                else ptstep=(dpmax-dpmin)/(k-1);
                aphi=dpmin+(i-ma)*ptstep;
                ma+=k;
              }
              if (ma>i) break;
            }
          }
          if (atheta<0.1) aphi=0.0f;
          curvepar[i*12]=kmin; curvepar[i*12+1]=kmax;
          curvepar[i*12+2]=kstep;
          curvepar[i*12+3]=atheta; curvepar[i*12+4]=aphi;
          curvepar[i*12+5]=1.0f; curvepar[i*12+6]=0.0f;
          k=i*nkout; curvepar[i*12+7]=(float)k;
          for (j=0;(j<nkout)&&(j+k<npoint);++j)
          { pdata[(j+k)*12]=kmin+j*kstep;
            pdata[(j+k)*12+1]=atheta; pdata[(j+k)*12+2]=aphi;
          }
        }
        else if ((scanuni==2)||(scanuni==6))
        { j=i/ndphi; k=i%ndphi;
          curvepar[i*12]=dtmin; curvepar[i*12+1]=dtmax;
          curvepar[i*12+2]=dtstep;
          curvepar[i*12+3]=akout=kmin+j*kstep;
          curvepar[i*12+4]=aphi=dpmin+k*dpstep;
          curvepar[i*12+5]=1.0f; curvepar[i*12+6]=0.0f;
          k=i*ndtheta; curvepar[i*12+7]=(float)k;
          for (j=0;(j<ndtheta)&&(j+k<npoint);++j)
          { pdata[(j+k)*12+1]=dtmin+j*dtstep;
            pdata[(j+k)*12]=akout; pdata[(j+k)*12+2]=aphi;
          }
        }
        else if ((scanuni==3)||(scanuni==7))
        { j=i/ndtheta; k=i%ndtheta;
          curvepar[i*12]=dpmin; curvepar[i*12+1]=dpmax;
          curvepar[i*12+2]=dpstep;
          curvepar[i*12+3]=akout=kmin+j*kstep;
          curvepar[i*12+4]=atheta=dtmin+k*dtstep;
          curvepar[i*12+5]=1.0f; curvepar[i*12+6]=0.0f;
          k=i*ndphi; curvepar[i*12+7]=(float)k;
          for (j=0;(j<ndphi)&&(j+k<npoint);++j)
          { pdata[(j+k)*12+2]=dpmin+j*dpstep;
            pdata[(j+k)*12]=akout; pdata[(j+k)*12+1]=atheta;
          }
        }
        else if ((scanuni==4)||(scanuni==8))
        { curvepar[i*12]=0.0f; curvepar[i*12+1]=0.0f;
          curvepar[i*12+2]=0.0f;
          curvepar[i*12+3]=akout=kmin+(float)i*kstep;
          curvepar[i*12+4]=0.0f;
          curvepar[i*12+5]=1.0f; curvepar[i*12+6]=0.0f;
          k=i*ndangle; curvepar[i*12+7]=(float)k;
          for (j=0;j<ndtheta;++j)
          { atheta=(dtmin+j*dtstep);
            if ((atheta<0.1)&&(k<npoint))
            { pdata[k*12]=akout; pdata[k*12+1]=atheta;
              pdata[k*12+2]=0.0f; ++k;
            }
            else
            { ptstep=dpstep/(float)sin(atheta*radian);
              ma=(int)((dpmax-dpmin)/ptstep+1.5);
              if (ma>1) ptstep=(dpmax-dpmin)/(ma-1);
              for (mb=0;(mb<ma)&&(k<npoint);++mb)
              { pdata[k*12]=akout; pdata[k*12+1]=atheta;
                pdata[k*12+2]=dpmin+mb*ptstep; ++k;
              }
            }
          }
        }
      }
      for (i=0;i<npoint;++i)
      { pdata[i*12+10]=pdata[i*12+1]; pdata[i*12+11]=pdata[i*12+2];
        pdata[i*12+8]=0.0f;
      }
    }
  }

  return(error);
} //end of Pdintensity::loadparameter

/*
----------------------------------------------------------------------
loadmode in loadintensity
     1   load calfile without chi
     2   load calfile
     3   load expfile
     4   load calfile and expfile with identical parameters
    11   load calfile without chi, using defined amin,amax,astep
    12   load calfile, using defined amin,amax,astep
    13   load expfile, using defined amin,amax,astep
    14   load calfile and expfile, using defined amin,amax,astep
 111-438 load expfile, with datatype=loadmode, and using defined amin,
           amax,astep
         (amin,amax,astep) is one of (kmin,kmax,kstep), (dtmin,dtmax,
           dtstep) and (dpmin,dpmax,dpstep)
----------------------------------------------------------------------
*/
int Pdintensity::loadintensity(char *calfile,char *expfile,int loadmode,
  float *amin,float *amax,float *astep,float *bmin,float *cmin)
{ int i,j,k,m,ndata,scanmode,scanuni,scanten,scanhun,begrow,memsize,
    linenum;
  char *membuf,*commentbuf;
  float xa,xb,xc,xd,emin,emax,estep;
  float *xdata,*ydata,*zdata,*udata,*vdata;

  membuf=commentbuf=NULL; memsize=256;
  xdata=ydata=zdata=udata=vdata=NULL;
  if ((error==101)||(error==0))
  { error=0;
    membuf=new char [memsize]; commentbuf=new char [100*80];
    if ((!membuf)||(!commentbuf)) error=102;
  }

  scanmode=datatype=scanuni=scanten=scanhun=ndata=linenum=0;
  emin=emax=estep=0.0f;
  if ((error==0)&&((loadmode<1)||((loadmode>3)&&(loadmode<11))||
    ((loadmode>14)&&(loadmode<111))||(loadmode>=500))) error=901;
  else if ((error==0)&&(loadmode<10))
  { emin=0.0f; emax=1000.0f; estep=0.01f;
  }
  else if ((error==0)&&((*astep<1.0e-3)||(*amax<=*amin)))
  { emax=emin=*amin; estep=0.0f;
  }
  else if (error==0)
  { emin=*amin; emax=*amax; estep=*astep;
  }

  if ((error==0)&&((loadmode==1)||(loadmode==2)||(loadmode==4)||
    (loadmode==11)||(loadmode==12)||(loadmode==14)))
  { std::ifstream fin(calfile,ios::in);
    if (!fin) error=201+(datatype<<16);
    else
    { fin >> xa >> xb >> xc; skipendline(fin);
      datatype=(int)xa; begrow=(int)xb; linenum=(int)xc;
      if (linenum>30000) linenum=30000;
      if ((!fin)||(datatype<=100)||(datatype>=500))
        error=207+(datatype<<16);
      else if (begrow>1)
      { scanmode=datatype; scanuni=scanmode%10;
        scanten=(scanmode/10)%10; scanhun=(scanmode/100)%10;
        for (i=2;i<begrow;++i) skipendline(fin);
      }
      if ((loadmode<100)&&((scanuni==4)||(scanuni==8))) loadmode%=10;
    }
    if ((error==0)&&(linenum>0)) ncurve=1;
    else if (error==0)
    { fin >> xa; skipendline(fin); ncurve=(int)xa;
    }
    if (ncurve<1) ncurve=1;
    else if (ncurve>100) ncurve=100;

    k=0;
    for (i=0;(error==0)&&(i<ncurve);++i)
    { if ((linenum<=0)&&(scanuni!=4)&&(scanuni!=8))
      { fin >> xa >> xb >> xc >> xc >> xd; skipendline(fin);
        ndata=(int)xb;
      }
      else if (linenum<=0)
      { fin >> xa >> xb >> xc >> xd; skipendline(fin);
        ndata=(int)xb;
      }
      else
      { ndata=linenum; xd=1.0f;
      }
      if (xd<1.0e-3)
      { for (j=0;j<ndata;++j) skipendline(fin);
        --ncurve; --i;
        if (ncurve<1) error=701;
      }
      else if (error==0)
      { fin >> xa; xc=xa;
        if ((scanuni==4)||(scanuni==8)) fin >> xa;
        skipendline(fin);
        if (ndata>1)
        { fin >> xb; xd=xb;
          if ((scanuni==4)||(scanuni==8)) fin >> xb;
          skipendline(fin);
          if ((fin)&&(scanuni!=4)&&(scanuni!=8)&&(xb<=xa))
            error=207+(datatype<<16);
          else if ((fin)&&((xd<xc)||((xd==xc)&&(xb<=xa))))
            error=207+(datatype<<16);
        }
        else xb=xa;
        if (emin<xa) emin=xa;
        if ((scanuni!=4)&&(scanuni!=8)&&(estep<xb-xa)) estep=xb-xa;
        for (j=2;j<ndata;++j)
        { xa=xb; xc=xd; fin >> xb; xd=xb;
          if ((scanuni==4)||(scanuni==8)) fin >> xb;
          skipendline(fin);
          if ((fin)&&(scanuni!=4)&&(scanuni!=8)&&(xb<=xa))
            error=207+(datatype<<16);
          else if ((fin)&&((xd<xc)||((xd==xc)&&(xb<=xa))))
            error=207+(datatype<<16);
        }
        if (emax>xb) emax=xb;
        if (!fin) error=207+(datatype<<16);
        k+=ndata;
      }
    }
    npoint=k;
  }

  if ((error==0)&&((loadmode==3)||(loadmode==13)||(loadmode==4)||
    (loadmode==14)||(loadmode>100)))
  { std::ifstream fin(expfile,ios::in);
    if (!fin) error=201+(datatype<<16);
    else
    { fin >> xa >> xb >> xc; skipendline(fin);
      j=(int)xa; begrow=(int)xb; linenum=(int)xc;
      if ((loadmode==3)||(loadmode==13)) datatype=j;
      else if (loadmode>100) datatype=loadmode;
      k=((datatype%100)/10)-((j%100)/10); m=(datatype%10)-(j%10);
      if (linenum>30000) linenum=30000;
      if ((!fin)||(j<=100)||(j>=500))
        error=207+(datatype<<16);
      else if ((error==0)&&((k!=0)||((m!=0)&&(m!=4)&&(m!=-4))))
        error=613+(datatype<<16);
      else if ((error==0)&&(m==-4)) datatype+=4;
      if ((error==0)&&(begrow>1))
      { scanmode=datatype; scanuni=scanmode%10;
        scanten=(scanmode/10)%10; scanhun=(j/100)%10;
        for (i=2;i<begrow;++i) skipendline(fin);
      }
    }
    if ((error==0)&&(linenum>0)) j=1;
    else if (error==0)
    { fin >> xa; skipendline(fin); j=(int)xa;
    }
    else j=1;
    if (j<1) j=1;
    else if (j>100) j=100;
    if ((loadmode==3)||(loadmode==13)||(loadmode>100)||(j<ncurve))
      ncurve=j;

    k=0;
    for (i=0;(error==0)&&(i<ncurve);++i)
    { if ((linenum<=0)&&(scanuni!=4)&&(scanuni!=8))
      { fin >> xa >> xb >> xc >> xc >> xd; skipendline(fin);
        ndata=(int)xb;
      }
      else if (linenum<=0)
      { fin >> xa >> xb >> xc >> xd; skipendline(fin);
        ndata=(int)xb;
      }
      else
      { ndata=linenum; xd=1.0f;
      }
      if (xd<1.0e-3)
      { for (j=0;j<ndata;++j) skipendline(fin);
        --ncurve; --i;
        if (ncurve<1) error=701;
      }
      else if (error==0)
      { fin >> xa; xc=xa;
        if ((scanuni==4)||(scanuni==8)) fin >> xa;
        skipendline(fin);
        if (ndata>1)
        { fin >> xb; xd=xb;
          if ((scanuni==4)||(scanuni==8)) fin >> xb;
          skipendline(fin);
          if ((fin)&&(scanuni!=4)&&(scanuni!=8)&&(xb<=xa))
            error=207+(datatype<<16);
          else if ((fin)&&((xd<xc)||((xd==xc)&&(xb<=xa))))
            error=207+(datatype<<16);
        }
        else xb=xa;
        if (emin<xa) emin=xa;
        if ((loadmode<100)&&(scanuni!=4)&&(scanuni!=8)&&(estep<xb-xa))
          estep=xb-xa;
        for (j=2;j<ndata;++j)
        { xa=xb; xc=xd; fin >> xb; xd=xb;
          if ((scanuni==4)||(scanuni==8)) fin >> xb;
          skipendline(fin);
          if ((fin)&&(scanuni!=4)&&(scanuni!=8)&&(xb<=xa))
            error=207+(datatype<<16);
          else if ((fin)&&((xd<xc)||((xd==xc)&&(xb<=xa))))
            error=207+(datatype<<16);
        }
        if (emax>xb) emax=xb;
        if (!fin) error=207+(datatype<<16);
        k+=ndata;
      }
    }
    if ((loadmode==3)||(loadmode==13)) npoint=k;
    else if ((loadmode>100)&&((scanuni==4)||(scanuni==8))) npoint=k;
    else if ((loadmode==4)&&(npoint!=k)) error=613+(datatype<<16);
  }

  if ((error==0)&&(loadmode>10)&&((emax<emin)||(estep<1.0e-3)))
  { emax=emin; estep=0.0f; ndata=1; npoint=ncurve;
  }
  else if ((error==0)&&(loadmode>10))
  { ndata=(int)((emax-emin)/estep+1.001); npoint=ncurve*ndata;
  }

  if (error==0) k=sizeof(int); else k=0;
  if ((error==0)&&(npoint>30000)) error=617;
  else if ((error==0)&&(npoint>1000)&&(k<4)) error=617;
  else if ((error==0)&&(ncurve>=100)) error=704;
  else if (error==0)
  { if ((scanuni==4)||(scanuni==8)) k=4000; else k=256;
    if (curvepar) delete [] curvepar;
    if (pdata) delete [] pdata;
    xdata=new float [k]; ydata=new float [k];
    zdata=new float [k]; udata=new float [k];
    vdata=new float [k];
    curvepar=new float [ncurve*12];
    pdata=new float [npoint*12];
    if ((!xdata)||(!ydata)||(!zdata)||(!udata)||(!vdata)||
      (!curvepar)||(!pdata))
      error=102;
    else for (i=0;i<npoint*12;++i) pdata[i]=0.0f;
  }

  if ((error==0)&&((loadmode==1)||(loadmode==2)||(loadmode==4)||
    (loadmode==11)||(loadmode==12)||(loadmode==14)))
  { std::ifstream fin(calfile,ios::in);
    if (!fin) error=201+(datatype<<16);
    else
    { fin >> xa >> xb >> xc; skipendline(fin);
      j=(int)xa; begrow=(int)xb; linenum=(int)xc;
      if (linenum>30000) linenum=30000;
      if ((!fin)||(j<=100)||(j>=500))
        error=207+(datatype<<16);
      else if (begrow>1)
      { scanmode=datatype; scanuni=scanmode%10;
        scanten=(scanmode/10)%10; scanhun=(j/100)%10;
        k=0; commentbuf[0]='\0';
        for (i=2;i<begrow;++i)
        { getendline(fin,membuf,memsize);
          k+=stringlength(membuf)+3;
          if (k<100*80) stringappend(commentbuf,membuf,0,"\n");
        }
        commentbuf[100*80-1]='\0';
        if (!fin) error=207+(datatype<<16);
        else error=loadcomment(commentbuf,100*80);
      }
      else error=207+(datatype<<16);
    }
    if ((error==0)&&(linenum<=0))
    { fin >> xa >> xa; fin >> xa >> xb >> xc >> xd;
      skipendline(fin);
      nkout=(int)xa; ndtheta=(int)xb; ndphi=(int)xc; ndangle=(int)xd;
    }
    else nkout=ndtheta=ndphi=ndangle=1;
    if ((scanuni==1)||(scanuni==5))
    { nkout=ndata; ndangle=ncurve;
    }
    else if ((scanuni==2)||(scanuni==6))
    { ndtheta=ndata; ndangle=ndtheta*ndphi;
    }
    else if ((scanuni==3)||(scanuni==7))
    { ndphi=ndata; ndangle=ndtheta*ndphi;
    }
    for (i=0;(error==0)&&(i<ncurve);++i)
    { if ((linenum<=0)&&(scanuni!=4)&&(scanuni!=8))
      { fin >> xa >> xb >> curvepar[i*12+3] >> curvepar[i*12+4]
          >> curvepar[i*12+5] >> curvepar[i*12+6];
        skipendline(fin);
        m=(int)xb;
      }
      else if (linenum<=0)
      { fin >> xa >> xb >> curvepar[i*12+3]
          >> curvepar[i*12+5] >> curvepar[i*12+6];
        skipendline(fin);
        m=(int)xb; curvepar[i*12+4]=0.0f;
      }
      else
      { m=linenum;
        curvepar[i*12+3]=*bmin; curvepar[i*12+4]=*cmin;
        curvepar[i*12+5]=1.0f; curvepar[i*12+6]=0.0f;
      }
      if ((error==0)&&(ndata<10)) error=701;
      else if ((error==0)&&(scanuni!=4)&&(scanuni!=8)&&
        ((ndata>256)||(m>256))) error=702;
      else if ((error==0)&&((ndata>4000)||(m>4000))) error=702;
      else if ((error==0)&&(!fin)) error=207+(datatype<<16);
      else if ((error==0)&&(linenum<=0)&&(curvepar[i*12+5]<1.0e-3))
      { for (j=0;j<m;++j) skipendline(fin);
        --i;
      }
      else if (error==0)
      { for (j=0;(error==0)&&(j<m);++j)
        { fin >> xdata[j];
          if ((scanuni==4)||(scanuni==8)) fin >> udata[j];
          if (scanhun>3)
          { fin >> vdata[j]; ydata[j]=zdata[j]=0.0f;
          }
          else if (((loadmode==1)||(loadmode==11))&&
            (scanmode>200)&&(scanten<3))
          { fin >> ydata[j]; zdata[j]=vdata[j]=0.0f;
          }
          else if ((loadmode==1)||(loadmode==11))
          { fin >> ydata[j] >> zdata[j]; vdata[j]=0.0f;
          }
          else fin >> ydata[j] >> zdata[j] >> vdata[j];
          skipendline(fin);
          if (!fin) error=207+(datatype<<16);
        }
        if ((error==0)&&(loadmode<10))
        { if (i==0) k=0;
          curvepar[i*12+7]=(float)k;
          for (j=0;(j<m)&&(j+k<npoint);++j)
          { if ((scanuni==1)||(scanuni==5))
            { pdata[(j+k)*12]=xdata[j];
              pdata[(j+k)*12+1]=curvepar[i*12+3];
              pdata[(j+k)*12+2]=curvepar[i*12+4];
            }
            else if ((scanuni==2)||(scanuni==6))
            { pdata[(j+k)*12+1]=xdata[j];
              pdata[(j+k)*12]=curvepar[i*12+3];
              pdata[(j+k)*12+2]=curvepar[i*12+4];
            }
            else if ((scanuni==3)||(scanuni==7))
            { pdata[(j+k)*12+2]=xdata[j];
              pdata[(j+k)*12]=curvepar[i*12+3];
              pdata[(j+k)*12+1]=curvepar[i*12+4];
            }
            else if ((scanuni==4)||(scanuni==8))
            { pdata[(j+k)*12+1]=xdata[j]; pdata[(j+k)*12+2]=udata[j];
              pdata[(j+k)*12]=curvepar[i*12+3];
            }
            pdata[(j+k)*12+5]=ydata[j]; pdata[(j+k)*12+6]=zdata[j];
            pdata[(j+k)*12+7]=vdata[j];
          }
          k+=m;
        }
        else if (error==0)
        { Polation pdcurve; pdcurve.init(m);
          k=i*ndata; curvepar[i*12+7]=(float)k;
          for (j=0;(error==0)&&(j<ndata)&&(j+k<npoint);++j)
          { if ((scanuni==1)||(scanuni==5))
            { pdata[(j+k)*12]=xa=emin+j*estep;
              pdata[(j+k)*12+1]=curvepar[i*12+3];
              pdata[(j+k)*12+2]=curvepar[i*12+4];
            }
            else if ((scanuni==2)||(scanuni==6))
            { pdata[(j+k)*12+1]=xa=emin+j*estep;
              pdata[(j+k)*12]=curvepar[i*12+3];
              pdata[(j+k)*12+2]=curvepar[i*12+4];
            }
            else if ((scanuni==3)||(scanuni==7))
            { pdata[(j+k)*12+2]=xa=emin+j*estep;
              pdata[(j+k)*12]=curvepar[i*12+3];
              pdata[(j+k)*12+1]=curvepar[i*12+4];
            }
          }
          if (error==0) error=pdcurve.loadcurve(xdata,ydata);
          xa=pdcurve.getmemory();
          if ((error==0)&&(basemem<xa)) basemem=xa;
          for (j=0;(error==0)&&(j<ndata)&&(j+k<npoint);++j)
          { xa=emin+j*estep;
            pdata[(j+k)*12+5]=pdcurve.fspline(xa);
          }
          if (error==0) error=pdcurve.loadcurve(xdata,zdata);
          for (j=0;(error==0)&&(j<ndata)&&(j+k<npoint);++j)
          { xa=emin+j*estep;
            pdata[(j+k)*12+6]=pdcurve.fspline(xa);
          }
          if (error==0) error=pdcurve.loadcurve(xdata,vdata);
          for (j=0;(error==0)&&(j<ndata)&&(j+k<npoint);++j)
          { xa=emin+j*estep;
            pdata[(j+k)*12+7]=pdcurve.fspline(xa);
          }
        }
      }
    }
  }

  if ((error==0)&&((loadmode==3)||(loadmode==4)||(loadmode==13)||
    (loadmode==14)||(loadmode>100)))
  { std::ifstream fin(expfile,ios::in);
    if (!fin) error=201+(datatype<<16);
    else
    { fin >> xa >> xb >> xc; skipendline(fin);
      j=(int)xa; begrow=(int)xb; linenum=(int)xc;
      if (linenum>30000) linenum=30000;
      if ((!fin)||(j<=100)||(j>=500))
        error=207+(datatype<<16);
      else if (begrow>1)
      { scanmode=datatype; scanuni=scanmode%10;
        scanten=(scanmode/10)%10; scanhun=(j/100)%10;
        for (i=2;i<begrow;++i) skipendline(fin);
      }
      else error=207+(datatype<<16);
    }
    if ((error==0)&&(linenum<=0))
    { fin >> xa >> xa; fin >> xa >> xb >> xc >> xd;
      skipendline(fin);
      nkout=(int)xa; ndtheta=(int)xb; ndphi=(int)xc; ndangle=(int)xd;
    }
    else nkout=ndtheta=ndphi=ndangle=1;
    if ((scanuni==1)||(scanuni==5))
    { nkout=ndata; ndangle=ncurve;
    }
    else if ((scanuni==2)||(scanuni==6))
    { ndtheta=ndata; ndangle=ndtheta*ndphi;
    }
    else if ((scanuni==3)||(scanuni==7))
    { ndphi=ndata; ndangle=ndtheta*ndphi;
    }
    for (i=0;(error==0)&&(i<ncurve);++i)
    { if ((linenum<=0)&&(scanuni!=4)&&(scanuni!=8))
      { fin >> xa >> xb >> curvepar[i*12+3] >> curvepar[i*12+4]
          >> curvepar[i*12+5] >> curvepar[i*12+6];
        skipendline(fin);
        m=(int)xb;
      }
      else if (linenum<=0)
      { fin >> xa >> xb >> curvepar[i*12+3]
          >> curvepar[i*12+5] >> curvepar[i*12+6];
        skipendline(fin);
        m=(int)xb; curvepar[i*12+4]=0.0f;
      }
      else
      { m=linenum;
        curvepar[i*12+3]=*bmin; curvepar[i*12+4]=*cmin;
        curvepar[i*12+5]=1.0f; curvepar[i*12+6]=0.0f;
      }
      if ((error==0)&&(ndata<10)) error=701;
      else if ((error==0)&&(scanuni!=4)&&(scanuni!=8)&&
        ((ndata>256)||(m>256))) error=702;
      else if ((error==0)&&((ndata>4000)||(m>4000))) error=702;
      else if ((error==0)&&(!fin)) error=207+(datatype<<16);
      else if ((error==0)&&(linenum<=0)&&(curvepar[i*12+5]<1.0e-3))
      { for (j=0;j<m;++j) skipendline(fin);
        --i;
      }
      else if (error==0)
      { for (j=0;(error==0)&&(j<m);++j)
        { fin >> xdata[j];
          if ((scanuni==4)||(scanuni==8)) fin >> udata[j];
          if (scanhun>3)
          { fin >> vdata[j]; ydata[j]=zdata[j]=0.0f;
          }
          else if (((loadmode==1)||(loadmode==11))&&
            (scanmode>200)&&(scanten<3))
          { fin >> ydata[j]; zdata[j]=vdata[j]=0.0f;
          }
          else if ((loadmode==1)||(loadmode==11))
          { fin >> ydata[j] >> zdata[j]; vdata[j]=0.0f;
          }
          else fin >> ydata[j] >> zdata[j] >> vdata[j];
          skipendline(fin);
          if (!fin) error=207+(datatype<<16);
        }
        if ((error==0)&&(loadmode<10))
        { if (i==0) k=0;
          curvepar[i*12+7]=(float)k;
          for (j=0;(j<m)&&(j+k<npoint);++j)
          { if ((scanuni==1)||(scanuni==5))
            { pdata[(j+k)*12]=xdata[j];
              pdata[(j+k)*12+1]=curvepar[i*12+3];
              pdata[(j+k)*12+2]=curvepar[i*12+4];
            }
            else if ((scanuni==2)||(scanuni==6))
            { pdata[(j+k)*12+1]=xdata[j];
              pdata[(j+k)*12]=curvepar[i*12+3];
              pdata[(j+k)*12+2]=curvepar[i*12+4];
            }
            else if ((scanuni==3)||(scanuni==7))
            { pdata[(j+k)*12+2]=xdata[j];
              pdata[(j+k)*12]=curvepar[i*12+3];
              pdata[(j+k)*12+1]=curvepar[i*12+4];
            }
            else if ((scanuni==4)||(scanuni==8))
            { pdata[(j+k)*12+1]=xdata[j]; pdata[(j+k)*12+2]=udata[j];
              pdata[(j+k)*12]=curvepar[i*12+3];
            }
            if ((loadmode==3)||(loadmode==4))
              pdata[(j+k)*12+8]=vdata[j];
          }
          k+=m;
        }
        else if (error==0)
        { Polation pdcurve; pdcurve.init(m);
          error=pdcurve.loadcurve(xdata,vdata);
          xa=pdcurve.getmemory();
          if ((error==0)&&(basemem<xa)) basemem=xa;
          k=i*ndata; curvepar[i*12+7]=(float)k;
          for (j=0;(error==0)&&(j<ndata)&&(j+k<npoint);++j)
          { if ((scanuni==1)||(scanuni==5))
            { pdata[(j+k)*12]=xa=emin+j*estep;
              pdata[(j+k)*12+1]=curvepar[i*12+3];
              pdata[(j+k)*12+2]=curvepar[i*12+4];
            }
            else if ((scanuni==2)||(scanuni==6))
            { pdata[(j+k)*12+1]=xa=emin+j*estep;
              pdata[(j+k)*12]=curvepar[i*12+3];
              pdata[(j+k)*12+2]=curvepar[i*12+4];
            }
            else if ((scanuni==3)||(scanuni==7))
            { pdata[(j+k)*12+2]=xa=emin+j*estep;
              pdata[(j+k)*12]=curvepar[i*12+3];
              pdata[(j+k)*12+1]=curvepar[i*12+4];
            }
            pdata[(j+k)*12+8]=pdcurve.fspline(xa);
          }
        }
      }
    }
  }

  if (error==0)
  { for (i=0;i<npoint;++i)
    { pdata[i*12+10]=pdata[i*12+1]; pdata[i*12+11]=pdata[i*12+2];
      if ((loadmode==1)||(loadmode==2)||(loadmode==11)||
        (loadmode==12)) pdata[i*12+8]=0.0f;
    }
  }

  scanuni=datatype%10;
  if ((error==0)&&(scanuni>4)&&((loadmode==2)||(loadmode==12)))
    error=chicalc(12);
  else if ((error==0)&&(scanuni>4)&&((loadmode==3)||(loadmode==13)||
    (loadmode>100)))
    error=chicalc(22);
  else if ((error==0)&&(scanuni>4)&&((loadmode==4)||(loadmode==14)))
    error=chicalc(32);

  if (error==0)
  { *amin=emin; *amax=emax; *astep=estep;
  }

  if (membuf) delete [] membuf; if (commentbuf) delete [] commentbuf;
  if (xdata) delete [] xdata; if (ydata) delete [] ydata;
  if (zdata) delete [] zdata; if (udata) delete [] udata;
  if (vdata) delete [] vdata;
  return(error);
} //end of Pdintensity::loadintensity

