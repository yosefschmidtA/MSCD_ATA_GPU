#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

#include "polation.h"
#include "curvefit.h"
#include "cartesia.h"
#include "userutil.h"
#include "pdinten.h"

int Pdintensity::getnumcurve()
{ int k;
  if (error==0) k=ncurve; else k=0;
  return(k);
} //end of Pdintensity::getnumcurve

int Pdintensity::getnumpoint()
{ int k;
  if (error==0) k=npoint; else k=0;
  return(k);
} //end of Pdintensity::getnumpoint

float Pdintensity::getmemory()
{ float xa;
  if (error==0)
    xa=basemem+sizeof(Pdintensity)+(256+100*80)*sizeof(char)+
      (npoint*13+ncurve*12+4000*4+256*3+80)*sizeof(float);
  else xa=-1.0f;
  return(xa);
} //end of Pdintensity::getmemory

/*
----------------------------------------------------------------------
pdata   0   1      2    3       4     5       6      7       8
        k   theta  phi  ptheta  pphi  intena  intenb chical  chiexp

        9-11
        reserved

curvepar     0     1     2     3     4     5       6       7      8
scanuni=1,5  kmin  kmax  kstep theta phi   weightc weightk numbeg numend
scanuni=2,6  tmin  tmax  tstep k     phi   weightc weightk numbeg numend
scanuni=3,7  pmin  pmax  pstep k     theta weightc weightk numbeg numend
scanuni=4,8  0.0   0.0   0.0   k     0.0   weightc weightk numbeg numend
----------------------------------------------------------------------
*/
int Pdintensity::loadpoint(int pnumber,float kout,float dtheta,
  float dphi,float ltheta,float lphi,float intena,float intenb,
  float chical,float chiexp)
{ if ((error==0)&&(pnumber>=0)&&(pnumber<npoint))
  { pdata[pnumber*12]=kout; pdata[pnumber*12+1]=dtheta;
    pdata[pnumber*12+2]=dphi; pdata[pnumber*12+3]=ltheta;
    pdata[pnumber*12+4]=lphi;
    pdata[pnumber*12+5]=intena; pdata[pnumber*12+6]=intenb;
    pdata[pnumber*12+7]=chical; pdata[pnumber*12+8]=chiexp;
  }
  return(error);
} //end of Pdintensity::loadpoint

int Pdintensity::setbeamangle(float ltheta,float lphi)
{ int i,scanten;
  float atheta,aphi;
  const float radian=(float)(3.14159265/180.0);

  scanten=(datatype/10)%10;
  for (i=0;(error==0)&&(i<npoint);++i)
  { atheta=pdata[i*12+1]; aphi=pdata[i*12+2];
    if ((scanten==2)||(scanten==4)||(scanten==6)||(scanten==8))
    { Cartesia aold((float)sin(ltheta*radian)*(float)cos(lphi*radian),
        (float)sin(ltheta*radian)*(float)sin(lphi*radian),
        (float)cos(ltheta*radian));
      Cartesia anew;
      anew=aold.euler(aphi,atheta,0.0f);
      pdata[i*12+3]=anew.theta(); pdata[i*12+4]=anew.phi();
    }
    else
    { pdata[i*12+3]=ltheta; pdata[i*12+4]=lphi;
    }
  }

  return(error);
} //end of Pdintensity::setbeamangle

int Pdintensity::getpoint(int pnumber,float *kout,float *dtheta,
  float *dphi,float *ltheta,float *lphi,float *intena,
  float *intenb,float *chical,float *chiexp)
{ if (pnumber<0) pnumber=0; else if (pnumber>npoint-1) pnumber=npoint-1;
  *kout=pdata[pnumber*12]; *dtheta=pdata[pnumber*12+1];
  *dphi=pdata[pnumber*12+2]; *ltheta=pdata[pnumber*12+3];
  *lphi=pdata[pnumber*12+4];
  *intena=pdata[pnumber*12+5]; *intenb=pdata[pnumber*12+6];
  *chical=pdata[pnumber*12+7]; *chiexp=pdata[pnumber*12+8];
  return(error);
} //end of Pdintensity::getpoint

int Pdintensity::loadcomment(char *commentbuf,int bufsize)
{ int i;
  if (error==0)
  { commentsize=stringlength(commentbuf);
    if (commentsize>bufsize) commentsize=bufsize;
    if (comment) delete [] comment;
    comment=new char [commentsize];
    if (!comment) error=102;
    else
    { for (i=0;i<commentsize-1;++i) comment[i]=commentbuf[i];
      if ((comment[commentsize-2]=='\n')||
        (comment[commentsize-2]=='\r')) comment[commentsize-2]='\0';
      comment[commentsize-1]='\0';
    }
    if (commentsize<5)
    { if (comment) delete [] comment;
      commentsize=0; comment=NULL;
    }
  }

  return(error);
} //end of Pdintensity::loadcomment

/*
----------------------------------------------------------------------
chicommand
  11   calculate chi only for calculated intensity
  12   calculate normalization only for calculated chi data
  13   calculate chi and normalization for calculated intensity
  22   calculate normalization only for experimental chi data
  23   same as 22
  31   same as 11
  32   calculate normalization only for both calculated and experimental
         chi data
  33   calculate chi and normalization for calculated intensity, and
         normalization for experimental chi data
----------------------------------------------------------------------
*/
int Pdintensity::chicalc(int chicommand)
{ int i,j,k,m,p,ndata,nfit,scanhun,scanten,scanuni;
  float xa,xb,xc,sa,sb,atheta,patheta;
  float *xdata,*ydata,*zdata,*ychi;

  xdata=ydata=zdata=ychi=NULL;
  scanuni=datatype%10; scanten=(datatype/10)%10;
  scanhun=(datatype/100)%10;
  k=chicommand/10; m=chicommand%10;
  if ((scanuni<1)||(scanuni>8)) chicommand=0;
  else if ((chicommand<11)||(chicommand>33)) chicommand=0;
  else if ((scanhun<1)||(scanhun>3)) chicommand=0;
  else if ((m<1)&&(m>3)) chicommand=0;
  else if (chicommand==21) chicommand=0;
  else if (chicommand==31) chicommand=11;
  else if (chicommand==23) chicommand=22;
  if (((m==2)||(m==3))&&((datatype%10)<=4)) datatype+=4;
  if ((m==1)&&((datatype%10)>4)) datatype-=4;

  if (error==0)
  { xdata=new float [256]; ydata=new float [256];
    zdata=new float [256]; ychi=new float [256];
    if ((!xdata)||(!ydata)||(!zdata)||(!ychi)) error=102;
  }

  if ((error==0)&&((scanhun==1)||(scanten>2))&&((chicommand==11)||
    (chicommand==13)||(chicommand==31)||(chicommand==33)))
  { for (j=0;j<npoint;++j)
    { xa=pdata[j*12+5]; xb=pdata[j*12+6];
      sa=(float)fabs(xa+xb); sb=(float)fabs(xb);
      if ((sa<1.0e-20)||(sb<1.0e-20)) xc=0.0f;
      else if (scanten<3) xc=(xa-xb)/xb;
      else xc=(xa-xb)/(xa+xb);
      pdata[j*12+7]=xc;
    }
  }

  if ((error==0)&&(chicommand>0)&&(scanuni!=4)&&(scanuni!=8))
  { for (i=0;i<ncurve;++i)
    { k=(int)curvepar[i*12+7];
      if (i<ncurve-1) ndata=(int)curvepar[(i+1)*12+7]-k;
      else ndata=npoint-k;
      if (ndata>256) ndata=256;
      if ((scanuni==1)||(scanuni==5))
      { for (j=0;j<ndata;++j) xdata[j]=pdata[(j+k)*12];
        nfit=(int)((xdata[ndata-1]-xdata[0])*0.5+2);
        if (nfit>5) nfit=5; else if (nfit<3) nfit=3;
      }
      else if ((scanuni==2)||(scanuni==6))
      { for (j=0;j<ndata;++j) xdata[j]=pdata[(j+k)*12+1];
        if (xdata[ndata-1]-xdata[0]>50) nfit=4; else nfit=3;
      }
      else
      { for (j=0;j<ndata;++j) xdata[j]=pdata[(j+k)*12+2];
        nfit=1;
      }

      if ((error==0)&&(scanhun>1)&&(scanten<3)&&((chicommand==11)||
        (chicommand==13)||(chicommand==31)||(chicommand==33)))
      { for (j=0;j<ndata;++j) ydata[j]=xa=pdata[(j+k)*12+5];
        error=chicurve(ndata,nfit,xdata,ydata,ychi);
        for (j=0;(error==0)&&(j<ndata);++j)
        { xa=ydata[j]; xb=ychi[j]; sb=(float)fabs(1.0+xb);
          if (sb<1.0e-20) xc=0.0f;
          else xc=xa/(1.0f+xb);
          pdata[(j+k)*12+6]=xc; pdata[(j+k)*12+7]=xb;
        }
      }

      if ((error==0)&&((chicommand==12)||(chicommand==13)||
        (chicommand==32)||(chicommand==33)))
      { for (j=0;j<ndata;++j) ychi[j]=pdata[(j+k)*12+7];
        error=normcurve(ndata,ychi,ychi);
        if (error==0) for (j=0;j<ndata;++j) pdata[(j+k)*12+7]=ychi[j];
      }

      if ((error==0)&&((chicommand==22)||(chicommand==32)||
        (chicommand==33)))
      { for (j=0;j<ndata;++j) ychi[j]=pdata[(j+k)*12+8];
        error=normcurve(ndata,ychi,ychi);
        if (error==0) for (j=0;j<ndata;++j) pdata[(j+k)*12+8]=ychi[j];
      }
    }
  }
  else if ((error==0)&&(chicommand>0))
  { if ((scanhun>1)&&(scanten<3)&&(chicommand==11)||
      (chicommand==13)||(chicommand==31)||(chicommand==33))
    { for (i=0;i<ncurve;++i)
      { k=(int)curvepar[i*12+7];
        if (i<ncurve-1) ndata=(int)curvepar[(i+1)*12+7]-k;
        else ndata=npoint-k;
        j=p=0; patheta=sa=0.0f;
        while ((j<ndata)&&(p<256))
        { atheta=pdata[(j+k)*12+1];
          if (j==0)
          { patheta=atheta; sa=pdata[(j+k)*12+5];
            ++j; m=1;
          }
          else if (fabs(atheta-patheta)<0.01)
          { sa+=pdata[(j+k)*12+5];
            ++j; ++m;
          }
          if ((j==ndata-1)||(fabs(atheta-patheta)>=0.01))
          { xdata[p]=patheta; ydata[p]=sa/m;
            patheta=atheta; sa=pdata[(j+k)*12+5];
            ++p; ++j; m=1;
          }
        }
        m=ndata; ndata=p; nfit=0;
        if ((error==0)&&(ndata<10)) error=701+(datatype<<16);
        else if ((error==0)&&(xdata[ndata-1]-xdata[0]>50.0)) nfit=4;
        else if (error==0) nfit=3;
        if (error==0) error=chicurve(ndata,nfit,xdata,ydata,ychi);
        p=ndata; ndata=m;
        for (j=0;(error==0)&&(j<ndata);++j)
        { atheta=pdata[(j+k)*12+1];
          for (m=0;m<p;++m) if (fabs(atheta-xdata[m])<0.1) break;
          xa=0.0f;
          if (m>=p) error=901;
          else if (1.0+ychi[m]>1.0e-10) xa=ydata[m]/(1.0f+ychi[m]);
          if ((error==0)&&(xa>1.0e-10)) xb=pdata[(j+k)*12+5]/xa-1.0f;
          else xb=0.0f;
          if (xb>10.0) xb=10.0f; else if (xb<-10.0) xb=-10.0f;
          pdata[(j+k)*12+6]=xa; pdata[(j+k)*12+7]=xb;
        }
      }
    }

    if ((error==0)&&((chicommand==12)||(chicommand==13)||
      (chicommand==32)||(chicommand==33)))
    { for (i=0;i<ncurve;++i)
      { k=(int)curvepar[i*12+7];
        if (i<ncurve-1) ndata=(int)curvepar[(i+1)*12+7]-k;
        else ndata=npoint-k;
        sa=0.0f; for (j=0;j<ndata;++j) sa+=pdata[(j+k)*12+7];
        if (ndata>0) sa/=ndata;
        sb=0.0f;
        for (j=0;j<ndata;++j)
          sb+=(pdata[(j+k)*12+7]-sa)*(pdata[(j+k)*12+7]-sa);
        if (ndata>0) sb=(float)sqrt(sb/ndata);
        if (sb<1.0e-10) sb=1.0e-10f;
        for (j=0;j<ndata;++j)
          pdata[(j+k)*12+7]=(pdata[(j+k)*12+7]-sa)/sb;
      }
    }

    if ((error==0)&&((chicommand==22)||(chicommand==32)||
      (chicommand==33)))
    { for (i=0;i<ncurve;++i)
      { k=(int)curvepar[i*12+7];
        if (i<ncurve-1) ndata=(int)curvepar[(i+1)*12+7]-k;
        else ndata=npoint-k;
        sa=0.0f; for (j=0;j<ndata;++j) sa+=pdata[(j+k)*12+8];
        if (ndata>0) sa/=ndata;
        sb=0.0f;
        for (j=0;j<ndata;++j)
          sb+=(pdata[(j+k)*12+8]-sa)*(pdata[(j+k)*12+8]-sa);
        if (ndata>0) sb=(float)sqrt(sb/ndata);
        if (sb<1.0e-10) sb=1.0e-10f;
        for (j=0;j<ndata;++j)
          pdata[(j+k)*12+8]=(pdata[(j+k)*12+8]-sa)/sb;
      }
    }
  }
  if (xdata) delete [] xdata; if (ydata) delete [] ydata;
  if (zdata) delete [] zdata; if (ychi) delete [] ychi;

  return(error);
} //end of Pdintensity::chicalc

int Pdintensity::chicurve(int ndata,int nfit,float *xdata,float *ydata,
  float *ychi)
{ int i,j,k,ta,tb,tc,fitmath,trymax,trynum;
  float sa,sb,tolerance,reliable;
  float *xafit,*abest,*afit,*safit,*ytemp;

  xafit=abest=afit=safit=ytemp=NULL;
  fitmath=1; trymax=100; tolerance=0.1f;

  if (error==0)
  { xafit=new float [nfit]; abest=new float [nfit];
    afit=new float [nfit]; safit=new float [nfit];
    ytemp=new float [ndata];
    if ((!xafit)||(!abest)||(!afit)||(!safit)||(!ytemp)) error=102;
  }
  if ((error==0)&&(nfit>1))
  { tc=(int)((ndata-1.0)/(nfit-1.0)/2.0);
    for (i=0;i<nfit;++i)
    { k=(int)(i*(ndata-1.0)/(nfit-1.0));
      ta=k-tc; tb=k+tc;
      if (ta<0) ta=0;
      if (tb>ndata-1) tb=ndata-1;
      sa=0.0f; for (j=ta;j<tb;++j) sa+=ydata[j];
      afit[i]=sa/(tb-ta+1.0f); xafit[i]=xdata[k]; safit[i]=1.0f;
    }
    Polation aspline;
    aspline.init(nfit);
    if (error==0) error=aspline.loadcurve(xafit,afit);
    sa=aspline.getmemory();
    if ((error==0)&&(basemem<sa)) basemem=sa;
    for (i=0;(error==0)&&(i<ndata);++i)
      ytemp[i]=aspline.fspline(xdata[i]);
  }
  if ((error==0)&&(nfit>1))
  { sa=afit[0];
    for (i=0;i<nfit;++i) if (sa>afit[i]) sa=afit[i];
    if (sa<1.0e-20) sa=1.0e-20f;
    sb=ydata[0];
    for (i=0;i<ndata;++i) if (sb>ydata[i]) sb=ydata[i];
    if (sb<sa/10.0) sb=sa/10.0f;
    for (i=0;i<ndata;++i)
    { if (ytemp[i]<sb) ytemp[i]=ydata[i]/sb;
      else ytemp[i]=ydata[i]/ytemp[i];
    }

    tc=(int)((ndata-1.0)/(nfit-1.0)/2.0);
    for (i=0;i<nfit;++i)
    { k=(int)(i*(ndata-1.0)/(nfit-1.0));
      ta=k-tc; tb=k+tc;
      if (ta<0) ta=0;
      if (tb>ndata-1) tb=ndata-1;
      sa=0.0f; for (j=ta;j<tb;++j) sa+=ytemp[j];
      afit[i]=sa/(tb-ta+1.0f); xafit[i]=xdata[k]; safit[i]=1.0f;
    }
  }

  if ((error==0)&&(nfit>1))
  { Curvefit pdchicurve;
    pdchicurve.init(ndata,nfit);
    error=pdchicurve.loadcurve(fitmath,trymax,tolerance,xdata,ytemp,
      afit,safit,fitspline);
    sa=pdchicurve.getmemory();
    if ((error==0)&&(basemem<sa)) basemem=sa;
    if (error==0) error=pdchicurve.dofit(ychi,abest,&reliable,&trynum);
  }
  else if (error==0)
  { sa=0.0f; for (i=0;i<ndata;++i) sa+=ydata[i];
    if (ndata>0) sa/=ndata;
    for (i=0;i<ndata;++i)
    { ytemp[i]=ydata[i]; ychi[i]=sa;
    }
  }

  if (error==0)
  { for (i=0;i<ndata;++i)
    { if (ychi[i]<1.0e-10) ychi[i]=1.0e-10f;
        ychi[i]=ytemp[i]/ychi[i]-1.0f;
    }
  }

  if (xafit) delete [] xafit; if (abest) delete [] abest;
  if (afit) delete [] afit; if (safit) delete [] safit;
  if (ytemp) delete [] ytemp;

  return(error);
} //end of Pdintensity::chicurve

int Pdintensity::normcurve(int ndata,float *ydata,float *ychi)
{ int i;
  float xa,xb;

  if (error==0)
  { xa=0.0f; for (i=0;i<ndata;++i) xa+=ydata[i];
    if (ndata>0) xa/=ndata;
    xb=0.0f; for (i=0;i<ndata;++i) xb+=(ydata[i]-xa)*(ydata[i]-xa);
    if (ndata>0) xb=(float)sqrt(xb/ndata);
    if (xb<1.0e-10) xb=1.0e-10f;
    for (i=0;i<ndata;++i) ychi[i]=(ydata[i]-xa)/xb;
  }

  return(error);
} //end of Pdintensity::normcurve

int Pdintensity::reliability(float *areliable,float *breliable,
  float *creliable,float *dreliable,float *ereliable,
  float *freliable,float *greliable,float *hreliable)
{ int i,j,k,ndata;
  float xa,xb,xc,xd,xe,xf,xg,xh,ya,yb,yc,yd,ye,yf,yg,yh,yi,
    zd,zg,zh;
  
  xa=xb=xc=xd=xe=xf=xg=xh=0.0f; ya=yb=yc=yd=ye=yf=yg=yh=yi=0.0f;
  zd=zg=zh=0.0f;
  if (error==0)
  { for (i=0;i<ncurve;++i)
    { k=(int)curvepar[i*12+7];
      if (i<ncurve-1) ndata=(int)curvepar[(i+1)*12+7]-k;
      else ndata=npoint-k;
      for (j=0;j<ndata;++j)
      { zd=xd; zg=xg; zh=xh;
        xg=pdata[(j+k)*12+7];
        xh=pdata[(j+k)*12+8];
        xd=xg-xh; xe=xg*xg; xf=xh*xh;
        xa+=xd*xd*curvepar[i*12+5];
        xb+=(xe-xf)*curvepar[i*12+5];
        xc+=(xe+xf)*curvepar[i*12+5];

        ya+=(float)fabs(xd)*curvepar[i*12+5];
        yb+=(float)fabs(xh)*curvepar[i*12+5];
        yc+=xf*curvepar[i*12+5];

        if (j>0)
        { yd+=(float)fabs(xd-zd)*curvepar[i*12+5];
          ye+=(float)fabs(xh-zh)*curvepar[i*12+5];
          yf+=(xd-zd)*(xd-zd)*curvepar[i*12+5];
          yg+=(xh-zh)*(xh-zh)*curvepar[i*12+5];
          if (((xg>zg)&&(xh<zh))||((xg<zg)&&(xh>zh)))
            yh+=curvepar[i*12+5];
          yi+=curvepar[i*12+5];
        }
      }
    }
  }
  if ((error==0)&&(xc>1.0e-20)&&(areliable)&&(breliable))
  { *areliable=xa/xc; *breliable=xb/xc;
    if (*areliable>10.0) *areliable=10.0f;
    if (*breliable>10.0) *breliable=10.0f;
    else if (*breliable<-10.0) *breliable=-10.0f;
  }
  else if ((areliable)&&(breliable))
    *areliable=*breliable=0.0f;

  if ((error==0)&&(creliable)&&(dreliable)&&(ereliable)&&(freliable)
    &&(greliable)&&(hreliable))
  { if (yb>1.0e-20) *creliable=ya/yb; else *creliable=0.0f;
    if (yc>1.0e-20) *dreliable=xa/yc; else *dreliable=0.0f;
    if (ye>1.0e-20) *ereliable=yd/ye; else *ereliable=0.0f;
    if (yg>1.0e-20) *freliable=yf/yg; else *freliable=0.0f;
    if (yi>1.0e-20) *greliable=yh/yi; else *greliable=0.0f;
    *hreliable=(*creliable+*dreliable+*ereliable+*freliable+
      *greliable)*0.2f;
  }

  return(error);
} //end of Pdintensity::reliability

int Pdintensity::savecurve(char *filename,char *usermessage)
{ int i,j,k,ndata,begrow,linenum,scanhun,scanten,scanuni;
  float akout,atheta,aphi;

  if (error==0)
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { scanuni=datatype%10; scanten=(datatype/10)%10;
      scanhun=(datatype/100)%10; i=0; begrow=3;
      if (comment)
      { while (comment[i]!='\0')
        { if ((comment[i]=='\n')||(comment[i]=='\r')) ++begrow;
          ++i;
        }
      }
      else
      { while (usermessage[i]!='\0')
        { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
          ++i;
        }
      }

      if (ncurve<2) linenum=npoint;
      else linenum=0;
      if ((!comment)&&(linenum>0)) begrow+=7;
      else if (!comment) begrow+=5;
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(linenum,6); fileout.space(5);
      if (linenum>0)
        fileout.string("datakind beginning-row linenumbers",0,1);
      else fileout.string("datakind beginning-row multi-curves",0,1);
      if (comment) fileout.string(comment,0,1);
      else
      { fileout.string(usermessage,0,1);
        fileout.string("   angle-resolved photoemission extended fine");
        fileout.string(" structure (ARPEFS)",0,1);
        if (((scanuni==1)||(scanuni==5))&&((ndtheta==1)||(ndphi==1)))
        { fileout.string("   photoemission energy scan curves",0,16+1);
          fileout.string(
            "     parameters: curve point theta phi weightc weighte",
            0,1);
        }
        else if ((scanuni==1)||(scanuni==5))
        { fileout.string("   photoemission energy scan hologram",
            0,16+1);
          fileout.string(
            "     parameters: curve point theta phi weightc weighte",
            0,1);
        }
        else if ((scanuni==2)||(scanuni==6))
        { fileout.string("   photoemission polar scan curves",0,16+1);
          fileout.string(
            "     parameters: curve point wavenum phi weightc weighte",
            0,1);
        }
        else if ((scanuni==3)||(scanuni==7))
        { fileout.string("   photoemission azimuthal scan curves",
            0,16+1);
          fileout.string(
          "     parameters: curve point wavenum theta weightc weighte",
            0,1);
        }
        else
        { fileout.string("   photoelectron diffraction curves",0,16+1);
          fileout.string(
            "     parameters: curve point wavenum    weightc weighte",
            0,1);
        }
        fileout.string("     columns: ");
        if ((scanuni==1)||(scanuni==5)) fileout.string("wavenum ");
        else if ((scanuni==2)||(scanuni==6)) fileout.string("theta ");
        else if ((scanuni==3)||(scanuni==7)) fileout.string("phi ");
        else if ((scanuni==4)||(scanuni==8))
          fileout.string("theta phi ");
        if (scanten<3) fileout.string("intensity background chical");
        else if (scanten<5) fileout.string("intlcp intrcp asymcal");
        else if (scanten<7) fileout.string("intmup intmdown asymcal");
        else if (scanten<9) fileout.string("intsup intsdown asymcal");
        fileout.newline();
      }
      if ((!comment)||(linenum==0))
      { fileout.integer(ncurve,6); fileout.integer(npoint,6);
        fileout.integer(nkout,5); fileout.integer(ndtheta,5);
        fileout.integer(ndphi,5); fileout.integer(ndangle,5);
        fileout.string("  ncurve npoint nk ntheta nphi nangle",0,1);
      }

      if ((scanuni==1)||(scanuni==5))
      { for (i=0;i<ncurve;++i)
        { k=(int)curvepar[i*12+7];
          if (i<ncurve-1) ndata=(int)curvepar[(i+1)*12+7]-k;
          else ndata=npoint-k;
          if (ndata>1)
          { atheta=pdata[(k+1)*12+1]; aphi=pdata[(k+1)*12+2];
          }
          else
          { atheta=pdata[k*12+1]; aphi=pdata[k*12+2];
          }
          if ((!comment)||(linenum==0))
          { fileout.integer(i+1,4); fileout.integer(ndata,6);
            fileout.floating(atheta,7.1f,256);
            fileout.floating(aphi,7.1f,256);
            fileout.floating(curvepar[i*12+5],7.1f,256);
            fileout.floating(curvepar[i*12+6],7.1f,256);
            fileout.space(3); fileout.charfill('-',28,1);
          }
          for (j=0;j<ndata;++j)
          { fileout.floating(pdata[(j+k)*12],9.2f,256);
            if (scanhun<4)
            { fileout.floating(pdata[(j+k)*12+5],14.4f,512);
              fileout.floating(pdata[(j+k)*12+6],14.4f,512);
            }
            fileout.floating(pdata[(j+k)*12+7],14.4f,512+1);
          }
        }
      }
      else if ((scanuni==2)||(scanuni==6))
      { for (i=0;i<ncurve;++i)
        { k=(int)curvepar[i*12+7];
          if (i<ncurve-1) ndata=(int)curvepar[(i+1)*12+7]-k;
          else ndata=npoint-k;
          akout=pdata[k*12];
          if (ndata>1) aphi=pdata[(k+1)*12+2];
          else aphi=pdata[k*12+2];
          if ((!comment)||(linenum==0))
          { fileout.integer(i+1,4); fileout.integer(ndata,6);
            fileout.floating(akout,7.2f,256);
            fileout.floating(aphi,7.1f,256);
            fileout.floating(curvepar[i*12+5],7.1f,256);
            fileout.floating(curvepar[i*12+6],7.1f,256);
            fileout.space(3); fileout.charfill('-',28,1);
          }
          for (j=0;j<ndata;++j)
          { fileout.floating(pdata[(j+k)*12+1],9.1f,256);
            if (scanhun<4)
            { fileout.floating(pdata[(j+k)*12+5],14.4f,512);
              fileout.floating(pdata[(j+k)*12+6],14.4f,512);
            }
            fileout.floating(pdata[(j+k)*12+7],14.4f,512+1);
          }
        }
      }
      else if ((scanuni==3)||(scanuni==7))
      { for (i=0;i<ncurve;++i)
        { k=(int)curvepar[i*12+7];
          if (i<ncurve-1) ndata=(int)curvepar[(i+1)*12+7]-k;
          else ndata=npoint-k;
          akout=pdata[k*12];
          if (ndata>1) atheta=pdata[(k+1)*12+1];
          else atheta=pdata[k*12+1];
          if ((!comment)||(linenum==0))
          { fileout.integer(i+1,4); fileout.integer(ndata,6);
            fileout.floating(akout,7.2f,256);
            fileout.floating(atheta,7.1f,256);
            fileout.floating(curvepar[i*12+5],7.1f,256);
            fileout.floating(curvepar[i*12+6],7.1f,256);
            fileout.space(3); fileout.charfill('-',28,1);
          }
          for (j=0;j<ndata;++j)
          { fileout.floating(pdata[(j+k)*12+2],9.1f,256);
            if (scanhun<4)
            { fileout.floating(pdata[(j+k)*12+5],14.4f,512);
              fileout.floating(pdata[(j+k)*12+6],14.4f,512);
            }
            fileout.floating(pdata[(j+k)*12+7],14.4f,512+1);
          }
        }
      }
      else if ((scanuni==4)||(scanuni==8))
      { for (i=0;i<ncurve;++i)
        { k=(int)curvepar[i*12+7];
          if (i<ncurve-1) ndata=(int)curvepar[(i+1)*12+7]-k;
          else ndata=npoint-k;
          akout=pdata[k*12];
          if ((!comment)||(linenum==0))
          { fileout.integer(i+1,4); fileout.integer(ndata,6);
            fileout.floating(akout,7.2f,256);
            fileout.floating(curvepar[i*12+5],14.1f,256);
            fileout.floating(curvepar[i*12+6],7.1f,256);
            fileout.space(3); fileout.charfill('-',28,1);
          }
          for (j=0;j<ndata;++j)
          { fileout.floating(pdata[(j+k)*12+1],9.1f,256);
            fileout.floating(pdata[(j+k)*12+2],9.1f,256);
            if (scanhun<4)
            { fileout.floating(pdata[(j+k)*12+5],14.4f,512);
              fileout.floating(pdata[(j+k)*12+6],14.4f,512);
            }
            fileout.floating(pdata[(j+k)*12+7],14.4f,512+1);
          }
        }
      }
      if (error==0) error=fileout.geterror();
    }
    else error+=(datatype<<16);
  }

  return(error);
} //end of Pdintensity::savecurve

//make the emission angle by deviation
int Pdintensity::makemission(float adtheta,float adphi,float devstep)
{ int i;
  float atheta,aphi,btheta,bphi;
  Cartesia pdetecold,pdetecnew;

  if (error==0)
  { for (i=0;i<npoint;++i)
    { atheta=pdata[i*12+10]; aphi=pdata[i*12+11];
      if (devstep<9)
      { pdetecold.loadthetaphi(adtheta,adphi);
        pdetecnew=pdetecold.euler(aphi,atheta,0.0f);
        btheta=pdetecnew.theta(); bphi=pdetecnew.phi();
      }
      else
      { btheta=atheta+adtheta; bphi=aphi+adphi;
        if (btheta<0.0)
        { btheta=-btheta; bphi=bphi+180.0f;
        }
      }
      while (bphi>360.0) bphi-=360.0f;
      while (bphi<0.0) bphi+=360.0f;
      if (bphi>360.0-1.0e-3) bphi=0.0f;
      pdata[i*12+1]=btheta; pdata[i*12+2]=bphi;
    }
  }

  return(error);
} //end of Pdintensity::makemission

int Pdintensity::getfitdata(float *xdata,float *ydata)
{ int i,j;
  float swap;

  if (error==0)
  { for (i=0;i<npoint;++i)
    { xdata[i]=(float)i; ydata[i]=pdata[i*12];
    }
    for (i=0;i<npoint;++i)
    { for (j=i+1;j<npoint;++j)
      { if (ydata[i]>ydata[j])
        { swap=ydata[i]; ydata[i]=ydata[j]; ydata[j]=swap;
          swap=xdata[i]; xdata[i]=xdata[j]; xdata[j]=swap;
        }
      }
    }
    for (i=0;i<npoint;++i)
    { j=(int)xdata[i]; ydata[i]=pdata[j*12+8];
    }
  }
  return(error);
} //end of Pdintensity::getfitdata

int Pdintensity::getcurvepar(int icurve,int *ndata,float *akout,
  float *adtheta,float *adphi,float *weightc,float *weightk)
{ int j,k,scanuni;

  scanuni=datatype%10;
  if (error==0)
  { j=(int)curvepar[icurve*12+7];
    if (icurve<ncurve-1) k=(int)curvepar[(icurve+1)*12+7];
    else k=npoint;
    if (icurve<ncurve-1)
      *ndata=(int)(curvepar[(icurve+1)*12+7]-curvepar[icurve*12+7]);
    else if (icurve==ncurve-1)
      *ndata=npoint-(int)curvepar[icurve*12+7];
    *weightc=curvepar[icurve*12+5]; *weightk=curvepar[icurve*12+6];
  }
  if ((error==0)&&((scanuni==1)||(scanuni==5)))
  { *akout=curvepar[icurve*12];
    *adtheta=curvepar[icurve*12+3]; *adphi=curvepar[icurve*12+4];
    *adtheta=pdata[(k-1)*12+1]; *adphi=pdata[(k-1)*12+2];
  }
  else if ((error==0)&&((scanuni==2)||(scanuni==6)))
  { *adtheta=curvepar[icurve*12];
    *akout=curvepar[icurve*12+3]; *adphi=curvepar[icurve*12+4];
    *adphi=pdata[(k-1)*12+2];
  }
  else if ((error==0)&&((scanuni==3)||(scanuni==7)))
  { *adphi=curvepar[icurve*12];
    *akout=curvepar[icurve*12+3]; *adtheta=curvepar[icurve*12+4];
    *adtheta=pdata[(k-1)*12+1];
  }
  else if ((error==0)&&((scanuni==4)||(scanuni==8)))
  { *akout=curvepar[icurve*12+3];
    *adtheta=*adphi=0.0f;
  }
  else if (error==0) *akout=*adtheta=*adphi=0.0f;

  return(error);
} //end of Pdintensity::getcurvepar

