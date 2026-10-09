#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "cartesia.h"
#include "phase.h"
#include "radmat.h"
#include "meanpath.h"
#include "pdinten.h"
#include "fcomplex.h"
#include "vibrate.h"
#include "msfuncs.h"
#include "rotamat.h"
#include "pdchifit.h"
#include "jobtime.h"
#include "userutil.h"
#include "mscdrun.h"

Fcomplex Mscdrun::matrixelement(int ali,int alf,int am,float akin)
{ float xa,xb,xc,xd;
  Fcomplex cxa;

  if (error==0)
  { if ((ali==0)&&(alf==ali+1))
    { am=0;
      xb=(float)(((ali+1.0)*(ali+1.0)-am*am)/(ali+ali+1.0)/
        (ali+ali+3.0));
      xc=1.0f; xd=0.0f;
    }
    else if (ali<1) xb=xc=xd=0.0f;
    else if (alf==ali+1)
    { xb=(float)(((ali+1.0)*(ali+1.0)-am*am)/(ali+ali+1.0)/
        (ali+ali+3.0));
      xc=radmatrix->famphase(akin,1);
      xd=radmatrix->famphase(akin,2);
    }
    else if (alf==ali-1)
    { xb=(float)((ali*ali-am*am)/(ali+ali-1.0)/
        (ali+ali+1.0));
      xc=radmatrix->famphase(akin,3);
      xd=radmatrix->famphase(akin,4);
    }
    else xb=xc=xd=0.0f;
    if (xb>0.0) xa=xc*(float)sqrt(xb); else xa=0.0f;
    xd-=(float)(alf*90.0);
    cxa=xa*expix->fexpix(xd);
  }
  else cxa=0.0f;

  return(cxa);
} //end of Mscdrun::matrixelement

float Mscdrun::fitphotoemission(int afitmath,int ndata,
  int nafit,int *yafit,float *afit,float *safit,float *xdata,
  float *ydata,float *ymod,float *dyda)
{ int i,j,k;
  float xc,xd,reliable;
  float astep[20];

  if ((error==0)&&(nafit>20)) error=642;
  else if (error==0)
  { for (k=0;k<nafit;++k)
    { if ((afitmath==1)&&(yafit[k]>=0))
      { j=(int)fitvars[nafit*2+k]; i=(int)fitvars[nafit*3+k];
        if ((i>=totlayer)&&(j==2)&&(afit[k]<15.0)) astep[k]=0.5f;
        else if ((i>=totlayer)&&(j==2)&&(afit[k]<30.0))
          astep[k]=-afit[k]*0.25f;
        else if ((i>=totlayer)&&(j==2)) astep[k]=-afit[k]*0.75f;
        else if ((i>=totlayer)&&(j==3)&&(afit[k]>tsample))
          astep[k]=-0.1f*afit[k]*afit[k]/(tsample+0.1f*afit[k]);
        else if ((i>=totlayer)&&(j==3)&&(afit[k]>tsample/1.5))
          astep[k]=0.25f*afit[k]*afit[k]/(tsample-0.25f*afit[k]);
        else if ((i>=totlayer)&&(j==3))
          astep[k]=tsample/1.25f-afit[k];
        else astep[k]=afit[k]*0.01f+0.001f;
        if (fabs(astep[k])<1.0e-3) astep[k]=1.0e-3f;
      }
      else if ((afitmath==1)&&(yafit[k]>=-k))
      { j=yafit[k];
        afit[k]=afit[-j-1]*safit[k]; astep[k]=astep[-j-1]*safit[k];
      }
    }
    for (k=0;(error==0)&&(afitmath==1)&&(k<nafit);++k)
    { if (yafit[k]>=0)
      { afit[k]+=astep[k];
        for (j=k+1;j<nafit;++j)
          if (-yafit[j]==k+1) afit[j]=afit[k]*safit[j];
        for (j=0;j<nafit;++j) fitvars[j]=afit[j];
        error=intensity(afitmath,afit,xdata,ydata,ymod);
        for (i=0;i<ndata;++i) dyda[k*ndata+i]=ymod[i];
        afit[k]-=astep[k];
        for (j=k+1;j<nafit;++j)
          if (-yafit[j]==k+1) afit[j]=afit[k]*safit[j];
      }
    }
    if (error==0) error=intensity(afitmath,afit,xdata,ydata,ymod);
    if ((error==0)&&(afitmath==1))
    { for (k=0;k<nafit;++k)
      { if (yafit[k]>=0)
        { for (i=0;i<ndata;++i)
            dyda[k*ndata+i]=(dyda[k*ndata+i]-ymod[i])/astep[k];
        }
        else
        { for (i=0;i<ndata;++i) dyda[k*ndata+i]=0.0f;
        }
      }
    }
  }

  xc=xd=0.0f;
  if (error==0)
  { for (i=0;i<ndata;++i)
    { xc+=(ymod[i]-ydata[i])*(ymod[i]-ydata[i]);
      xd+=ymod[i]*ymod[i]+ydata[i]*ydata[i];
    }
  }
  if ((error==0)&&(xd>1.0e-10))
  { reliable=xc/xd;
    if (reliable>10.0) reliable=10.0f;
  }
  else if (error==0) reliable=0.0f;
  else reliable=(float)error;

  return(reliable);
} //end of Mscdrun::fitphotoemission

int Mscdrun::getfitdata(float *xdata,float *ydata,float *afit,
  float *safit)
{ int k;
  if (error==0) error=pdintensity->getfitdata(xdata,ydata);
  for (k=0;(error==0)&&(k<nfit);++k)
  { afit[k]=fitvars[k]; safit[k]=fitvars[nfit+k];
  }
  return(error);
} //end of Mscdrun::getfitdata

int Mscdrun::getnumpoint()
{ return(npoint);
} //end of Mscdrun::getnumpoint

int Mscdrun::getnumcurve()
{ return(ncurve);
} //end of Mscdrun::getnumcurve

int Mscdrun::getnumfit()
{ return(nfit);
} //end of Mscdrun::getnumfit

int Mscdrun::fitarpefs(Mscdrun *mscdrun)
{ int k;
  float xa,reliable;
  float *xdata,*ydata,*ybest,*abest,*afit,*safit;

  xdata=ydata=ybest=abest=afit=safit=NULL;
  if (error==0)
  { xdata=new float [npoint]; ydata=new float [npoint];
    ybest=new float [npoint]; abest=new float [nfit];
    afit=new float [nfit]; safit=new float [nfit];
  }
  if ((error==0)&&((!xdata)||(!ydata)||(!ybest)||(!abest)||(!afit)||
    (!safit))) error=102;
  if (error==0) error=pdintensity->getfitdata(xdata,ydata);
  for (k=0;(error==0)&&(k<nfit);++k)
  { afit[k]=fitvars[k]; safit[k]=fitvars[nfit+k];
  }
  trynum=0;
  if (error==0)
  { Pdchifit arpefs;
    arpefs.init(npoint,nfit);
    if (error==0) error=arpefs.loadcurve(fitmath,trymax,ftolerance,
      xdata,ydata,afit,safit,mscdrun,fitpdintensity);
    xa=arpefs.getmemory();
    if (basemem<xa) basemem=xa;
    if (error==0) error=arpefs.dofit(ybest,abest,&reliable,&trynum);
  }

  if (xdata) delete [] xdata; if (ydata) delete [] ydata;
  if (ybest) delete [] ybest; if (abest) delete [] abest;
  if (afit) delete [] afit; if (safit) delete [] safit;

  return(error);
} //end of Mscdrun::fitarpefs

int Mscdrun::savecurve(char *filename,char *usermessage)
{ int i,j,k,ndata,begrow,linenum,scanuni,scanten,afitmath;
  float xa,akout,adtheta,adphi,intena,intenb,chical,chiexp,weightc,
    weightk,areliable,breliable;
  char timestring[15];

  scanuni=datatype%10; scanten=(datatype/10)%10;
  if (error==0)
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { i=0; begrow=20;
      while (usermessage[i]!='\0')
      { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
        ++i;
      }
      if (ncurve<2) linenum=npoint;
      else linenum=0;
      if (linenum>0) begrow+=2;
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(linenum,6); fileout.space(5);
      if (linenum>0)
        fileout.string("datakind beginning-row linenumbers",0,1);
      else fileout.string("datakind beginning-row multi-curves",0,1);
      fileout.string(usermessage,0);
      fileout.string(" angle-resolved photoemission extended fine");
      fileout.string(" structure (ARPEFS)",0,1);
      fileout.string(" multiple scattering calculation of ");
      fileout.string(sysname,0,1);
      fileout.string(" calculated by ");
      fileout.string(username); fileout.string(" on ");
      jobtime->endtime(timestring,13);
      fileout.string(timestring,0,1);
      fileout.string("   initial angular momentum (li) = ",0,16);
      fileout.integer(linitial); fileout.string("   msorder= ");
      fileout.integer(msorder,2); fileout.string("   raorder=");
      fileout.integer(raorder,2,1);
      fileout.string(
        "   photon polarization angle (polar,azimuth) =(");
      fileout.floating(ltheta,6.1f,256);
      fileout.floating(lphi,6.1f,256);
      fileout.string(" ) (deg)",0,1);
      fileout.string("   radius, depth and lattice constant=",0,16);
      fileout.floating(radius,6.1f,256); fileout.string(",");
      fileout.floating(depth,6.1f,256); fileout.string(" and ");
      fileout.floating(lattice,6.2f,256);
      fileout.string(" angstrom",0,1);
      fileout.string("   cluster size="); fileout.integer(natoms,4);
      fileout.string(" atoms and spacings=");
      if (totlayer>0) fileout.floating(layorig[2]*lattice,6.2f,256);
      if (totlayer>1) fileout.floating(layorig[4+2]*lattice,6.2f,256);
      if (totlayer>2) fileout.floating(layorig[8+2]*lattice,6.2f,256);
      if (totlayer>3) fileout.floating(layorig[12+2]*lattice,6.2f,256);
      if (totlayer>3) fileout.string(" angs",0,1);
      else fileout.string(" angstrom",0,1);
      fileout.string("   inner potential=");
      fileout.floating(vinner,5.1f,256);
      fileout.string(" V  debye and sample temperature=");
      fileout.integer((int)tdebye,4); fileout.string(",");
      fileout.integer((int)tsample,4); fileout.string(" K",0,1);
      if (valence>1.0-1.0e-3)
      { fileout.string("   number of valence electrons=");
        fileout.integer((int)valence,4);
        fileout.string("     bandgap energy=");
        fileout.floating(bandgap,7.2f,256);
        fileout.string(" eV",0,1);
      }
      else if (valence>0.1)
      { fileout.string("   electron attenuation length=");
        fileout.floating(bandgap,8.4f,256);
        fileout.string("*energy^");
        fileout.floating(valence,5.3f,256);
        fileout.string(" angstrom",0,1);
      }
      else
      { fileout.string("   electron wave attenuation due to inelastic");
        fileout.string(" process not considered",0,1);
      }
      fileout.string("   density of bulk=");
      fileout.floating(density,7.2f,256);
      fileout.string(" g/cm3     molecular weight=");
      fileout.floating(mweight,7.2f,256);
      fileout.string(" amu",0,1);
      if (katoms<2)
      { fileout.string("   effective weight=");
        fileout.floating(aweight[0],7.1f,256);
        fileout.string(" amu",0,1);
      }
      else
      { fileout.string("   effective weight for kind 1-");
        fileout.integer(katoms); fileout.string(" =");
        for (j=0;j<katoms;++j) fileout.floating(aweight[j],7.1f,256);
        fileout.string(" amu",0,1);
      }
      fileout.string("   half aperture angle= ");
      fileout.floating(accepang,6.1f,256); fileout.string(" deg");
      fileout.space(12); fileout.string("pathcut= ");
      if ((pathcut<1.0e-10)||(pathcut>1.0e-4))
        fileout.floating(pathcut,10.4f,256+1);
      else fileout.floating(pathcut,10.1f,512+1);


      if (((scanuni==1)||(scanuni==5))&&((ndtheta==1)||(ndphi==1)))
        fileout.string("   photoemission energy scan curves",0,16+1);
      else if ((scanuni==1)||(scanuni==5))
        fileout.string("   photoemission energy scan hologram",0,16+1);
      else if ((scanuni==2)||(scanuni==6))
        fileout.string("   photoemission polar scan curves",0,16+1);
      else if ((scanuni==3)||(scanuni==7))
        fileout.string("   photoemission azimuthal scan curves",0,16+1);
      else
        fileout.string("   photoelectron diffraction curves",0,16+1);
      fileout.string(
        "     parameters: curve point theta phi weightc weighte",0,1);
      fileout.string("     columns: ");
      if ((scanuni==1)||(scanuni==5)) fileout.string("wavevec ");
      else if ((scanuni==2)||(scanuni==6)) fileout.string("theta ");
      else if ((scanuni==3)||(scanuni==7)) fileout.string("phi ");
      else if ((scanuni==4)||(scanuni==8)) fileout.string("theta phi ");
      if (scanten<3) fileout.string("intensity background chical");
      else if (scanten<5) fileout.string("intlcp intrcp asymcal");
      else if (scanten<7) fileout.string("intmup intmdown asymcal");
      else if (scanten<9) fileout.string("intsup intsdown asymcal");
      if ((fitmode>0)&&(scanten<3)) fileout.string(" chiexp");
      else if ((fitmode>0)&&(scanten<9)) fileout.string(" asymexp");
      fileout.newline();
      fileout.integer(ncurve,6); fileout.integer(npoint,6);
      fileout.integer(nkout,5); fileout.integer(ndtheta,5);
      fileout.integer(ndphi,5); fileout.integer(ndangle,5);
      fileout.string("  ncurve npoint nk ntheta nphi nangle",0,1);
      if (error==0) error=fileout.geterror();
    }
    k=ndata=0;
    for (i=0;(error==0)&&(i<ncurve);++i)
    { if (i>0) k+=ndata;
      error=pdintensity->getcurvepar(i,&ndata,&akout,&adtheta,&adphi,
        &weightc,&weightk);
      if ((error==0)&&((scanuni==1)||(scanuni==5)))
      { fileout.integer(i+1,4); fileout.integer(ndata,6);
        fileout.floating(adtheta,7.1f,256);
        fileout.floating(adphi,7.1f,256);
        fileout.floating(weightc,7.1f,256);
        fileout.floating(weightk,7.1f,256);
        fileout.space(3); fileout.charfill('-',28,1);
      }
      else if ((error==0)&&((scanuni==2)||(scanuni==6)))
      { fileout.integer(i+1,4); fileout.integer(ndata,6);
        fileout.floating(akout,7.2f,256);
        fileout.floating(adphi,7.1f,256);
        fileout.floating(weightc,7.1f,256);
        fileout.floating(weightk,7.1f,256);
        fileout.space(3); fileout.charfill('-',28,1);
      }
      else if ((error==0)&&((scanuni==3)||(scanuni==7)))
      { fileout.integer(i+1,4); fileout.integer(ndata,6);
        fileout.floating(akout,7.2f,256);
        fileout.floating(adtheta,7.1f,256);
        fileout.floating(weightc,7.1f,256);
        fileout.floating(weightk,7.1f,256);
        fileout.space(3); fileout.charfill('-',28,1);
      }
      else if ((error==0)&&((scanuni==4)||(scanuni==8)))
      { fileout.integer(i+1,4); fileout.integer(ndata,6);
        fileout.floating(akout,7.2f,256);
        fileout.floating(weightc,14.1f,256);
        fileout.floating(weightk,7.1f,256);
        fileout.space(3); fileout.charfill('-',28,1);
      }
      for (j=0;(error==0)&&(j<ndata);++j)
      { error=pdintensity->getpoint(j+k,&akout,&adtheta,&adphi,
          &xa,&xa,&intena,&intenb,&chical,&chiexp);
        if ((error==0)&&((scanuni==1)||(scanuni==5)))
        { fileout.floating(akout,7.2f,256);
          fileout.floating(intena,14.4f,512);
          fileout.floating(intenb,14.4f,512);
          fileout.floating(chical,14.4f,512);
        }
        else if ((error==0)&&((scanuni==2)||(scanuni==6)))
        { fileout.floating(adtheta,7.1f,256);
          fileout.floating(intena,14.4f,512);
          fileout.floating(intenb,14.4f,512);
          fileout.floating(chical,14.4f,512);
        }
        else if ((error==0)&&((scanuni==3)||(scanuni==7)))
        { fileout.floating(adphi,7.1f,256);
          fileout.floating(intena,14.4f,512);
          fileout.floating(intenb,14.4f,512);
          fileout.floating(chical,14.4f,512);
        }
        else if ((error==0)&&((scanuni==4)||(scanuni==8)))
        { fileout.floating(adtheta,7.1f,256);
          fileout.floating(adphi,7.1f,256);
          fileout.floating(intena,14.4f,512);
          fileout.floating(intenb,14.4f,512);
          fileout.floating(chical,14.4f,512);
        }
        if ((error==0)&&(fitmode>0))
          fileout.floating(chiexp,14.4f,512);
        if (error==0) fileout.newline();
      }
      if (error==0) error=fileout.geterror();
    }

    if ((error==0)&&(fitmode>0))
      error=pdintensity->reliability(&areliable,&breliable);
    if ((error==0)&&(fitmode>0))
    { fileout.string("fitted parameters ( nfit = ",0,16);
      fileout.integer(nfit,2); fileout.string("  nfollow = ");
      fileout.integer(nfollow,2); fileout.string("  )",0,1);
      fileout.string("   r-factora = ",0,16);
      fileout.floating(areliable,9.4f,256);
      fileout.string("   r-factorb = ");
      fileout.floating(breliable,9.4f,256+2);
      fileout.integer(0,5);
      fileout.floating(vinner,9.2f,256);
      fileout.floating(tdebye,9.1f,256);
      fileout.floating(lattice,9.3f,256); fileout.space(11);
      fileout.string("general vinner tdebye lattice",0,1);
      for (i=0;i<totlayer;++i)
      { fileout.integer(i+1,5);
        fileout.floating(layorig[i*4+2],9.5f,256);
        fileout.floating(layorig[i*4],9.5f,256);
        fileout.floating(layorig[i*4+3]*laycell[i*4],9.5f,256);
        fileout.floating(layorig[i*4+3]*laycell[i*4+2],9.5f,256);
        fileout.string("  layer spacing length units",0,1);
      }
      fileout.space(43); fileout.string("unit: lattice constant",0,2);
      fileout.string("   r-factora=sum((chic-chie)*(chic-chie))");
      fileout.string("/(sum(chic*chic+chie*chie)",0,1);
      fileout.string("   r-factorb=sum(chic*chic-chie*chie)");
      fileout.string("/sum(chic*chic+chie*chie)",0,1);
      if (error==0) error=fileout.geterror();
    }

    if ((error==0)&&(nfit>0)&&(fitmode>1))
    { fileout.string("fitting history (",0,16);
      fileout.integer(trynum,3); fileout.string(" trials )",0,2);
      if (trynum<trymax+50) k=trynum; else k=trymax+50;
      for (i=0;i<k;++i)
      { afitmath=(int)fithist[i*(nfit+10)+nfit+2];
        fileout.space(3); fileout.integer(i+1,3);
        fileout.string("  factors = ");
        fileout.floating(fithist[i*(nfit+10)+nfit],9.4f,256);
        fileout.floating(fithist[i*(nfit+10)+nfit+1],9.4f,256);
        fileout.space(20);
        if (afitmath>2) fileout.string("Netsearch",0,1);
        else if (afitmath==2) fileout.string("Downhill",0,1);
        else if (afitmath==1) fileout.string("Marquardt",0,1);
        else fileout.newline();
        fileout.space(6); fileout.string("  fitvars = ");
        for (j=0;j<nfit;++j)
          fileout.floating(fithist[i*(nfit+10)+j],13.7f,256);
        fileout.newline();
      }
      if (error==0) error=fileout.geterror();
    }

    if (error==0)
    { error=jobtime->timestamp(fileout);
      if (error!=0) error+=(datatype<<16);
    }
  }

  if ((error==0)&&(mype==0)&&(tusage)&&(dispmode>4))
  { Textout conout;
    conout.string(
      "processor time distribution (in seconds or percent)",0,32+2);
    conout.string(" pid computation   sending  receiving      idle");
    conout.string(" comput send receiv  idle",0,1);
    for (i=0;i<numpe+1;++i)
    { if (i<numpe)
      { k=(int)tusage[i*10]; conout.integer(k,4);
      }
      else conout.string(" sum");
      for (j=1;j<9;++j)
      { if (j<5) conout.floating(tusage[i*10+j],11.3f,512);
        else conout.floating(tusage[i*10+j],6.1f,256);
      }
      conout.newline();
    }
    conout.string("Job ",0,16); conout.integer(jobnum+1,3);
    conout.string(" of "); conout.integer(jobtotal,3);
    conout.string(" total time of "); conout.integer(numpe,3);
    if (numpe<2) conout.string(" processor = ");
    else conout.string(" processors = ");
    conout.floating(tusage[numpe*10],12.4f,512);
    conout.string(" seconds",0,1);
  }
  if ((error==0)&&(mype==0)&&(displog>0)&&(flogout))
  { flogout->string(
      "processor time distribution (in seconds or percent)",0,32+2);
    flogout->string(" pid computation   sending  receiving      idle");
    flogout->string(" comput send receiv  idle",0,1);
    for (i=0;i<numpe+1;++i)
    { if (i<numpe)
      { k=(int)tusage[i*10]; flogout->integer(k,4);
      }
      else flogout->string(" sum");
      for (j=1;j<9;++j)
      { if (j<5) flogout->floating(tusage[i*10+j],11.3f,512);
        else flogout->floating(tusage[i*10+j],6.1f,256);
      }
      flogout->newline();
    }
    flogout->string("Job ",0,16); flogout->integer(jobnum+1,3);
    flogout->string(" of "); flogout->integer(jobtotal,3);
    flogout->string(" total time of "); flogout->integer(numpe,3);
    if (numpe<2) flogout->string(" processor = ");
    else flogout->string(" processors = ");
    flogout->floating(tusage[numpe*10],12.4f,512);
    flogout->string(" seconds",0,1);
  }

  return(error);
} //end of Mscdrun::savecurve

float fitpdintensity(int fitmath,int ndata,int nfit,int *yafit,
  float *afit,float *safit,float *xdata,float *ydata,float *ymod,
  float *dyda,Mscdrun *mscdrun)
{ float reliable;
  reliable=mscdrun->fitphotoemission(fitmath,ndata,nfit,yafit,
    afit,safit,xdata,ydata,ymod,dyda);
  return(reliable);
} //end of fitpdintensity

