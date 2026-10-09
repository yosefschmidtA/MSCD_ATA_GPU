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

int Mscdrun::summation(float akin,float *xdetec,float *polaron,
  float *suminten,float *bakinten,Fcomplex *asum,Fcomplex *bsum,
  Fcomplex *csum)
{ int ia,ib,ic,id,j,k,m,p,q,t,ali,am,alf,mevadd,megadd,evedim,
    eegdim;
  float xa,xb,xc,emiter;
  Fcomplex cxa,cxb,cxc,dsum,esum;
  int lamda[64];
  Fcomplex algam[256];

  for (j=0;(error==0)&&(j<32);++j)
  { if ((j==0)||(j==3)||(j==10)) k=0;
    else if (j<3) k=3-j*2;
    else if (j<6) k=18-j*4;
    else if (j<8) k=13-j*2;
    else if (j<10) k=51-j*6;
    else if (j<13) k=46-j*4;
    else if (j<15) k=108-j*8;
    else k=0;

    if (j==10) m=2;
    else if ((j==3)||(j==6)||(j==7)||(j==11)||(j==12)) m=1;
    else m=0;

    lamda[j]=k; lamda[32+j]=m;
  }


  ali=linitial;
  for (ib=0;(error==0)&&(ib<natoms);++ib)
  { for (ic=0;ic<natoms;++ic)
    { for (j=0;j<radim;++j)
      { id=ib*natoms*radim+ic*radim+j;
        if ((msorder>0)&&(ic!=ib)) asum[id]=devendetec[id];
        else asum[id]=0.0f;
      }
    }
  }

  for (m=msorder;(error==0)&&(m>=2);--m)
  { for (ib=0;ib<natoms;++ib)
    { for (ic=0;ic<natoms;++ic)
      { if (ic==ib) continue;
        for (j=0;j<radim;++j)
        { id=ib*natoms*radim+ic*radim+j;
          bsum[id]=asum[id];
        }
      }
    }
    for (ia=0;ia<natoms;++ia)
    { emiter=patom[ia*12+7];
      if ((m==2)&&(emiter==0)) continue;
      for (ib=0;ib<natoms;++ib)
      { if (ib==ia) continue;
        for (j=0;j<radim;++j)
        { id=ia*natoms*radim+ib*radim+j;
          csum[j]=devendetec[id];
        }

        if (tevencut[(m-1)*natoms*natoms+ia*natoms+ib]==0) continue;
        for (ic=0;ic<natoms;++ic)
        { id=ia*natoms*natoms+ib*natoms+ic;
          evedim=tevendim[id];
          if ((sizeint<4)&&(m>5)) evedim>>=12;
          else if (m>8) evedim>>=24;
          else evedim>>=(m-2)*4;
          evedim&=15;
          if ((ic==ib)||(evedim<1)) continue;
          k=tevenadd[id]; eegdim=(int)tevenpar[k*10+5];
          megadd=ib*natoms*radim+ic*radim;
          mevadd=(int)tevenpar[k*10+6];

          if (evedim<1) continue;
          else if (evedim==1)
            csum[0]+=bsum[megadd]*tevenelem[mevadd];
          else if ((evedim<16)&&(eegdim<16))
          { xa=talpha[id]; xb=tgamma[id];
            for (j=0;j<evedim;++j)
            { for (k=0;k<evedim;++k)
              { t=j*evedim+k;
                p=lamda[j]; q=lamda[k];
                if ((p==0)&&(q==0))
                  algam[t]=1.0;
                else if ((p==0)&&(q<0))
                  algam[t]=conj(algam[t-1]);
                else if ((p<0)&&(q==0))
                  algam[t]=conj(algam[t-evedim]);
                else if ((p<0)&&(q>0))
                  algam[t]=conj(algam[t-evedim+1]);
                else if ((p<0)&&(q<0))
                  algam[t]=conj(algam[t-evedim-1]);
                else
                { xc=-p*xb-q*xa;
                  algam[t]=expix->fexpix(xc);
                }
                csum[j]+=algam[t]*bsum[megadd+k]*
                  tevenelem[mevadd+j*eegdim+k];
              }
            }
          }
          else error=901;
        }
        for (j=0;j<radim;++j)
          asum[ia*natoms*radim+ib*radim+j]=csum[j];
      }
    }
  }

  for (ia=0;ia<natoms;++ia)
  { emiter=patom[ia*12+7];
    if (emiter==0) continue;
    for (am=-ali;am<=ali;++am)
    { m=am+ali+1;
      dsum=esum=0.0f;
      for (alf=ali-1;alf<=ali+1;alf+=2)
      { if ((alf<0)||(am<-alf)||(am>alf)||((finals==1)&&(alf==ali-1))||
          ((finals==2)&&(alf==ali+1))) continue;
        cxa=0.0f;
        for (ib=0;ib<natoms;++ib)
        { if ((msorder<1)||(finals==3)||(ib==ia)) continue;
          error=onevenemit(ia,ib,alf,am,akin,xdetec,polaron,csum);
          id=ia*natoms*radim+ib*radim;
          for (j=0;j<radim;++j) cxa+=asum[id+j]*csum[j];
        }
        if (finals==4) cxb=0.0f;
        else cxb=onemidetec(akin,ia,alf,am,xdetec,polaron);
        cxc=matrixelement(ali,alf,am,akin);
        dsum+=cxa*cxc; esum+=cxb*cxc;
      }
      *suminten+=emiter*norm(dsum+esum); *bakinten+=emiter*norm(esum);
    }
  }

  return(error);
} //end of Mscdrun::summation

int Mscdrun::intensity(int afitmath,float *afit,float *xdata,
  float *ydata,float *ymod)
{ int i,j,k,n,scanuni;
  float xa,xb,xc,ya,akin,akout,bkout,adtheta,adphi,suminten,
    bakinten,netinten,altheta,alphi,bdtheta,bdphi,areliable,breliable;
  float xdetec[6],polaron[10];
  Fcomplex *asum,*bsum,*csum;

  asum=bsum=csum=NULL;
  if ((error==0)&&(mype==0)&&((!xdata)||(!ydata)||(!ymod))) error=901;
  else if ((error==0)&&(mype==0)&&(nfit>0)&&(!fithist)) error=901;
  else if (error==0)
  { asum=new Fcomplex [natoms*natoms*radim];
    bsum=new Fcomplex [natoms*natoms*radim];
    csum=new Fcomplex [radim];
    if ((!asum)||(!bsum)||(!csum)) error=102;
    else if (mype==0) error=dispintensity(1);
  }
  
  if (mype==0)
  { if (error==0) for (i=0;i<npoint;++i) pdnum[i]=(int)xdata[i];
    if (error==0) for (n=0;n<nfit;++n) fitvars[n]=afit[n];
    if ((error==0)&&(nfit>0)&&(fitmode==2))
    { for (n=0;n<nfit;++n)
        if (trynum<trymax+50) fithist[trynum*(nfit+10)+n]=fitvars[n];
      bdtheta=fitvars[0]; bdphi=fitvars[1];
      error=pdintensity->makemission(bdtheta,bdphi,(float)devstep);
      if (error==0) error=pdintensity->setbeamangle(ltheta,lphi);
    }
    else if ((error==0)&&(nfit>0)&&(fitmode>3))
    { for (n=0;n<nfit;++n)
      { if (trynum<trymax+50) fithist[trynum*(nfit+10)+n]=fitvars[n];
        j=(int)fitvars[nfit*2+n]; i=(int)fitvars[nfit*3+n];
        layorig[i*4+j]=fitvars[n];
      }
      lattice=layorig[totlayer*4]; vinner=layorig[totlayer*4+2];
      tdebye=layorig[totlayer*4+3];

      if ((fitmode==4)||(fitmode==6))
        error=vibrate->loadparameter(density,mweight,tdebye,tsample);
    }
    if ((error==0)&&((fitmode==5)||(fitmode==6)))
      error=makeatoms(0.0f);
    if ((trynum==0)||(fitmode==5)||(fitmode==6))
    { if (error==0) error=maketripar();
      if (error==0) error=makedblpar();
      if (error==0) error=allrotation();
    }

    if (error==0) error=dispintensity(2);
    if (error==0) error=dispintensity(3,-1);
  }
  if ((error==0)&&(numpe>1)&&(mype==0)) error=sendsup(1);
  else if ((error==0)&&(numpe>1))
  { error=receivesup(1);
    if ((error==0)&&(nfit>0)&&(fitmode==2))
    { bdtheta=fitvars[0]; bdphi=fitvars[1];
      error=pdintensity->makemission(bdtheta,bdphi,(float)devstep);
      if (error==0) error=pdintensity->setbeamangle(ltheta,lphi);
    }
    if ((error==0)&&((fitmode==4)||(fitmode==6)))
      error=vibrate->loadparameter(density,mweight,tdebye,tsample);
  }

  bkout=akin=0.0f; k=scanuni=0;
  for (i=pdbeg;(error==0)&&(i<pdend);++i)
  { j=pdnum[i];
    error=pdintensity->getpoint(j,&akout,&adtheta,&adphi,&altheta,
      &alphi,&xa,&xb,&xc,&ya);
    if (i==pdbeg) bkout=akout+1.0f;
    if ((error==0)&&(bkout!=akout))
    { bkout=akout; akin=kinside(akout,vinner);
      error=alltrievent(0,akin);
    }
    suminten=bakinten=0.0f;
    for (k=0;(error==0)&&(k<5);++k)
    { thetainside(k,akin,akout,adtheta,adphi,altheta,alphi,xdetec,
        polaron);
      if (k==0) error=alldblevent(akin,xdetec);
      if (error==0) error=allevendetec(akin,xdetec);
      if (error==0) error=summation(akin,xdetec,polaron,
        &suminten,&bakinten,asum,bsum,csum);
      if ((k==0)&&(accepang>=1.0e-3))
      { suminten+=suminten; bakinten+=bakinten;
      }
      else if (accepang<1.0e-3) break;
    }
    if ((error==0)&&(k>0))
    { suminten/=float(k+1.0); bakinten/=float(k+1.0);
    }
    if (bakinten<1.0e-10) netinten=0.0f;
    else netinten=suminten/bakinten-1.0f;
    if (netinten>10.0) netinten=10.0f;
    else if (netinten<-10.0) netinten=-10.0f;
    if (error==0) error=pdintensity->loadpoint(j,akout,adtheta,adphi,
      altheta,alphi,suminten,bakinten,netinten,ya);

    if ((error==0)&&(mype==0)) error=dispintensity(4,i,0,0.0f,0.0f,
      akout,adtheta,adphi,suminten);
  }
  if ((error==0)&&(numpe>1)&&(mype==0)) error=receivesup(2);
  else if ((error==0)&&(numpe>1)) error=sendsup(2);

  if (mype==0)
  { if (error==0) scanuni=datatype%10;
    if ((error==0)&&(scanuni<5)) k=11;
    else if (error==0) k=13;
    if ((error==0)&&(scanmode>200)) error=pdintensity->chicalc(k);
    for (i=0;(error==0)&&(i<npoint);++i)
    { j=pdnum[i];
      error=pdintensity->getpoint(j,&akout,&adtheta,&adphi,&altheta,
        &alphi,&xa,&xb,&xc,&ya);
      if (error==0)
      { suminten=xa; netinten=xc; ymod[i]=netinten;
      }
      if ((error==0)&&((i<pdbeg)||(i>=pdend)))
        error=dispintensity(4,i,0,0.0f,0.0f,akout,adtheta,adphi,
          suminten);
    }

    if ((error==0)&&(fitmode>0))
      error=pdintensity->reliability(&areliable,&breliable);
    else areliable=breliable=0.0f;
    if ((error==0)&&(nfit>0)&&(trynum<trymax+50))
    { fithist[trynum*(nfit+10)+nfit]=areliable;
      fithist[trynum*(nfit+10)+nfit+1]=breliable;
      fithist[trynum*(nfit+10)+nfit+2]=(float)afitmath;
    }
    if (error==0) dispintensity(5,0,afitmath,areliable,breliable);
    if (error==0) ++trynum;
  }

  if (asum) delete [] asum; if (bsum) delete [] bsum;
  if (csum) delete [] csum;

  return(error);
}
//end of Mscdrun::intensity

int Mscdrun::dispintensity(int messagenum,int pointnum,int afitmath,
  float areliable,float breliable,float akout,float adtheta,
  float adphi,float suminten)
{ int i;
  float job;

  if ((error==0)&&(npoint>0))
    job=(float)((pointnum+1.0)*100.0/npoint);
  else job=0.0f;
  if ((error==0)&&(messagenum==1)&&(dispmode>1)&&(trynum==0))
  { Textout conout;
    if (msorder>1) conout.string("Multiple scattering",0,16);
    else if (msorder==1) conout.string("Single scattering",0,16);
    else conout.string("Photoemission",0,16);
    conout.string(" calculation (msorder=");
    conout.integer(msorder); conout.string(" raorder=");
    conout.integer(raorder); conout.string(")",0,1);
  }

  if ((error==0)&&(displog>0)&&(messagenum==1)&&(trynum==0)&&(flogout))
  { if (msorder>1) flogout->string("Multiple scattering",0,16);
    else if (msorder==1) flogout->string("Single scattering",0,16);
    else flogout->string("Photoemission",0,16);
    flogout->string(" calculation (msorder=");
    flogout->integer(msorder); flogout->string(" raorder=");
    flogout->integer(raorder); flogout->string(")",0,1);
    if (error==0) error=flogout->geterror();
  }

  if ((error==0)&&(messagenum==2)&&(dispmode>2)&&(fitmode==2))
  { Textout conout;
    conout.space(8,16); conout.floating(fitvars[0],8.1f,256);
    conout.floating(fitvars[1],8.1f,256); conout.space(11);
    conout.string("deviation theta and phi",0,1);
  }
  else if ((error==0)&&(messagenum==2)&&(dispmode>2)&&(fitmode>3))
  { Textout conout;
    conout.newline();
    if ((layfit[totlayer*4]!=0.0)||(layfit[totlayer*4+2]!=0.0)||
      (layfit[totlayer*4+3]!=0.0))
    { conout.integer(0,8);
      conout.floating(layorig[totlayer*4+2],8.2f,256);
      conout.floating(layorig[totlayer*4+3],8.1f,256);
      conout.floating(layorig[totlayer*4],8.3f,256);
      conout.string(" general vinner tdebye lattice",0,1);
    }
    for (i=0;i<totlayer;++i)
    { if ((layfit[i*4]!=0.0)||(layfit[i*4+2]!=0.0)||
        (layfit[i*4+3]!=0.0))
      { conout.integer(i+1,8);
        conout.floating(layorig[i*4+2],8.4f,256);
        conout.floating(layorig[i*4],8.4f,256);
        conout.floating(layorig[i*4+3],8.4f,256);
        conout.string("   layer spacing length units",0,1);
      }
    }
  }
  if ((error==0)&&(displog>0)&&(messagenum==2)&&(fitmode==2)&&(flogout))
  { flogout->space(8,16); flogout->floating(fitvars[0],8.1f,256);
    flogout->floating(fitvars[1],8.1f,256); flogout->space(11);
    flogout->string("deviation theta and phi",0,1);
    if (error==0) error=flogout->geterror();
  }
  else if ((error==0)&&(displog>0)&&(messagenum==2)&&(fitmode>3)&&
    (flogout))
  { flogout->newline();
    if ((layfit[totlayer*4]!=0.0)||(layfit[totlayer*4+2]!=0.0)||
      (layfit[totlayer*4+3]!=0.0))
    { flogout->integer(0,8);
      flogout->floating(layorig[totlayer*4+2],8.2f,256);
      flogout->floating(layorig[totlayer*4+3],8.1f,256);
      flogout->floating(layorig[totlayer*4],8.3f,256);
      flogout->string(" general vinner tdebye lattice",0,1);
    }
    for (i=0;i<totlayer;++i)
    { if ((layfit[i*4]!=0.0)||(layfit[i*4+2]!=0.0)||
        (layfit[i*4+3]!=0.0))
      { flogout->integer(i+1,8);
        flogout->floating(layorig[i*4+2],8.4f,256);
        flogout->floating(layorig[i*4],8.4f,256);
        flogout->floating(layorig[i*4+3],8.4f,256);
        flogout->string("   layer spacing length units",0,1);
      }
    }
    if (error==0) error=flogout->geterror();
  }

  if ((error==0)&&(messagenum==3)&&(dispmode>3)&&(fitmode==0))
  { Textout conout;
    conout.floating(job,8.2f,256+16); conout.string("% of job ");
    conout.integer(jobnum+1,3); conout.string(" of ");
    conout.integer(jobtotal,3); conout.string(" completed.",0,1);
  }
  else if ((error==0)&&(messagenum==3)&&(dispmode>3))
  { Textout conout;
    conout.floating(job,8.2f,256+16); conout.string("% of try ");
    conout.integer(trynum+1,3); conout.string(", job ");
    conout.integer(jobnum+1,3); conout.string(" of ");
    conout.integer(jobtotal,3); conout.string(" completed.",0,1);
  }
  if ((error==0)&&(displog>0)&&(messagenum==3)&&(fitmode==0)&&(flogout))
  { flogout->floating(job,8.2f,256+16); flogout->string("% of job ");
    flogout->integer(jobnum+1,3); flogout->string(" of ");
    flogout->integer(jobtotal,3); flogout->string(" completed.",0,1);
    if (error==0) error=flogout->geterror();
  }
  else if ((error==0)&&(displog>0)&&(messagenum==3)&&(flogout))
  { flogout->floating(job,8.2f,256+16); flogout->string("% of try ");
    flogout->integer(trynum+1,3); flogout->string(", job ");
    flogout->integer(jobnum+1,3); flogout->string(" of ");
    flogout->integer(jobtotal,3); flogout->string(" completed.",0,1);
    if (error==0) error=flogout->geterror();
  }

  if ((error==0)&&(messagenum==4)&&(dispmode>3))
  { Textout conout;
    conout.string("number",6,16); conout.string("kout",15);
    conout.string("theta",15); conout.string("phi",15);
    conout.string("intensity",15,1);
    conout.integer(pointnum+1,6);
    conout.floating(akout,15.2f,256);
    conout.floating(adtheta,15.2f,256);
    conout.floating(adphi,15.2f,256);
    conout.floating(suminten,15.7f,256+1);
  }
  if ((error==0)&&(messagenum==4)&&(dispmode>3)&&(fitmode<2))
  { Textout conout;
    conout.floating(job,8.2f,256); conout.string("% of job ");
    conout.integer(jobnum+1,3); conout.string(" of ");
    conout.integer(jobtotal,3); conout.string(" completed.",0,1);
  }
  else if ((error==0)&&(messagenum==4)&&(dispmode>3))
  { Textout conout;
    conout.floating(job,8.2f,256); conout.string("% of try ");
    conout.integer(trynum+1,4); conout.string(", job ");
    conout.integer(jobnum+1,3); conout.string(" of ");
    conout.integer(jobtotal,3); conout.string(" completed.",0,1);
  }
  if ((error==0)&&(messagenum==4)&&(displog>0)&&(flogout))
  { flogout->string("number",6,16); flogout->string("kout",15);
    flogout->string("theta",15); flogout->string("phi",15);
    flogout->string("intensity",15,1);
    flogout->integer(pointnum+1,6);
    flogout->floating(akout,15.2f,256);
    flogout->floating(adtheta,15.2f,256);
    flogout->floating(adphi,15.2f,256);
    flogout->floating(suminten,15.7f,256+1);
    if (error==0) error=flogout->geterror();
  }
  if ((error==0)&&(displog>0)&&(messagenum==4)&&(fitmode<2)&&(flogout))
  { flogout->floating(job,8.2f,256); flogout->string("% of job ");
    flogout->integer(jobnum+1,3); flogout->string(" of ");
    flogout->integer(jobtotal,3); flogout->string(" completed.",0,1);
    if (error==0) error=flogout->geterror();
  }
  else if ((error==0)&&(displog>0)&&(messagenum==4)&&(flogout))
  { flogout->floating(job,8.2f,256); flogout->string("% of try ");
    flogout->integer(trynum+1,4); flogout->string(", job ");
    flogout->integer(jobnum+1,3); flogout->string(" of ");
    flogout->integer(jobtotal,3); flogout->string(" completed.",0,1);
    if (error==0) error=flogout->geterror();
  }

  if ((error==0)&&(messagenum==5)&&(dispmode>2)&&(fitmode==0))
  { Textout conout;
    conout.string("job ",0,16); conout.integer(jobnum+1,3);
    conout.string(" of "); conout.integer(jobtotal,3);
    conout.string(" completed.",0,1);
  }
  else if ((error==0)&&(messagenum==5)&&(dispmode>2)&&(fitmode==1))
  { Textout conout;
    conout.string("job ",0,16); conout.integer(jobnum+1,3);
    conout.string(" of "); conout.integer(jobtotal,3);
    conout.string(" completed, factors = ");
    conout.floating(areliable,9.4f,256);
    conout.floating(breliable,9.4f,256+1);
  }
  else if ((error==0)&&(messagenum==5)&&(dispmode>2)&&(afitmath>2))
  { Textout conout;
    conout.string("TRY ",0,16); conout.integer(trynum+1,4);
    conout.string(", job "); conout.integer(jobnum+1,3);
    conout.string(" of "); conout.integer(jobtotal,3);
    conout.string(" completed, factors = ");
    conout.floating(areliable,9.4f,256);
    conout.floating(breliable,9.4f,256+1);
  }
  else if ((error==0)&&(messagenum==5)&&(dispmode>2)&&(afitmath>1))
  { Textout conout;
    conout.string("Try ",0,16); conout.integer(trynum+1,4);
    conout.string(", job "); conout.integer(jobnum+1,3);
    conout.string(" of "); conout.integer(jobtotal,3);
    conout.string(" completed, factors = ");
    conout.floating(areliable,9.4f,256);
    conout.floating(breliable,9.4f,256+1);
  }
  else if ((error==0)&&(messagenum==5)&&(dispmode>2)&&(afitmath==1))
  { Textout conout;
    conout.string("try ",0,16); conout.integer(trynum+1,4);
    conout.string(", job "); conout.integer(jobnum+1,3);
    conout.string(" of "); conout.integer(jobtotal,3);
    conout.string(" completed, factors = ");
    conout.floating(areliable,9.4f,256);
    conout.floating(breliable,9.4f,256+1);
  }
  else if ((error==0)&&(messagenum==5)&&(dispmode>2))
  { Textout conout;
    conout.string("job ",0,16); conout.integer(jobnum+1,3);
    conout.string(" of "); conout.integer(jobtotal,3);
    conout.string(" completed, factors = ");
    conout.floating(areliable,9.4f,256);
    conout.floating(breliable,9.4f,256+1);
  }

  if ((error==0)&&(messagenum==5)&&(dispmode>4)&&(fitmode>1))
    waitenter();

  if ((error==0)&&(displog>0)&&(messagenum==5)&&(fitmode==0)&&(flogout))
  { flogout->string("job ",0,16); flogout->integer(jobnum+1,3);
    flogout->string(" of "); flogout->integer(jobtotal,3);
    flogout->string(" completed.",0,1);
    if (error==0) error=flogout->geterror();
  }
  else if ((error==0)&&(displog>0)&&(messagenum==5)&&(fitmode==1)&&
    (flogout))
  { flogout->string("job ",0,16); flogout->integer(jobnum+1,3);
    flogout->string(" of "); flogout->integer(jobtotal,3);
    flogout->string(" completed, factors = ");
    flogout->floating(areliable,9.4f,256);
    flogout->floating(breliable,9.4f,256+1);
    if (error==0) error=flogout->geterror();
  }
  else if ((error==0)&&(displog>0)&&(messagenum==5)&&(afitmath>2)&&
    (flogout))
  { flogout->string("TRY ",0,16); flogout->integer(trynum+1,4);
    flogout->string(", job "); flogout->integer(jobnum+1,3);
    flogout->string(" of "); flogout->integer(jobtotal,3);
    flogout->string(" completed, factors = ");
    flogout->floating(areliable,9.4f,256);
    flogout->floating(breliable,9.4f,256+1);
    if (error==0) error=flogout->geterror();
  }
  else if ((error==0)&&(displog>0)&&(messagenum==5)&&(afitmath>1)&&
    (flogout))
  { flogout->string("Try ",0,16); flogout->integer(trynum+1,4);
    flogout->string(", job "); flogout->integer(jobnum+1,3);
    flogout->string(" of "); flogout->integer(jobtotal,3);
    flogout->string(" completed, factors = ");
    flogout->floating(areliable,9.4f,256);
    flogout->floating(breliable,9.4f,256+1);
    if (error==0) error=flogout->geterror();
  }
  else if ((error==0)&&(displog>0)&&(messagenum==5)&&(afitmath>1)&&
    (flogout))
  { flogout->string("try ",0,16); flogout->integer(trynum+1,4);
    flogout->string(", job "); flogout->integer(jobnum+1,3);
    flogout->string(" of "); flogout->integer(jobtotal,3);
    flogout->string(" completed, factors = ");
    flogout->floating(areliable,9.4f,256);
    flogout->floating(breliable,9.4f,256+1);
    if (error==0) error=flogout->geterror();
  }

  return(error);
} //end of Mscdrun::dispintensity

//make the coordinates of each atom
int Mscdrun::makeatoms(float disturb)
{ int i,j,k,n;
  float xa,ya;
  const float radian=(float)(3.14159265/180.0);

  k=0;
  for (n=0;(error==0)&&(n<nfit);++n)
  { ++k;
    j=(int)fitvars[nfit*2+n]; i=(int)fitvars[nfit*3+n];
    layorig[i*4+j]+=(float)(k*disturb*0.1/lattice);
  }

  for (i=0;(error==0)&&(i<totlayer);++i)
  { xa=laycell[i*4]*lattice*(float)cos(laycell[i*4+1]*radian);
    ya=laycell[i*4]*lattice*(float)sin(laycell[i*4+1]*radian);
    laxcell[i*4]=xa*layorig[i*4+3];
    laxcell[i*4+1]=ya*layorig[i*4+3];
    xa=laycell[i*4+2]*lattice*(float)cos(laycell[i*4+3]*radian);
    ya=laycell[i*4+2]*lattice*(float)sin(laycell[i*4+3]*radian);
    laxcell[i*4+2]=xa*layorig[i*4+3];
    laxcell[i*4+3]=ya*layorig[i*4+3];
    xa=layorig[i*4]*lattice*(float)cos(layorig[i*4+1]*radian);
    ya=layorig[i*4]*lattice*(float)sin(layorig[i*4+1]*radian);
    laxorig[i*4]=xa;
    laxorig[i*4+1]=ya;
    if (i==0) laxorig[i*4+2]=0.0f-layorig[i*4+2]*lattice;
    else laxorig[i*4+2]=laxorig[(i-1)*4+2]-layorig[i*4+2]*lattice;
  }

  for (n=0;(error==0)&&(n<natoms);++n)
  { i=(int)patom[n*12+3]; j=(int)patom[n*12+4]; k=(int)patom[n*12+5];
    patom[n*12]=i*laxcell[k*4]+j*laxcell[k*4+2]+laxorig[k*4];
    patom[n*12+1]=i*laxcell[k*4+1]+j*laxcell[k*4+3]+laxorig[k*4+1];
    patom[n*12+2]=laxorig[k*4+2];
  }
  return(error);
} //Mscdrun::makeatoms

//make triple vertex event parameters
int Mscdrun::maketripar()
{ int ia,ib,ic,j;
  float xa,ya,za,xb,yb,zb,lenga,lengb,vlena,vlenb,vlenc,cosab;

  vlenc=kmin/100.0f;
  for (j=0;(error==0)&&(j<ntrieven);++j)
  { ia=(int)tevenpar[j*10+7]; ib=(int)tevenpar[j*10+8];
    ic=(int)tevenpar[j*10+9];
    xa=patom[ib*12]-patom[ia*12]; ya=patom[ib*12+1]-patom[ia*12+1];
    za=patom[ib*12+2]-patom[ia*12+2];
    lenga=(float)sqrt(xa*xa+ya*ya+za*za);
    xb=patom[ic*12]-patom[ib*12]; yb=patom[ic*12+1]-patom[ib*12+1];
    zb=patom[ic*12+2]-patom[ib*12+2];
    lengb=(float)sqrt(xb*xb+yb*yb+zb*zb);
    if ((lenga<1.0e-3)||(lengb<1.0e-3)) error=605;
    else
    { cosab=(xb*xa+yb*ya+zb*za)/lenga/lengb;
      vlena=1.0f/lenga; vlenb=1.0f/lengb;
      if (cosab>1.0) cosab=1.0f; else if (cosab<-1.0) cosab=-1.0f;
      if (ntrieven>1000)
      { vlena=round(vlena,vlenc); vlenb=round(vlenb,vlenc);
        lenga=round(lenga,0.001f); cosab=round(cosab,0.001f);
      }
      tevenpar[j*10]=lenga; tevenpar[j*10+1]=vlena;
      tevenpar[j*10+2]=vlenb; tevenpar[j*10+3]=cosab;
    }
  }
  return(error);
} //end of Mscdrun::maketripar

//make the double vertex event parameters
int Mscdrun::makedblpar()
{ int ia,ib,j;
  float xa,ya,za,lenga;

  for (j=0;(error==0)&&(j<ndbleven);++j)
  { ia=(int)devenpar[j*7+5]; ib=(int)devenpar[j*7+6];
    xa=patom[ib*12]-patom[ia*12]; ya=patom[ib*12+1]-patom[ia*12+1];
    za=patom[ib*12+2]-patom[ia*12+2];
    lenga=(float)sqrt(xa*xa+ya*ya+za*za);
    devenpar[j*7]=xa; devenpar[j*7+1]=ya; devenpar[j*7+2]=za;
    devenpar[j*7+3]=lenga;
  }
  return(error);
} //end of Mscdrun::makedblpar

float Mscdrun::kinside(float akout,float avinner)
{ float akin;
  const float kcoeff=0.512332f;
  if ((error==0)&&(avinner<0.1)) akin=akout;
  else akin=(float)sqrt(akout*akout+kcoeff*kcoeff*avinner);
  return(akin);
} //end of Mscdrun::kinside

float Mscdrun::thetainside(int accepnum,float akin,float akout,
  float thetaout,float phiout,float altheta,float alphi,
  float *xdetec,float *polaron)
{ float xa,xb,thetain,bdtheta,bdphi,cdtheta,cdphi;
  const float radian=(float)(3.14159265/180.0);
  Cartesia ydeteca,ydetecb;

  if (accepnum==0) bdtheta=bdphi=0.0f;
  else
  { bdtheta=accepang; bdphi=(accepnum-1.0f)*90.0f;
  }
  ydeteca.loadthetaphi(bdtheta,bdphi);
  ydetecb=ydeteca.euler(phiout,thetaout,0.0f);
  cdtheta=ydetecb.theta(); cdphi=ydetecb.phi();
  if ((error==0)&&((akin<fabs(akout)+1.0e-3)||(cdtheta>89.0)))
  { xa=(float)sin(cdtheta*radian); thetain=cdtheta;
  }
  else if (error==0)
  { xa=akout*(float)sin(cdtheta*radian)/akin;
    thetain=(float)asin(xa)/radian;
  }
  else xa=thetain=0.0f;
  if (error==0)
  { xdetec[0]=xa*(float)cos(cdphi*radian);
    xdetec[1]=xa*(float)sin(cdphi*radian);
    xdetec[2]=(float)cos(thetain*radian);
    xb=(float)sin(altheta*radian);
    polaron[0]=xb*(float)cos(alphi*radian);
    polaron[1]=xb*(float)sin(alphi*radian);
    polaron[2]=(float)cos(altheta*radian);
    polaron[3]=altheta; polaron[4]=alphi;
    polaron[5]=-(float)sin(alphi*radian);
    polaron[6]=(float)cos(alphi*radian);
    polaron[7]=0.0f;
    polaron[8]=90.0f; polaron[9]=alphi+90.0f;
  }

  return(thetain);
} //end of Mscdrun::thetainside

