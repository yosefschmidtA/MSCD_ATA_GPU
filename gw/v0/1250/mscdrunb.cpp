#include <iostream>
using namespace std;
#include <fstream>
#include <iomanip>
#include <cmath>

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

int Mscdrun::paradisp()
{ int i,j;
  float akout,adtheta,adphi,altheta,alphi,unknown;

  Textout conout;
  if ((error==0)&&(dispmode>4))
  { conout.newline();
    if (katoms>0)
    { conout.string("  psfile: "); conout.string(psfile,20);
    }
    if (katoms>1)
    { conout.space(10); conout.string(psfile+100,20);
    }
    conout.newline();
    for (i=2;i<katoms;++i)
    { conout.space(10);
      if (i<katoms-1) conout.string((psfile+i*100),20);
      else conout.string((psfile+i*100),20,1);
    }
    if (linitial>0)
    { conout.string("  rmfile: "); conout.string(rmfile,20);
    }
    conout.string("  pefile: "); conout.string(foutname,20,2);
    conout.integer(scanmode,8); conout.integer(dispmode+displog*10,8);
    conout.floating(ftolerance,8.3f,256); conout.space(11);
    conout.string("scanmode dispmode tolerance",0,2);
    conout.integer(linitial,8); conout.integer(lnum,8);
    conout.integer(msorder,8); conout.integer(raorder,8);
    conout.string("   linitial lnum msorder raorder",0,1);
    conout.integer(totlayer,8); conout.integer(finals,8);
    conout.integer(fitmath,8); conout.integer(trymax,8);
    conout.string("   totlayer finals fitmath trymax",0,1);
    conout.floating(kmin,8.2f,256); conout.floating(kmax,8.2f,256);
    conout.floating(kstep,8.2f,256); conout.space(11);
    conout.string("kmin kmax kstep",0,1);
    conout.floating(dtmin,8.1f,256); conout.floating(dtmax,8.1f,256);
    conout.floating(dtstep,8.2f,256); conout.space(11);
    conout.string("dthetamin dthetamax dthetastep",0,1);
    conout.floating(dpmin,8.1f,256); conout.floating(dpmax,8.1f,256);
    conout.floating(dpstep,8.2f,256); conout.space(11);
    conout.string("dphimin dphimax dphistep",0,1);
    conout.floating(ltheta,8.1f,256); conout.floating(lphi,8.1f,256);
    conout.integer(beampol,8); conout.space(11);
    conout.string("ltheta lphi beampol",0,1);
    conout.floating(mtheta,8.1f,256); conout.floating(mphi,8.1f,256);
    conout.floating(accepang,8.1f,256); conout.space(11);
    conout.string("mtheta mphi accepang",0,1);
    conout.floating(radius,8.2f,256); conout.floating(depth,8.2f,256);
    conout.floating(lattice,8.2f,256); conout.space(11);
    conout.string("radius depth lattice",0,1);
    conout.floating(valence,8.2f,256);
    conout.floating(bandgap,8.2f,256);
    conout.floating(density,8.2f,256);
    conout.floating(mweight,8.2f,256);
    conout.string("   valence bandgap density mweight",0,1);
    for (i=0;i<katoms;++i)
    { conout.floating(aweight[i],8.2f,256);
    }
    conout.space((4-katoms)*8+3);
    conout.string("effective weight for kind 1-");
    conout.integer(katoms,0,1);
    for (i=0;i<katoms;++i)
    { conout.floating(magamp[i],8.2f,256);
    }
    conout.space((4-katoms)*8+3);
    conout.string("magnetization amplitude for kind 1-");
    conout.integer(katoms,0,1);
    conout.floating(vinner,8.2f,256);
    conout.floating(tdebye,8.1f,256);
    conout.floating(tsample,8.1f,256);
    if ((pathcut<1.0e-10)||(pathcut>1.0e-4))
      conout.floating(pathcut,10.4f,256);
    else conout.floating(pathcut,10.1f,512);
    conout.string(" vinner tdebye tsample pathcut",0,1);
    conout.floating(layfit[totlayer*4+2],8.3f,256);
    conout.floating(layfit[totlayer*4+3],8.3f,256);
    conout.floating(layfit[totlayer*4],8.3f,256);
    conout.space(11);
    conout.string("fitting try for vinner, tdebye and lattice",0,1);
    waitenter();
  }
  if ((error==0)&&(dispmode>5))
  { for (i=0;i<totlayer;++i)
    { if ((i>0)&&(lakatom[i*4+3]<0)) conout.integer(-i-1,8);
      else conout.integer(i+1,8);
      for (j=0;j<3;++j) conout.integer(lakatom[i*4+j],8);
      conout.string("   layer katom emiter atoms",0,1);
      for (j=0;j<4;++j) conout.integer(layatom[i*4+j],8);
      conout.string("   latoma(x1,x2) latomb(y1,y2)",0,1);
      conout.floating(laycell[i*4],8.3f,256);
      conout.floating(laycell[i*4+1],8.1f,256);
      conout.floating(laycell[i*4+2],8.3f,256);
      conout.floating(laycell[i*4+3],8.1f,256);
      conout.string("   unita(len,ang) unitb(len,ang)",0,1);
      conout.floating(layorig[i*4],8.3f,256);
      conout.floating(layorig[i*4+1],8.3f,256);
      conout.floating(layorig[i*4+2],8.3f,256);
      conout.space(11);
      conout.string("origin(len ang spacing)",0,1);
      conout.floating(layfit[i*4+2],8.3f,256);
      conout.floating(layfit[i*4],8.3f,256);
      conout.floating(layfit[i*4+3],8.3f,256);
      conout.space(11);
      conout.string("fitting try for spacing, length and units",0,1);
      if ((i%3==2)||(i==totlayer-1)) waitenter();
      else conout.newline();
    }
  }
  if ((error==0)&&(dispmode>6))
  { for (i=0;i<natoms;++i)
    { conout.integer(i+1,6);
      for (j=6;j<7;++j) conout.integer((int)patom[i*12+j],6);
      conout.floating(patom[i*12+7],6.2f,256);
      for (j=0;j<3;++j)
      { conout.floating(patom[i*12+j],9.4f,256);
      }
      conout.string("     no katom emiter x y z",0,1);
      if (i==natoms-1)
      { conout.floating(nearest,9.4f,256+16);
        conout.floating(biggest,9.4f,256);
        conout.string("     nearest and farrest bonding length",0,1);
      }
      if ((i%15==14)||(i==natoms-1)) waitenter();
    }
  }
  if ((error==0)&&(dispmode>7))
  { conout.newline();
    for (i=0;i<npoint;++i)
    { error=pdintensity->getpoint(i,&akout,&adtheta,&adphi,&altheta,
        &alphi,&unknown,&unknown,&unknown,&unknown);
      conout.integer(i+1,6);
      conout.floating(akout,8.2f,256);
      conout.floating(adtheta,8.2f,256);
      conout.floating(adphi,8.2f,256);
      conout.floating(altheta,8.2f,256);
      conout.floating(alphi,8.2f,256);
      conout.string("   no k dtheta dphi ltheta lphi",0,1);
      if (((i%15)==14)||(i==npoint-1)) waitenter();
    }
  }

  if ((error==0)&&(dispmode>2))
  { conout.integer(linitial,8,16); conout.integer(lnum,8);
    conout.integer(msorder,8); conout.integer(raorder,8);
    conout.string("   linitial lnum msorder raorder",0,1);
    conout.integer(natoms,8); conout.integer(eatoms,8);
    conout.integer(katoms,8); conout.integer(totlayer,8);
    conout.string("   natoms emiters katoms nlayer",0,1);
    conout.integer(scanmode,8); conout.integer(nfit,8);
    conout.integer(nfollow,8); conout.integer(npoint,8);
    conout.string("   scanmode nfit nfollow npoint",0,2);
  }
  if ((error==0)&&(dispmode>4)) waitenter();

  return(error);
} //end of Mscdrun::paradisp

int Mscdrun::paradisplog()
{ int i,j;
  float akout,adtheta,adphi,altheta,alphi,unknown;

  if ((error==0)&&(displog>0)&&(flogout))
  { flogout->newline();
    if (katoms>0)
    { flogout->string("  psfile: "); flogout->string(psfile,20);
    }
    if (katoms>1)
    { flogout->space(10); flogout->string(psfile+100,20);
    }
    flogout->newline();
    for (i=2;i<katoms;++i)
    { flogout->space(10);
      if (i<katoms-1) flogout->string((psfile+i*100),20);
      else flogout->string((psfile+i*100),20,1);
    }
    if (linitial>0)
    { flogout->string("  rmfile: "); flogout->string(rmfile,20);
    }
    flogout->string("  pefile: "); flogout->string(foutname,20,2);
    flogout->integer(scanmode,8);
    flogout->integer(dispmode+displog*10,8);
    flogout->floating(ftolerance,8.3f,256); flogout->space(11);
    flogout->string("scanmode dispmode tolerance",0,2);
    flogout->integer(linitial,8); flogout->integer(lnum,8);
    flogout->integer(msorder,8); flogout->integer(raorder,8);
    flogout->string("   linitial lnum msorder raorder",0,1);
    flogout->integer(totlayer,8); flogout->integer(finals,8);
    flogout->integer(fitmath,8); flogout->integer(trymax,8);
    flogout->string("   totlayer finals fitmath trymax",0,1);
    flogout->floating(kmin,8.2f,256);
    flogout->floating(kmax,8.2f,256);
    flogout->floating(kstep,8.2f,256); flogout->space(11);
    flogout->string("kmin kmax kstep",0,1);
    flogout->floating(dtmin,8.1f,256);
    flogout->floating(dtmax,8.1f,256);
    flogout->floating(dtstep,8.2f,256); flogout->space(11);
    flogout->string("dthetamin dthetamax dthetastep",0,1);
    flogout->floating(dpmin,8.1f,256);
    flogout->floating(dpmax,8.1f,256);
    flogout->floating(dpstep,8.2f,256); flogout->space(11);
    flogout->string("dphimin dphimax dphistep",0,1);
    flogout->floating(ltheta,8.1f,256);
    flogout->floating(lphi,8.1f,256);
    flogout->integer(beampol,8); flogout->space(11);
    flogout->string("ltheta lphi beampol",0,1);
    flogout->floating(mtheta,8.1f,256);
    flogout->floating(mphi,8.1f,256);
    flogout->floating(accepang,8.1f,256); flogout->space(11);
    flogout->string("mtheta mphi accepang",0,1);
    flogout->floating(radius,8.2f,256);
    flogout->floating(depth,8.2f,256);
    flogout->floating(lattice,8.2f,256); flogout->space(11);
    flogout->string("radius depth lattice",0,1);
    flogout->floating(valence,8.2f,256);
    flogout->floating(bandgap,8.2f,256);
    flogout->floating(density,8.2f,256);
    flogout->floating(mweight,8.2f,256);
    flogout->string("   valence bandgap density mweight",0,1);
    for (i=0;i<katoms;++i)
    { flogout->floating(aweight[i],8.2f,256);
    }
    flogout->space((4-katoms)*8+3);
    flogout->string("effective weight for kind 1-");
    flogout->integer(katoms,0,1);
    for (i=0;i<katoms;++i)
    { flogout->floating(magamp[i],8.2f,256);
    }
    flogout->space((4-katoms)*8+3);
    flogout->string("magnetization amplitude for kind 1-");
    flogout->integer(katoms,0,1);
    flogout->floating(vinner,8.2f,256);
    flogout->floating(tdebye,8.1f,256);
    flogout->floating(tsample,8.1f,256);
    if ((pathcut<1.0e-10)||(pathcut>1.0e-4))
      flogout->floating(pathcut,10.4f,256);
    else flogout->floating(pathcut,10.1f,512);
    flogout->string(" vinner tdebye tsample pathcut",0,1);
    flogout->floating(layfit[totlayer*4+2],8.3f,256);
    flogout->floating(layfit[totlayer*4+3],8.3f,256);
    flogout->floating(layfit[totlayer*4],8.3f,256);
    flogout->space(11);
    flogout->string("fitting try for vinner, tdebye and lattice",0,2);
  }

  if ((error==0)&&(displog>0)&&(flogout))
  { for (i=0;i<totlayer;++i)
    { if ((i>0)&&(lakatom[i*4+3]<0)) flogout->integer(-i-1,8);
      else flogout->integer(i+1,8);
      for (j=0;j<3;++j) flogout->integer(lakatom[i*4+j],8);
      flogout->string("   layer katom emiter atoms",0,1);
      for (j=0;j<4;++j) flogout->integer(layatom[i*4+j],8);
      flogout->string("   latoma(x1,x2) latomb(y1,y2)",0,1);
      flogout->floating(laycell[i*4],8.3f,256);
      flogout->floating(laycell[i*4+1],8.1f,256);
      flogout->floating(laycell[i*4+2],8.3f,256);
      flogout->floating(laycell[i*4+3],8.1f,256);
      flogout->string("   unita(len,ang) unitb(len,ang)",0,1);
      flogout->floating(layorig[i*4],8.3f,256);
      flogout->floating(layorig[i*4+1],8.3f,256);
      flogout->floating(layorig[i*4+2],8.3f,256);
      flogout->space(11);
      flogout->string("origin(len ang spacing)",0,1);
      flogout->floating(layfit[i*4+2],8.3f,256);
      flogout->floating(layfit[i*4],8.3f,256);
      flogout->floating(layfit[i*4+3],8.3f,256);
      flogout->space(11);
      flogout->string("fitting try for spacing, length and units",0,2);
    }
  }

  if ((error==0)&&(displog>0)&&(flogout))
  { for (i=0;i<natoms;++i)
    { flogout->integer(i+1,6);
      for (j=6;j<7;++j) flogout->integer((int)patom[i*12+j],6);
      flogout->floating(patom[i*12+7],6.2f,256);
      for (j=0;j<3;++j)
      { flogout->floating(patom[i*12+j],9.4f,256);
      }
      flogout->string("     no katom emiter x y z",0,1);
      if (i==natoms-1)
      { flogout->floating(nearest,9.4f,256+16);
        flogout->floating(biggest,9.4f,256);
        flogout->string("     nearest and farrest bonding length",0,1);
      }
    }
  }

  if ((error==0)&&(displog>0)&&(flogout))
  { flogout->newline();
    for (i=0;i<npoint;++i)
    { error=pdintensity->getpoint(i,&akout,&adtheta,&adphi,&altheta,
        &alphi,&unknown,&unknown,&unknown,&unknown);
      flogout->integer(i+1,6);
      flogout->floating(akout,8.2f,256);
      flogout->floating(adtheta,8.2f,256);
      flogout->floating(adphi,8.2f,256);
      flogout->floating(altheta,8.2f,256);
      flogout->floating(alphi,8.2f,256);
      flogout->string("   no k dtheta dphi ltheta lphi",0,1);
    }
  }

  if ((error==0)&&(displog>0)&&(flogout))
  { flogout->integer(linitial,8,16); flogout->integer(lnum,8);
    flogout->integer(msorder,8); flogout->integer(raorder,8);
    flogout->string("   linitial lnum msorder raorder",0,1);
    flogout->integer(natoms,8); flogout->integer(eatoms,8);
    flogout->integer(katoms,8); flogout->integer(totlayer,8);
    flogout->string("   natoms emiters katoms nlayer",0,1);
    flogout->integer(scanmode,8); flogout->integer(nfit,8);
    flogout->integer(nfollow,8); flogout->integer(npoint,8);
    flogout->string("   scanmode nfit nfollow npoint",0,2);
  }
  if ((error==0)&&(displog>0)&&(flogout)) error=flogout->geterror();

  return(error);
} //end of Mscdrun::paradisplog

/*
----------------------------------------------------------------------
tevenpar   0        1            2             3       4
          lengtha  inverse_lena  inverse_lengb cosbeta atom_kind

          5        6       7     8     9
          eledim   memadd  ia    ib   ic
----------------------------------------------------------------------
*/
int Mscdrun::symtrivert()
{ int ia,ib,ic,j,k,m,n,p,q,tscatter,analyze;
  float xa,ya,za,xb,yb,zb,lena,lenb,vlenb,vlenc,awidth,astep,
    cosbeta,akind,search;
  int *tempadd;
  float *tempar;

  analyze=0;
  Textout conout;
  if (((error==0)&&(dispmode>1)&&(msorder>2)&&(natoms>20))||
    ((error==0)&&(dispmode>1)&&(msorder==2)&&(natoms>100))||
    ((error==0)&&(dispmode>4)))
    conout.string("Analyzing, please wait ...",0,1);

  tempadd=NULL; tempar=NULL; vlenc=kmin/100.0f;
  nsymm=natoms*natoms/20; if (nsymm<2) nsymm=2;
  if (natoms<15) m=10;
  else if (natoms<200) m=5-natoms/50;
  else m=1;
  if (msorder>2) mscatter=m*natoms*natoms*natoms;
  else if (msorder==2) mscatter=m*eatoms*natoms*natoms;
  else mscatter=ntrieven=nscorse=nsymm=0;

  ia=natoms+10; tscatter=0; awidth=astep=0.0f;
  while ((error==0)&&(nsymm>=2)&&(ia>natoms))
  { ++analyze;
    tscatter=mscatter/nsymm/3;
    awidth=round(1.0f/nearest,vlenc)+2.01f;
    astep=(nsymm-1.0f)/(katoms*awidth);
    if ((error==0)&&(mscatter>0))
    { tempadd=new int [nsymm]; tempar=new float [mscatter];
      if ((!tempadd)||(!tempar)) error=102;
      else for (j=0;j<nsymm;++j) tempadd[j]=0;
    }
    for (ia=0;(error==0)&&(ia<natoms);++ia)
    { if ((msorder==2)&&(patom[ia*12+7]==0)) continue;
      for (ib=0;ib<natoms;++ib)
      { if (ib==ia) continue;
        xa=patom[ib*12]-patom[ia*12]; ya=patom[ib*12+1]-patom[ia*12+1];
        za=patom[ib*12+2]-patom[ia*12+2];
        lena=(float)sqrt(xa*xa+ya*ya+za*za);
        akind=patom[ib*12+6]; lena=round(lena,0.001f);
        for (ic=0;ic<natoms;++ic)
        { if (ic==ib) continue;
          xb=patom[ic*12]-patom[ib*12];
          yb=patom[ic*12+1]-patom[ib*12+1];
          zb=patom[ic*12+2]-patom[ib*12+2];
          lenb=(float)sqrt(xb*xb+yb*yb+zb*zb);
          cosbeta=(xb*xa+yb*ya+zb*za)/(lena*lenb);

          vlenb=round(1.0f/lenb,vlenc);
          cosbeta=round(cosbeta,0.001f);
          search=(float)((akind-1.0)*awidth+vlenb+cosbeta+1.001);
          search=round(search,0.001f);
          m=(int)(search*astep)+1;
          if (m<1) m=1; else if (m>nsymm-1) m=nsymm-1;
          n=m*tscatter;
          for (j=n;j<n+tempadd[m];++j)
          { if ((lena==tempar[j*3])&&(vlenb==tempar[j*3+1])&&
              (cosbeta==tempar[j*3+2])) j=n+tempadd[m]+10;
          }
          if ((j==n+tempadd[m])&&(j*3+2<mscatter))
          { tempar[j*3]=lena; tempar[j*3+1]=vlenb;
            tempar[j*3+2]=cosbeta;
            ++tempadd[m];
          }
          if ((tempadd[m]>tscatter)&&(nsymm>=2))
          { if (((dispmode>1)&&(msorder>2)&&(natoms>20))||
              ((dispmode>1)&&(msorder==2)&&(natoms>100))||
              (dispmode>4))
              conout.string("Reanalyzing, please wait ...",0,1);
            if ((nsymm<10)&&(mscatter<=5*natoms*natoms*natoms))
              mscatter=mscatter*3/2;
            nsymm/=2; if (nsymm<2) error=621;
            ia=ib=ic=natoms+10;
          }
        }
      }
    }
  }

  k=0; if (tempadd) tempadd[0]=k;
  for (m=1;(error==0)&&(m<nsymm);++m)
  { n=m*tscatter;
    for (p=n;p<n+tempadd[m];++p)
    { for (q=p+1;q<n+tempadd[m];++q)
      { if (tempar[q*3+2]<tempar[p*3+2])
        { for (j=0;j<3;++j)
          { xa=tempar[p*3+j]; tempar[p*3+j]=tempar[q*3+j];
            tempar[q*3+j]=xa;
          }
        }
      }
    }
    for (p=n;p<n+tempadd[m];++p)
    { for (q=p+1;q<n+tempadd[m];++q)
      { if (tempar[q*3+2]!=tempar[p*3+2]) q=n+tempadd[m];
        else if (tempar[q*3+1]<tempar[p*3+1])
        { for (j=0;j<3;++j)
          { xa=tempar[p*3+j]; tempar[p*3+j]=tempar[q*3+j];
            tempar[q*3+j]=xa;
          }
        }
      }
    }
    for (p=n;p<n+tempadd[m];++p)
    { for (q=p+1;q<n+tempadd[m];++q)
      { if ((tempar[q*3+2]!=tempar[p*3+2])||
          (tempar[q*3+1]!=tempar[p*3+1])) q=n+tempadd[m];
        else if (tempar[q*3]<tempar[p*3])
        { for (j=0;j<3;++j)
          { xa=tempar[p*3+j]; tempar[p*3+j]=tempar[q*3+j];
            tempar[q*3+j]=xa;
          }
        }
      }
    }
    for (p=n;p<n+tempadd[m];++p)
    { if (tempar[p*3]==0.0) continue;
      for (q=p;q<n+tempadd[m];++q)
      { if (tempar[q*3]==0.0) continue;
        else if ((k==0)||(q==p)||((tempar[(k-1)*3+1]==tempar[q*3+1])&&
          (tempar[(k-1)*3+2]==tempar[q*3+2])))
        { tempar[k*3]=tempar[q*3]; tempar[k*3+1]=tempar[q*3+1];
          tempar[k*3+2]=tempar[q*3+2]; tempar[q*3]=0.0f;
          ++k;
        }
      }
    }
    tempadd[m]=k;
  }
  ntrieven=k;

  if (error==0)
  { if (tevenpar) delete [] tevenpar;
    tevenpar=new float [ntrieven*10];
    if (!tevenpar) error=102;
  }
  for (j=0;(error==0)&&(j<ntrieven);++j)
  { tevenpar[j*10]=tempar[j*3];
    tevenpar[j*10+1]=round(1.0f/tempar[j*3],vlenc);
    tevenpar[j*10+2]=tempar[j*3+1];
    tevenpar[j*10+3]=tempar[j*3+2];
    for (k=4;k<10;++k) tevenpar[j*10+k]=0.0f;
  }

  if (tempar) delete [] tempar;
  if (error==0)
  { if (tevenadd) delete [] tevenadd;
    tevenadd=new int [natoms*natoms*natoms];
    if (!tevenadd) error=102;
  }

  for (ia=0;(error==0)&&(ntrieven>0)&&(ia<natoms);++ia)
  { if ((msorder==2)&&(patom[ia*12+7]==0)) continue;
    for (ib=0;(error==0)&&(ib<natoms);++ib)
    { if (ib==ia) continue;
      xa=patom[ib*12]-patom[ia*12]; ya=patom[ib*12+1]-patom[ia*12+1];
      za=patom[ib*12+2]-patom[ia*12+2];
      lena=(float)sqrt(xa*xa+ya*ya+za*za);
      akind=patom[ib*12+6]; lena=round(lena,0.001f);
      for (ic=0;(error==0)&&(ic<natoms);++ic)
      { if (ic==ib) continue;
        xb=patom[ic*12]-patom[ib*12]; yb=patom[ic*12+1]-patom[ib*12+1];
        zb=patom[ic*12+2]-patom[ib*12+2];
        lenb=(float)sqrt(xb*xb+yb*yb+zb*zb);
        cosbeta=(xb*xa+yb*ya+zb*za)/(lena*lenb);

        vlenb=round(1.0f/lenb,vlenc);
        cosbeta=round(cosbeta,0.001f);
        search=(float)((akind-1.0)*awidth+vlenb+cosbeta+1.001);
        search=round(search,0.001f);
        m=(int)(search*astep);
        if (m<0) m=0; else if (m>nsymm-2) m=nsymm-2;
        for (j=tempadd[m];j<tempadd[m+1];++j)
        { if ((lena==tevenpar[j*10])&&(vlenb==tevenpar[j*10+2])&&
            (cosbeta==tevenpar[j*10+3]))
          { tevenadd[ia*natoms*natoms+ib*natoms+ic]=j;
            if (tevenpar[j*10+4]<0.5)
            { tevenpar[j*10+4]=(float)akind;
              tevenpar[j*10+7]=(float)ia;
              tevenpar[j*10+8]=(float)ib; tevenpar[j*10+9]=(float)ic;
            }
            j=tempadd[m+1]+10;
          }
        }
        if (j==tempadd[m+1]) error=621;
      }
    }
  }

  p=q=0;
  for (m=0;(error==0)&&(m<nsymm-1);++m)
  { n=tempadd[m+1]-tempadd[m];
    if (p<n)
    { p=n; ia=m;
    }
  }
  for (m=0;(error==0)&&(m<nsymm-1);++m)
  { n=tempadd[m+1]-tempadd[m];
    if ((q<n)&&(m!=ia))
    { q=n; ib=m;
    }
  }
  if ((error==0)&&(ntrieven>0))
  { xa=(float)(p*100.0/tscatter); xb=(float)(q*100.0/tscatter);
    k=0;
    for (j=0;j<ntrieven;++j)
    { tevenpar[j*10+3]=confine(tevenpar[j*10+3],-1.0f,1.0f);
      if ((j>0)&&(tevenpar[j*10+1]==tevenpar[(j-1)*10+1])&&
        (tevenpar[j*10+2]==tevenpar[(j-1)*10+2])&&
        (tevenpar[j*10+3]==tevenpar[(j-1)*10+3])&&
        (tevenpar[j*10+4]==tevenpar[(j-1)*10+4]))
        ++k;
    }
    nscorse=ntrieven-k;
    if ((msorder>2)&&(natoms>1))
      ya=(float)(ntrieven*100.0/natoms/(natoms-1.0)/(natoms-1.0));
    else if ((msorder==2)&&(eatoms>0)&&(natoms>1))
      ya=(float)(ntrieven*100.0/eatoms/(natoms-1.0)/(natoms-1.0));
    else ya=0.0f;
    yb=(float)(nscorse*100.0/ntrieven);
    if (natoms>1) za=(float)(mscatter/natoms/natoms/natoms);
    else za=0.0f;
    if ((error==0)&&(dispmode>4))
    { conout.integer(nsymm,8); conout.floating(za,8.1f,256);
      conout.floating(xa,8.2f,256); conout.floating(xb,8.2f,256);
      conout.string("   nsymm mscatter distribute percentages",0,1);
      conout.integer(ntrieven,8); conout.integer(nscorse,8);
      conout.floating(ya,8.1f,256); conout.floating(yb,8.1f,256);
      conout.string("   ntrieven nscorse their percents",0,1);
    }
  }
  else xa=xb=ya=yb=za=0.0f;
  if ((error==0)&&(ntrieven>0)&&(displog>0)&&(flogout))
  { flogout->string("Analyzed symmetries for ");
    flogout->integer(analyze); flogout->string(" times",0,1);
    flogout->integer(nsymm,8); flogout->floating(za,8.1f,256);
    flogout->floating(xa,8.2f,256); flogout->floating(xb,8.2f,256);
    flogout->string("   nsymm mscatter distribute percentages",0,1);
    flogout->integer(ntrieven,8); flogout->integer(nscorse,8);
    flogout->floating(ya,8.1f,256); flogout->floating(yb,8.1f,256);
    flogout->string("   ntrieven nscorse their percents",0,1);
    if (error==0) error=flogout->geterror();
  }

  if (tempadd) delete [] tempadd;

  return(error);
} //end of Mscdrun::symtrivert

/*
----------------------------------------------------------------------
devenpar   0   1   2   3        4          5   6
          xa  ya  za  lengtha  atom_kind  ia  ib
----------------------------------------------------------------------
*/
int Mscdrun::symdblvert()
{ int ia,ib,j,n;
  float xa,ya,za,lena,akind;
  float *tempar;

  tempar=NULL;
  if ((error==0)&&(msorder>1)) ndbleven=natoms*natoms;
  else if ((error==0)&&(msorder==1)) ndbleven=eatoms*natoms;
  else ndbleven=0;
  if ((error==0)&&(ndbleven>0))
  { if (devenadd) delete [] devenadd;
    devenadd=new int [natoms*natoms]; tempar=new float [ndbleven*7];
    if ((!devenadd)&&(!tempar)) error=102;
  }

  n=0;
  for (ia=0;(error==0)&&(ia<natoms);++ia)
  { if ((msorder<1)||((msorder==1)&&(patom[ia*12+7]==0))) continue;
    for (ib=0;(error==0)&&(ib<natoms);++ib)
    { if (ib==ia) continue;
      xa=patom[ib*12]-patom[ia*12]; ya=patom[ib*12+1]-patom[ia*12+1];
      za=patom[ib*12+2]-patom[ia*12+2];
      lena=(float)sqrt(xa*xa+ya*ya+za*za);
      akind=patom[ib*12+6];
      xa=round(xa,1.0e-4f); ya=round(ya,1.0e-4f);
      za=round(za,1.0e-4f); lena=round(lena,1.0e-4f);
      for (j=0;(error==0)&&(j<n);++j)
      { if ((xa==tempar[j*7])&&(ya==tempar[j*7+1])&&
          (za==tempar[j*7+2])&&(akind==tempar[j*7+4]))
        { devenadd[ia*natoms+ib]=j;
          j=n+10;
        }
      }
      if (j==n)
      { if (n<ndbleven-1)
        { devenadd[ia*natoms+ib]=n;
          tempar[n*7]=xa; tempar[n*7+1]=ya; tempar[n*7+2]=za;
          tempar[n*7+3]=lena; tempar[n*7+4]=(float)akind;
          tempar[n*7+5]=(float)ia; tempar[n*7+6]=(float)ib;
          ++n;
        }
        else error=621;
      }
    }
  }
  ndbleven=n;
  if ((error==0)&&(ndbleven>0))
  { if (devenpar) delete [] devenpar;
    devenpar=new float [ndbleven*7];
    if (!devenpar) error=102;
  }
  if ((error==0)&&(ndbleven>0))
  { for (j=0;j<ndbleven*7;++j) devenpar[j]=tempar[j];
  }
  if (tempar) delete [] tempar;

  if ((msorder>1)&&(natoms>1))
    xa=(float)(ndbleven*100.0/natoms/(natoms-1.0));
  else if ((msorder==1)&&(eatoms>0)&&(natoms>1))
    xa=(float)(ndbleven*100.0/eatoms/(natoms-1.0));
  else xa=0.0f;

  if ((error==0)&&(dispmode>4))
  { Textout conout;
    conout.integer(ndbleven,8); conout.floating(xa,8.2f,256);
    conout.space(19); conout.string("ndbleven percentage",0,1);
  }

  if ((error==0)&&(displog>0)&&(flogout))
  { flogout->integer(ndbleven,8); flogout->floating(xa,8.2f,256);
    flogout->space(19); flogout->string("ndbleven percentage",0,1);
    if (error==0) error=flogout->geterror();
  }

  return(error);
} //end of Mscdrun::symdblvert

