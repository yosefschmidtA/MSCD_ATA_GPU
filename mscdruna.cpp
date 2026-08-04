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

int Mscdrun::loadvars(char *ifinname,char *ifoutname,int ijobnum,
  int ijobtotal,int idispmode,int idisplog,Textout *iflogout)
{ if ((error==101)||(error==0))
  { datatype=741; jobnum=ijobnum; jobtotal=ijobtotal;
    dispmode=idispmode; displog=idisplog;
    flogout=iflogout;
    stringcopy(finname,ifinname,100);
    stringcopy(foutname,ifoutname,100);
    error=0;
  }
  return(error);
} //end of Mscdrun::loadvars

float Mscdrun::getmemory()
{ int i;
  float xa,xb;
  if (error==0)
  { xa=(float)((2*80+7*100)*sizeof(char)+101*8*sizeof(int)+
      (14*21+101*20+300*12+npoint*3+nfit*3)*sizeof(float)+
      (natoms*natoms*radim*2+radim)*sizeof(Fcomplex));
    xa+=natoms*natoms*sizeof(int);
    if (msorder>1) xa+=natoms*natoms*(natoms+msorder)*sizeof(int);
    xa+=(numpe+1)*10*sizeof(float);
    xa+=(ntrieven*10+ndbleven*7)*sizeof(float);
    if ((msorder>1)&&(raorder>0))
      xa+=natoms*natoms*natoms*sizeof(float);
    xa+=(ntrielem+ndbleven*radim+natoms*natoms*radim+
      (trymax+50)*(nfit+10))*sizeof(Fcomplex);
    for (i=0;i<katoms;++i) xa+=phaseshift[i].getmemory();
    if (radmatrix) xa+=radmatrix->getmemory();
    if (vibrate) xa+=vibrate->getmemory();
    if (evenmat) xa+=evenmat->getmemory();
    if (termmat) xa+=termmat->getmemory();
    if (expix) xa+=expix->getmemory();
    if (hanka) xa+=hanka->getmemory();
    if (hankb) xa+=hankb->getmemory();
    if (pdintensity) xa+=pdintensity->getmemory();
    xb=(float)(nsymm*sizeof(int)+mscatter*sizeof(float));
    if (xa<xb) xa=xb;
    xa+=basemem+sizeof(Mscdrun);
  }
  else xa=-1.0f;
  return(xa);
} //end of Mscdrun::getmemory

/*
----------------------------------------------------------------------
lakatom   0          1       2          3
          atom-kind  emiter  atomunits  layer-number

layfit    0             1             2              3
          lattice       ---           vinner         tdebye
          origin-vector ---           layer-spacing  scaling-factor

patom     0  1  2  3  4  5  6        7        8  9  10  11  12
          x  y  z  i  j  k  kindatom emitter  ----------------

finals
    0     real calculation for both channels
    1     high channel (li+1) only
    2     low channel (li-1) only
    3     reference only
    4     scattering term only
    5     same as 0
----------------------------------------------------------------------
*/
int Mscdrun::readparameter()
{ int i,j,k,m,n,p,layem,begrow,memsize,scanuni;
  float xa,xb,xc,xd,ya,za,ra,rb,radauto,depauto;
  char *membuf,*namebuf,*fnamebuf,*expfile;
  char *unknown="unknown";
  const float radian=(float)(3.14159265/180.0);

  membuf=namebuf=fnamebuf=expfile=NULL; memsize=256;
  if (error==0) sizeint=sizeof(int); else sizeint=sizeof(int);
  if (error==0)
  { membuf=new char [memsize]; namebuf=new char [memsize];
    fnamebuf=new char [memsize]; expfile=new char [memsize];
    if ((!membuf)||(!namebuf)||(!fnamebuf)||(!expfile)) error=102;
  }
  if (error==0)
  { std::ifstream fin(finname,std::ios::in);
    k=0;
    if (!fin) error=201+(datatype<<16);
    else
    { fin >> xa >> xb >> xc; k=(int)xa; begrow=(int)xb;
      for (i=1;i<begrow;++i) skipendline(fin);
    }
    if ((!fin)||(k!=datatype)) error=207+(datatype<<16);
    for (k=0;k<4;++k) lakatom[k]=-1;
    katoms=0; username[0]=sysname[0]=rmfile[0]=expfile[0]='\0';
    for (i=0;(error==0)&&(i<12);++i)
    { getendline(fin,membuf,memsize);
      spacehide(membuf);
      j=nextnonwhite(membuf);
      k=nextwhite(membuf,j);
      stringcopy(namebuf,membuf+j,k-j+1);
      j=nextnonwhite(membuf,k);
      k=nextwhite(membuf,j);
      stringcopy(fnamebuf,membuf+j,k-j+1);
      spaceback(fnamebuf);
      if (((namebuf[0]=='p')||(namebuf[0]=='P'))&&
        ((namebuf[1]=='s')||(namebuf[1]=='S')))
      { j=0;
        while ((j==lakatom[0])||(j==lakatom[1])||(j==lakatom[2])||
          (j==lakatom[3])) ++j;
        if ((namebuf[2]>='1')&&(namebuf[2]<='4'))
          k=namebuf[2]-'1';
        else k=j;
        if ((k==lakatom[0])||(k==lakatom[1])||(k==lakatom[2])||
          (k==lakatom[3])) k=j;
        if ((k>=0)&&(k<4))
        { lakatom[k]=k;
          if (katoms<k+1) katoms=k+1;
          stringcopy((psfile+k*100),fnamebuf,100);
        }
      }
      else if (((namebuf[0]=='u')||(namebuf[0]=='U'))&&
        ((namebuf[1]=='n')||(namebuf[1]=='N')))
        stringcopy(username,fnamebuf,40);
      else if (((namebuf[0]=='s')||(namebuf[0]=='S'))&&
        ((namebuf[1]=='n')||(namebuf[1]=='N')))
        stringcopy(sysname,fnamebuf,40);
      else if (((namebuf[0]=='r')||(namebuf[0]=='R'))&&
        ((namebuf[1]=='m')||(namebuf[1]=='M'))&&(!radmatrix))
        stringcopy(rmfile,fnamebuf,100);
      else if (((namebuf[0]=='e')||(namebuf[0]=='E'))&&
        ((namebuf[1]=='x')||(namebuf[1]=='X')))
        stringcopy(expfile,fnamebuf,40);
      else if (((namebuf[0]=='p')||(namebuf[0]=='P'))&&
        ((namebuf[1]=='d')||(namebuf[1]=='D')||(namebuf[1]=='e')||
        (namebuf[1]=='E')))
        i=12;
    }
    if ((error==0)&&(katoms<1)) error=201+(711<<16);
    else if (error==0)
    { for (i=0;i<katoms;++i)
      { if ((lakatom[i]<0)||(lakatom[i]>=katoms))
          error=201+(711<<16);
      }
    }
    if (error==0)
    { fin >> xa >> xb >> ftolerance; skipendline(fin);
      scanmode=(int)xa; k=(int)xb;
      fin >> xa >> xb >> xc >> xd; skipendline(fin);
      linitial=(int)xa; lnum=(int)xb; msorder=(int)xc;
      raorder=(int)xd;
      fin >> xa >> xb >> xc >> xd; skipendline(fin);
      totlayer=(int)xa; finals=(int)xb; fitmath=(int)xc;
      trymax=(int)xd;
      fin >> kmin >> kmax >> kstep; skipendline(fin);
      fin >> dtmin >> dtmax >> dtstep; skipendline(fin);
      fin >> dpmin >> dpmax >> dpstep; skipendline(fin);
      fin >> ltheta >> lphi >> xc; skipendline(fin);
      beampol=(int)xc;
      fin >> mtheta >> mphi >> accepang; skipendline(fin);
      fin >> radius >> depth >> lattice; skipendline(fin);
      fin >> valence >> bandgap >> density >> mweight;
      skipendline(fin);
      for (i=0;i<katoms;++i) fin >> aweight[i];
      skipendline(fin);
      fin >> ATA >> ATAweight1 >> ATAweight2;
      skipendline(fin);
      for (i=0;i<katoms;++i) fin >> magamp[i];
      skipendline(fin);
      fin >> vinner >> tdebye >> tsample >> pathcut;
      skipendline(fin);
      fin >> layfit[100*4+2] >> layfit[100*4+3] >> layfit[100*4];
      skipendline(fin); layfit[100*4+1]=0.0f;
      if ((error==0)&&(!fin)) error=207+(datatype<<16);
      else if ((error==0)&&(totlayer<1)) error=601;
    }
    for (i=0;(error==0)&&(i<totlayer)&&(i<100);++i)
    { fin >> xa >> xb >> xc >> xd; skipendline(fin);
      lakatom[i*4+3]=(int)xa; lakatom[i*4]=(int)xb;
      lakatom[i*4+1]=(int)xc; lakatom[i*4+2]=(int)xd;
      fin >> xa >> xb >> xc >> xd; skipendline(fin);
      layatom[i*4]=(int)xa; layatom[i*4+1]=(int)xb;
      layatom[i*4+2]=(int)xc; layatom[i*4+3]=(int)xd;
      fin >> laycell[i*4] >> laycell[i*4+1];
      skipendline(fin);
      fin >> laycell[i*4+2] >> laycell[i*4+3];
      skipendline(fin);
      fin >> layorig[i*4] >> layorig[i*4+1];
      skipendline(fin);
      fin >> layorig[i*4+2]; skipendline(fin);
      fin >> layfit[i*4+2] >> layfit[i*4] >> layfit[i*4+3];
      skipendline(fin);
      layorig[i*4+3]=1.0f; layfit[i*4+1]=0.0f;
      if (!fin) error=207+(datatype<<16);
      else if ((lakatom[i*4]<1)||(lakatom[i*4]>katoms))
        error=201+(711<<16);
    }
  }
  if (error==0)
  { scanmode=confine(scanmode,0,299); totlayer=confine(totlayer,1,100);
    msorder=confine(msorder,0,10); raorder=confine(raorder,-1,4);
    linitial=confine(linitial,0,3); finals=confine(finals,0,5);
    lnum=confine(lnum,0,60); fitmath=confine(fitmath,0,6);
    trymax=confine(trymax,0,10000);
    ltheta=confine(ltheta,0.0f,180.0f);
    lphi=confine(lphi,0.0f,360.0f);
    mtheta=confine(ltheta,0.0f,180.0f);
    mphi=confine(lphi,0.0f,360.0f);
    lattice=confine(lattice,0.1f,50.0f);
    valence=confine(valence,0.0f,100.0f);
    bandgap=confine(bandgap,0.0f,100.0f);
    density=confine(density,0.0f,100.0f);
    mweight=confine(mweight,0.0f,1000.0f);
    radius=confine(radius,0.0f,50.0f);
    depth=confine(depth,0.0f,50.0f);
    tsample=confine(tsample,0.0f,1000.0f);
    tdebye=confine(tdebye,1.0f,1000.0f);
    pathcut=confine(pathcut,0.0f,10.0f);
    vinner=confine(vinner,0.0f,20.0f);
    ftolerance=confine(ftolerance,0.0f,10.0f);
    for (i=0;i<katoms;++i)
      aweight[i]=confine(aweight[i],0.0f,1000.0f);
    for (i=0;i<katoms;++i) magamp[i]=confine(magamp[i],0.0f,1000.0f);
    if ((beampol<1)||(beampol>5)) beampol=1;
    if ((pathcut>1.0)&&(msorder>1)) msorder=1;
    if ((tsample<0.1)||(mweight<0.1)) tsample=0.0f;
    if (expfile[0]=='\0') ftolerance=0.0f;
    if ((sizeint<4)&&(lnum>20)) lnum=20;
    devstep=trymax/1000+1; trymax=trymax%1000;
    trymax=confine(trymax,0,200); devstep=confine(devstep,1,10);
    if (username[0]=='\0')
      stringcopy(username,"Yufeng Chen (LBNL)",40);
  }

  if (error==0)
  { fitmode=0; error=fitcheck();
  }

  if (error==0)
  { radauto=depauto=0.25f;
    if (radius<=radauto) radius=0.0f;
    if (depth<depauto) depth=(float)(radius*(1.5-dtmin/90.0));

    za=0.0f;
    for (i=0;i<totlayer;++i)
    { laycell[i*4]=(float)fabs(laycell[i*4]);
      laycell[i*4+2]=(float)fabs(laycell[i*4+2]);
      layorig[i*4]=(float)fabs(layorig[i*4]);
      layorig[i*4+2]=(float)fabs(layorig[i*4+2]);
      layorig[i*4+3]=(float)fabs(layorig[i*4+3]);
      xa=laycell[i*4]*lattice*(float)cos(laycell[i*4+1]*radian);
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
      if ((radius>radauto)&&(-laxorig[i*4+2]>depth)) totlayer=i;
    }

    layem=-1; m=0;
    for (i=0;i<totlayer;++i)
    { if ((radius<radauto)&&(lakatom[i*4+2]==0))
      { layatom[i*4]=-iabs(layatom[i*4]);
        layatom[i*4+1]=iabs(layatom[i*4+1]);
        layatom[i*4+2]=-iabs(layatom[i*4+2]);
        layatom[i*4+3]=iabs(layatom[i*4+3]);
      }
      else
      { layatom[i*4]=-10;
        layatom[i*4+1]=10;
        layatom[i*4+2]=-10;
        layatom[i*4+3]=10;
      }
      if ((lakatom[i*4+1]!=0)&&(layem<0)) layem=i;
      if (m<lakatom[i*4]) m=lakatom[i*4];
    }
    if (layem<0) layem=0;
    if (m<katoms) katoms=m;
    n=eatoms=0;
    for (k=0;k<totlayer;++k)
    { m=0;
      for (j=layatom[k*4+2];j<=layatom[k*4+3];++j)
      { for (i=layatom[k*4];i<=layatom[k*4+1];++i)
        { xa=laycell[k*4]*lattice; ya=laycell[k*4+2]*lattice;
          if (xa>ya) ra=xa*lakatom[k*4+2];
          else ra=ya*lakatom[k*4+2];
          xa=i*laxcell[k*4]+j*laxcell[k*4+2]+laxorig[k*4];
          ya=i*laxcell[k*4+1]+j*laxcell[k*4+3]+laxorig[k*4+1];
          za=laxorig[k*4+2];
          if (lattice>1.0)
          { xa-=laxorig[layem*4]; ya-=laxorig[layem*4+1];
          }
          rb=0.0f;
          if (n>=300) error=602;
          else if (radius>radauto)
            rb=(xa*xa+ya*ya)/(radius*radius)+za*za/(depth*depth);
          else if ((lakatom[k*4+2]>0)&&(ra>0.1))
            rb=(float)sqrt(xa*xa+ya*ya)/ra;
          else rb=0.0f;
          if (lattice>1.0)
          { xa+=laxorig[layem*4]; ya+=laxorig[layem*4+1];
          }
          if ((error==0)&&(rb<1.001)&&(k>0)&&(lakatom[k*4+3]<0))
          { for (p=0;p<n;++p)
            { if ((fabs(patom[p*12]-xa)<0.1)&&
                (fabs(patom[p*12+1]-ya)<0.1)&&
                (fabs(patom[p*12+2]-laxorig[(k-1)*4+2])<0.001))
              { patom[p*12]=xa; patom[p*12+1]=ya;
                patom[p*12+2]=za;
                patom[p*12+3]=(float)i; patom[p*12+4]=(float)j;
                patom[p*12+5]=(float)k;
                patom[p*12+6]=(float)lakatom[k*4];
                if ((i==0)&&(j==0)&&(patom[p*12+7]!=0)&&
                  (lakatom[k*4+1]==0))
                { patom[p*12+7]=0.0f; --eatoms;
                }
                else if ((i==0)&&(j==0)&&(patom[p*12+7]==0)&&
                  (lakatom[k*4+1]!=0))
                { patom[p*12+7]=(float)lakatom[k*4+1]; ++eatoms;
                }
                else patom[p*12+7]=(float)lakatom[k*4+1];
                p=n+10; --lakatom[(k-1)*4+2];
              }
            }
          }
          else if ((error==0)&&(rb<1.001)) p=0;
          else p=0;

          if ((error==0)&&(rb<1.001)&&(p<=n))
          { patom[n*12]=xa; patom[n*12+1]=ya; patom[n*12+2]=za;
            patom[n*12+3]=(float)i; patom[n*12+4]=(float)j;
            patom[n*12+5]=(float)k;
            patom[n*12+6]=(float)lakatom[k*4];
            if ((i==0)&&(j==0)) patom[n*12+7]=(float)lakatom[k*4+1];
            else patom[n*12+7]=0.0f;
            if (patom[n*12+7]!=0) ++eatoms;
            ++m; ++n;
          }
        }
      }
      lakatom[k*4+2]=m;
    }
    natoms=n;
    if ((error==0)&&(natoms<1)) error=601;
    else if ((error==0)&&(natoms>300)) error=602;
    else if ((error==0)&&(sizeint<4)&&(natoms>12)) error=602;
  }

  if (error==0)
  { eatoms=0;
    for (i=0;i<natoms;++i)
    { if (patom[i*12+7]!=0)
      { k=(int)patom[i*12+5];
        xa=(float)sin((laycell[k*4+3]-laycell[k*4+1])*radian);
        xa=(float)fabs(xa*laycell[k*4]*laycell[k*4+2]);
        if (xa>1.0e-5) xa=1.0f/xa;
        if (xa<1.0e-5) xa=0.0f;
        else ++eatoms;
        patom[i*12+7]=xa;
      }
    }
  }
  if (error==0)
  { nearest=100.0f; biggest=0.0f;
    for (i=0;i<natoms;++i)
    { for (j=i+1;j<natoms;++j)
      { xa=patom[i*12]-patom[j*12]; ya=patom[i*12+1]-patom[j*12+1];
        za=patom[i*12+2]-patom[j*12+2];
        ra=(float)sqrt(xa*xa+ya*ya+za*za);
        if (nearest>ra) nearest=ra;
        if (biggest<ra) biggest=ra;
      }
    }
  }
  if ((error==0)&&(nearest<1.0)) error=605;
  if ((error==0)&&(biggest-nearest<0.001)) biggest=nearest+nearest;

  if ((error==0)&&(natoms<2)&&(msorder>1)) msorder=1;
  if ((error==0)&&(linitial>0))
  { if (radmatrix) delete radmatrix; radmatrix=new Radialmatrix;
    if (linitial==0) stringcopy(fnamebuf,"1s",3);
    else if (linitial==1) stringcopy(fnamebuf,"2p",3);
    else if (linitial==2) stringcopy(fnamebuf,"3d",3);
    else if (linitial==3) stringcopy(fnamebuf,"4f",3);
    else fnamebuf[0]='\0';
    if (!radmatrix) error=102;
    else
    { radmatrix->init();
      error=radmatrix->loadcurve(rmfile,unknown,unknown,fnamebuf);
    }
    if (error==0) error=radmatrix->kconfine(&kmin,&kmax);
  }
  if (error==0)
  { if (phaseshift) delete phaseshift;
    phaseshift=new Phaseshift [katoms];
  }
  if ((error==0)&&(!phaseshift)) error=102;
  for (i=0;i<katoms;++i)
  { if (error==0) phaseshift[i].init();
    if (error==0) error=phaseshift[i].loadcurve((psfile+i*100),
      unknown,unknown,lnum);
  }
  for (i=0;i<katoms;++i)
  { if (error==0) error=phaseshift[i].kconfine(&kmin,&kmax,&lnum);
  }
  if (error==0)
  { if (meanpath) delete meanpath;
    meanpath=new Meanpath; meanpath->init();
  }
  if (error==0)
    error=meanpath->loadparameter(valence,bandgap,density,mweight);

  if (error==0)
  { fitmode=0; error=fitcheck();
  }
  if (error==0)
  { scanuni=scanmode%10;

    if ((error==0)&&(eatoms<1)) error=604;
    else if ((error==0)&&(npoint>30000)) error=617;
    else if ((error==0)&&(sizeint<4)&&(npoint>1000)) error=617;
    if ((error==0)&&(fitmode>0)&&((scanuni==1)||(scanuni==5)))
      error=pdintensity->loadintensity(expfile,expfile,scanmode,
        &kmin,&kmax,&kstep,&dtmin,&dpmin);
    else if ((error==0)&&(fitmode>0)&&((scanuni==2)||(scanuni==6)))
      error=pdintensity->loadintensity(expfile,expfile,scanmode,
        &dtmin,&dtmax,&dtstep,&kmin,&dpmin);
    else if ((error==0)&&(fitmode>0)&&((scanuni==3)||(scanuni==7)))
      error=pdintensity->loadintensity(expfile,expfile,scanmode,
        &dpmin,&dpmax,&dpstep,&kmin,&dtmin);
    else if (error==0)
      error=pdintensity->loadparameter(scanmode,kmin,kmax,kstep,
        dtmin,dtmax,dtstep,dpmin,dpmax,dpstep);
    if (error==0) error=pdintensity->setbeamangle(ltheta,lphi);
    if (error==0)
    { npoint=pdintensity->getnumpoint();
      ncurve=pdintensity->getnumcurve();
    }
  }

  if (error==0) error=fitcheck();

  if (error==0)
  { if ((msorder<1)||(lnum==0)||(finals==3)||
      ((linitial==0)&&(finals==2))) msorder=raorder=lnum=0;
    if (finals==0) finals=5;
    if (raorder>lnum-1) raorder=lnum-1;

    if (raorder==4) radim=15;
    else if (raorder==3) radim=10;
    else if (raorder==2) radim=6;
    else if (raorder==1) radim=3;
    else radim=1;

    if ((fitmath==0)&&(fitmode>=4)) fitmath=2;
    if ((trymax==0)&&(fitmath>4)) trymax=100;
    else if ((trymax==0)&&(fitmath>3)&&(nfit-nfollow<2)) trymax=10;
    else if ((trymax==0)&&(fitmath>3)&&(nfit-nfollow<3)) trymax=26;
    else if ((trymax==0)&&(fitmath>3)) trymax=126;
    else if ((trymax==0)&&(fitmath>2)&&(nfit-nfollow<2)) trymax=50;
    else if ((trymax==0)&&(fitmath>2)&&(nfit-nfollow<3)) trymax=75;
    else if ((trymax==0)&&(fitmath>2))
      trymax=125+25*(nfit-nfollow)*(nfit-nfollow);
    else if ((trymax==0)&&(fitmath>1)) trymax=(nfit-nfollow)*40+1;
    else if (trymax==0) trymax=(nfit-nfollow)*25;
    if (trymax>200) trymax=200;
  }

  if ((error==0)&&((numpe>npoint)||(numpe>128))) error=761;
  else if (error==0)
  { if (tusage) delete [] tusage;
    tusage=new float [(numpe+1)*10];
    if (tusage)
    { for (i=0;i<(numpe+1)*10;++i) tusage[i]=0.0f;
      tusage[mype*10]=(float)mype;
    }
    else error=102;
  }

  if (membuf) delete [] membuf; if (namebuf) delete [] namebuf;
  if (fnamebuf) delete [] fnamebuf; if (expfile) delete [] expfile;

  return(error);
} //end of Mscdrun::readparameter

/*
----------------------------------------------------------------------
layfit    0             1             2              3
          lattice       ---           vinner         tdebye
          origin-vector ---           layer-spacing  scaling-factor

layorig   0             1             2              3
          lattice       ---           vinner         tdebye
          origin-vector origin-angle  layer-spacing  scaling-factor

fitmath
    0     automatic (currently = 2)
    1     non-linear marquadt fitting only
    2     simplex downhill fitting, then non-linear marquadt
    3     grid search, then simplex downhill, then non-linear marquadt
    4     grid search only
    5     search emission angle deviation only

fitmode
    0     no fitting at all
    1     no fitting, but do calculate reliability factor
    2     fitting emission angle deviation
    3     fitting beam angle deviation (not implemented here)
    4     non-structural fitting
    5     structural fitting
    6     both non-structural and structural fitting
----------------------------------------------------------------------
*/

int Mscdrun::fitcheck()
{ int i,j,k,m,scanhun,scanten,scanuni;
  float atheta,ptstep;
  const float radian=(float)(3.14159265/180.0), localeps=1.0e-3f;

  scanhun=scanten=scanuni=0;
  if (error==0)
  { scanhun=(scanmode%1000)/100; scanten=(scanmode%100)/10;
    scanuni=scanmode%10;

    kmin=confine(kmin,3.0f,25.0f); kmax=confine(kmax,kmin,25.0f);
    if ((kstep<0.001)||(kstep>kmax-kmin))
    { kmax=kmin; kstep=0.0f;
    }
    accepang=confine(accepang,0.0f,5.0f);
    dtmin=confine(dtmin,0.0f,89.0f-accepang);
    dtmax=confine(dtmax,dtmin,89.0f-accepang);
    if ((dtstep<0.1)||(dtstep>dtmax-dtmin))
    { dtmax=dtmin; dtstep=0.0f;
    }
    dpmin=confine(dpmin,0.0f,360.0f);
    dpmax=confine(dpmax,dpmin,360.0f);
    if (dtmax<0.1) dpmax=dpmin=dpstep=0.0f;
    else if ((dpstep<0.1)||(dpstep>dpmax-dpmin))
    { dpmax=dpmin; dpstep=0.0f;
    }

    if (kstep>0.0) nkout=(int)((kmax-kmin)/kstep+1.001); else nkout=1;
    nkout=confine(nkout,1,256); kmax=kmin+kstep*(nkout-1);
    if (dtstep>0.0) ndtheta=(int)((dtmax-dtmin)/dtstep+1.001);
    else ndtheta=1;
    ndtheta=confine(ndtheta,1,256); dtmax=dtmin+dtstep*(ndtheta-1);
    if (dpstep>0.0) ndphi=(int)((dpmax-dpmin)/dpstep+1.001);
    else ndphi=1;
    ndphi=confine(ndphi,1,256); dpmax=dpmin+dpstep*(ndphi-1);

    if ((scanuni==2)||(scanuni==3)||(scanuni==6)||(scanuni==7)||
      (ndtheta==1)||(ndphi==1)) ndangle=ndtheta*ndphi;
    else
    { ndangle=0;
      for (i=0;i<ndtheta;++i)
      { atheta=(dtmin+i*dtstep);
        if (atheta<1.0e-3) ++ndangle;
        else
        { ptstep=dpstep/(float)sin(atheta*radian);
          k=(int)((dpmax-dpmin)/ptstep+1.5);
          ndangle+=k;
        }
      }
    }
    if (fitmode==0) npoint=nkout*ndangle;
    if ((fitmode>0)&&(scanuni!=1)&&(scanuni!=5)&&(nkout>ncurve))
      nkout=ncurve;
    if ((fitmode>0)&&(scanuni!=2)&&(scanuni!=6)&&(scanuni!=4)&&
      (scanuni!=8)&&(ndtheta>ncurve))
      ndtheta=ncurve;
    if ((fitmode>0)&&(scanuni!=3)&&(scanuni!=7)&&(scanuni!=4)&&
      (scanuni!=8)&&(ndphi>ncurve))
      ndphi=ncurve;
    if ((fitmode>0)&&((scanuni==1)||(scanuni==5))) ndangle=ncurve;

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
  }

  if (error==0)
  { layorig[totlayer*4]=lattice; layorig[totlayer*4+1]=0.0f;
    layorig[totlayer*4+2]=vinner; layorig[totlayer*4+3]=tdebye;

    for (i=0;i<totlayer+1;++i)
    { for (j=0;j<4;++j)
      { if (fitmath>4) layfit[i*4+j]=0.0f;
        else if (j==1) layfit[i*4+j]=0.0f;
        else if (layorig[i*4+j]<localeps) layfit[i*4+j]=0.0f;
        else if ((i==0)&&(layfit[i*4+j]<localeps)) layfit[i*4+j]=0.0f;
        else if ((i<totlayer)&&(fabs(layfit[i*4+j])<localeps))
          layfit[i*4+j]=0.0f;
        else if ((i==totlayer)&&(layfit[100*4+j]<localeps))
          layfit[i*4+j]=layfit[100*4+j]=0.0f;
        else if ((i==totlayer)&&(j==3)&&(tsample<1.0))
          layfit[i*4+j]=layfit[100*4+j]=0.0f;
        else if ((i==totlayer)&&(layfit[100*4+j]>0.5))
          layfit[i*4+j]=layfit[100*4+j]=0.5f;
        else if (i==totlayer) layfit[i*4+j]=layfit[100*4+j];
        else if (layfit[i*4+j]>0.5) layfit[i*4+j]=0.5f;
        else if (layfit[i*4+j]<-0.5) layfit[i*4+j]=-0.5f;
      }
    }

    nfit=nfollow=0;
    if ((ftolerance<1.0e-5)||(scanuni==4)||(scanuni==8)||
      (trymax==1))
    { ftolerance=0.0f; fitmath=fitmode=0; nfit=trymax=0;
      for (i=0;i<totlayer+1;++i)
      { for (j=0;j<4;++j) layfit[i*4+j]=0.0f;
      }
      for (j=0;j<4;++j) layfit[100*4+j]=0.0f;
    }
    else if ((ftolerance>2.5)||((trymax>1)&&(trymax<4)))
    { fitmath=0; fitmode=1; nfit=trymax=0;
      for (i=1;i<totlayer+1;++i)
      { for (j=0;j<4;++j) layfit[i*4+j]=0.0f;
      }
      for (j=0;j<4;++j) layfit[100*4+j]=0.0f;
    }
    else if (fitmath<5)
    { fitmode=1;
      for (j=0;j<4;++j)
      { i=totlayer;
        if (layfit[i*4+j]>0.0)
        { fitvars[nfit]=layorig[i*4+j];
          fitvars[20+nfit]=layfit[i*4+j];
          fitvars[40+nfit]=(float)j; fitvars[60+nfit]=(float)i;
          ++nfit;
          if (j<2) fitmode=5;
          else if (fitmode<5) fitmode=4;
          else if (fitmode==5) fitmode=6;
        }
      }
      for (i=0;i<totlayer;++i)
      { for (j=0;j<4;++j)
        { if ((nfit<10)&&(layfit[i*4+j]>0.0))
          { fitvars[nfit]=layorig[i*4+j];
            fitvars[20+nfit]=layfit[i*4+j];
            fitvars[40+nfit]=(float)j; fitvars[60+nfit]=(float)i;
            ++nfit;
            if (fitmode<4) fitmode=5;
            else if (fitmode==4) fitmode=6;
          }
          else if (layfit[i*4+j]<0.0)
          { for (k=i-1;k>=0;--k) if (layfit[k*4+j]>0.0) break;
            if (k>=0)
            { fitvars[nfit]=layorig[i*4+j];
              fitvars[40+nfit]=(float)j; fitvars[60+nfit]=(float)i;
              for (m=nfit;m>=0;--m)
                if ((fitvars[40+m]==j)&&(fitvars[60+m]==k)) break;
              if (m>=0)
              { layfit[i*4+j]=-(float)(m+1);
                fitvars[20+nfit]=layfit[i*4+j];
                ++nfit; ++nfollow;
              }
            }
          }
          else layfit[i*4+j]=0.0f;
        }
      }
    }
    else if ((fitmath>4)&&(scanuni!=1)&&(scanuni!=5))
    { nfit=0; fitmath=0; fitmode=1;
    }
    else if ((fitmath==5)&&(dtmax+accepang<85.0))
    { nfit=2; fitmode=2;
      fitvars[0]=fitvars[1]=0.0f;
      fitvars[20]=(float)devstep; fitvars[21]=89.0f-dtmax-accepang;
      for (i=2;i<4;++i) for (j=0;j<nfit;++j)
        fitvars[i*20+j]=0.0f;
    }
    else if (fitmath>4) fitmath=0;
    if (nfit==0) fitmath=0;
    if ((dispmode==0)&&(fitmode<2)) dispmode=4;
    else if (dispmode==0) dispmode=3;
  }
  if (error==0)
  { for (j=0;j<nfit;++j) fitvars[nfit+j]=fitvars[20+j];
    for (j=0;j<nfit;++j) fitvars[nfit*2+j]=fitvars[40+j];
    for (j=0;j<nfit;++j) fitvars[nfit*3+j]=fitvars[60+j];
  }

  return(error);
} //end of Mscdrun::fitcheck

