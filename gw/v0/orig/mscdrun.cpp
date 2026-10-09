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

Mscdrun::Mscdrun(int imype,int inumpe)
{ mype=imype; numpe=inumpe;
  error=107; init();
} //end of Mscdrun::Mscdrun

Mscdrun::~Mscdrun()
{ if (username) delete [] username; if (sysname) delete [] sysname;
  if (finname) delete [] finname; if (foutname) delete [] foutname;
  if (rmfile) delete [] rmfile; if (psfile) delete [] psfile;
  if (lakatom) delete [] lakatom; if (layatom) delete [] layatom;
  if (pdnum) delete [] pdnum; if (tevenadd) delete [] tevenadd;
  if (tevendim) delete [] tevendim; if (tevencut) delete [] tevencut;
  if (devenadd) delete [] devenadd;
  if (laycell) delete [] laycell; if (layorig) delete [] layorig;
  if (layfit)  delete [] layfit; if (laxcell) delete [] laxcell;
  if (laxorig) delete [] laxorig; if (fitvars) delete [] fitvars;
  if (aweight) delete [] aweight; if (magamp) delete [] magamp;
  if (patom) delete [] patom; if (tevenpar) delete [] tevenpar;
  if (devenpar) delete [] devenpar; if (talpha) delete [] talpha;
  if (tgamma) delete [] tgamma; if (fithist) delete [] fithist;
  if (tusage) delete [] tusage;
  if (tevenelem) delete [] tevenelem;
  if (devenelem) delete [] devenelem;
  if (devendetec) delete [] devendetec;
  if (phaseshift) delete [] phaseshift;
  if (radmatrix) delete radmatrix; if (meanpath) delete meanpath;
  if (vibrate) delete vibrate;
  if (evenmat) delete evenmat; if (termmat) delete termmat;
  if (expix) delete expix;
  if (hanka) delete hanka; if (hankb) delete hankb;
  if (pdintensity) delete pdintensity;
  if (jobtime) delete jobtime;
} //end of Mscdrun::~Mscdrun

void Mscdrun::init()
{ if ((error==107)||(error==103))
  { datatype=741; jobnum=jobtotal=dispmode=displog=0;
    flogout=NULL;
    username=sysname=finname=foutname=rmfile=psfile=NULL;
    lakatom=layatom=pdnum=tevenadd=tevendim=tevencut=devenadd=NULL;
    fitvars=laycell=layorig=layfit=laxcell=laxorig=aweight=magamp=
      patom=tevenpar=devenpar=talpha=tgamma=fithist=tusage=NULL;
    tevenelem=devenelem=devendetec=NULL;
    phaseshift=NULL; radmatrix=NULL; meanpath=NULL; vibrate=NULL;
    evenmat=termmat=NULL; expix=NULL; hanka=hankb=NULL;
    pdintensity=NULL;
    nsymm=mscatter=0; pdbeg=pdend=0; basemem=0.0f;
  }
  if (error==103)
  { username=new char [80]; sysname=new char [80];
    finname=new char [100]; foutname=new char [100];
    rmfile=new char [100]; psfile=new char [100*4];
    lakatom=new int [101*4]; layatom=new int [101*4];
    laycell=new float [101*4]; layorig=new float [101*4];
    layfit=new float [101*4]; laxcell=new float [101*4];
    laxorig=new float [101*4]; fitvars=new float [4*20];
    aweight=new float [4]; magamp=new float [4];
    patom=new float [300*12];
    pdintensity=new Pdintensity;
    if ((!username)||(!sysname)||(!finname)||(!foutname)||(!rmfile)||
      (!psfile)||(!lakatom)||(!layatom)||(!fitvars)||(!laycell)||
      (!layorig)||(!layfit)||(!laxcell)||(!laxorig)||(!aweight)||
      (!magamp)||(!patom)||(!pdintensity)) error=102;
    else pdintensity->init();
  }
  else if (error==106)
  { username=new char [80]; sysname=new char [80];
    finname=new char [100]; foutname=new char [100];
    rmfile=new char [100]; psfile=new char [100*4];
    lakatom=new int [101*4]; layatom=new int [101*4];
    laycell=new float [101*4]; layorig=new float [101*4];
    layfit=new float [101*4]; laxcell=new float [101*4];
    laxorig=new float [101*4]; fitvars=new float [4*20];
    aweight=new float [4]; magamp=new float [4];
    patom=new float [300*12];
    if ((!username)||(!sysname)||(!finname)||(!foutname)||(!rmfile)||
      (!psfile)||(!lakatom)||(!layatom)||(!fitvars)||(!laycell)||
      (!layorig)||(!layfit)||(!laxcell)||(!laxorig)||(!aweight)||
      (!magamp)||(!patom)||(!pdintensity)) error=102;

    if (pdnum)
    { pdnum=new int [npoint]; if (!pdnum) error=102;
    }
    if (tevenadd)
    { tevenadd=new int [natoms*natoms*natoms];
      if (!tevenadd) error=102;
    }
    if (tevendim)
    { tevendim=new int [natoms*natoms*natoms];
      if (!tevendim) error=102;
    }
    if (tevencut)
    { tevencut=new int [msorder*natoms*natoms];
      if (!tevencut) error=102;
    }
    if (devenadd)
    { devenadd=new int [natoms*natoms];
      if (!devenadd) error=102;
    }
    if (tevenpar)
    { tevenpar=new float [ntrieven*10];
      if (!tevenpar) error=102;
    }
    if (devenpar)
    { devenpar=new float [ndbleven*7];
      if (!devenpar) error=102;
    }
    if (talpha)
    { talpha=new float [natoms*natoms*natoms];
      if (!talpha) error=102;
    }
    if (tgamma)
    { tgamma=new float [natoms*natoms*natoms];
      if (!tgamma) error=102;
    }
    if (fithist)
    { fithist=new float [(trymax+50)*(nfit+10)];
      if (!fithist) error=102;
    }
    if (tusage)
    { tusage=new float [(numpe+1)*10];
      if (!tusage) error=102;
    }
    if (tevenelem)
    { tevenelem=new Fcomplex [ntrielem];
      if (!tevenelem) error=102;
    }
    if (devenelem)
    { devenelem=new Fcomplex [ndbleven*radim];
      if (!devenelem) error=102;
    }
    if (devendetec)
    { devendetec=new Fcomplex [natoms*natoms*radim];
      if (!devendetec) error=102;
    }
  }
  else if (error==107)
  { jobtime=new Jobtime;
    if (!jobtime) error=102;
    else jobtime->reset();
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Mscdrun::init

int Mscdrun::getlength()
{ int i,ka;
  if (error==0)
  { ka=sizeof(int)+sizeof(Mscdrun)+
      (2*80+7*100)*sizeof(char)+2*101*4*sizeof(int)+
      (5*101*4+4*20+2*4+300*12)*sizeof(float);

    if (pdnum) ka+=npoint*sizeof(int);
    if (tevenadd) ka+=natoms*natoms*natoms*sizeof(int);
    if (tevendim) ka+=natoms*natoms*natoms*sizeof(int);
    if (tevencut) ka+=msorder*natoms*natoms*sizeof(int);
    if (devenadd) ka+=natoms*natoms*sizeof(int);
    if (tevenpar) ka+=ntrieven*10*sizeof(float);
    if (devenpar) ka+=ndbleven*7*sizeof(float);
    if (talpha) ka+=natoms*natoms*natoms*sizeof(float);
    if (tgamma) ka+=natoms*natoms*natoms*sizeof(float);
    if (fithist) ka+=(trymax+50)*(nfit+10)*sizeof(float);
    if (tusage) ka+=(numpe+1)*10*sizeof(float);
    if (tevenelem) ka+=ntrielem*sizeof(Fcomplex);
    if (devenelem) ka+=ndbleven*radim*sizeof(Fcomplex);
    if (devendetec) ka+=natoms*natoms*radim*sizeof(Fcomplex);

    for (i=0;i<katoms;++i) ka+=phaseshift[i].getlength();
    if (radmatrix) ka+=radmatrix->getlength();
    if (meanpath) ka+=meanpath->getlength();
    if (vibrate) ka+=vibrate->getlength();
    if (evenmat) ka+=evenmat->getlength();
    if (termmat) ka+=termmat->getlength();
    if (expix) ka+=expix->getlength();
    if (hanka) ka+=hanka->getlength();
    if (hankb) ka+=hankb->getlength();
    if (pdintensity) ka+=pdintensity->getlength();
  }
  else ka=0;
  return(ka);
} //end of Mscdrun::getlength

int Mscdrun::paexport(char *dest,int length)
{ int i,ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Mscdrun);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=80;
    if (ka>=0) ka=memorysend(dest,username,ka,kb,length);
    kb=80;
    if (ka>=0) ka=memorysend(dest,sysname,ka,kb,length);
    kb=100;
    if (ka>=0) ka=memorysend(dest,finname,ka,kb,length);
    kb=100;
    if (ka>=0) ka=memorysend(dest,foutname,ka,kb,length);
    kb=100;
    if (ka>=0) ka=memorysend(dest,rmfile,ka,kb,length);
    kb=100*4;
    if (ka>=0) ka=memorysend(dest,psfile,ka,kb,length);
    kb=101*4*sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)lakatom,ka,kb,length);
    kb=101*4*sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)layatom,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)laycell,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)layorig,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)layfit,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)laxcell,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)laxorig,ka,kb,length);
    kb=4*20*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)fitvars,ka,kb,length);
    kb=4*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)aweight,ka,kb,length);
    kb=4*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)magamp,ka,kb,length);
    kb=300*12*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)patom,ka,kb,length);

    kb=npoint*sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)pdnum,ka,kb,length);
    kb=natoms*natoms*natoms*sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)tevenadd,ka,kb,length);
    kb=natoms*natoms*natoms*sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)tevendim,ka,kb,length);
    kb=msorder*natoms*natoms*sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)tevencut,ka,kb,length);
    kb=natoms*natoms*sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)devenadd,ka,kb,length);
    kb=ntrieven*10*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)tevenpar,ka,kb,length);
    kb=ndbleven*7*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)devenpar,ka,kb,length);
    kb=natoms*natoms*natoms*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)talpha,ka,kb,length);
    kb=natoms*natoms*natoms*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)tgamma,ka,kb,length);
    kb=(trymax+50)*(nfit+10)*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)fithist,ka,kb,length);
    kb=(numpe+1)*10*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)tusage,ka,kb,length);
    kb=ntrielem*sizeof(Fcomplex);
    if (ka>=0) ka=memorysend(dest,(char *)tevenelem,ka,kb,length);
    kb=ndbleven*radim*sizeof(Fcomplex);
    if (ka>=0) ka=memorysend(dest,(char *)devenelem,ka,kb,length);
    kb=natoms*natoms*radim*sizeof(Fcomplex);
    if (ka>=0) ka=memorysend(dest,(char *)devendetec,ka,kb,length);
    if (ka<0) error=901;

    for (i=0;(error==0)&&(phaseshift)&&(i<katoms);++i)
    { kb=phaseshift[i].getlength();
      if (ka+kb<=length) error=phaseshift[i].paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(radmatrix))
    { kb=radmatrix->getlength();
      if (ka+kb<=length) error=radmatrix->paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(meanpath))
    { kb=meanpath->getlength();
      if (ka+kb<=length) error=meanpath->paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(vibrate))
    { kb=vibrate->getlength();
      if (ka+kb<=length) error=vibrate->paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(evenmat))
    { kb=evenmat->getlength();
      if (ka+kb<=length) error=evenmat->paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(termmat))
    { kb=termmat->getlength();
      if (ka+kb<=length) error=termmat->paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(expix))
    { kb=expix->getlength();
      if (ka+kb<=length) error=expix->paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(hanka))
    { kb=hanka->getlength();
      if (ka+kb<=length) error=hanka->paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(hankb))
    { kb=hankb->getlength();
      if (ka+kb<=length) error=hankb->paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(pdintensity))
    { kb=pdintensity->getlength();
      if (ka+kb<=length) error=pdintensity->paexport(dest+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Mscdrun::paexport

int Mscdrun::paimport(char *source,int length)
{ int i,ka,kb,ma,mb;
  void *temp;

  ma=mype; mb=numpe; temp=(void *)jobtime;
  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Mscdrun);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=80;
    if (ka>=0) ka=memoryrec(username,source,ka,kb,length);
    kb=80;
    if (ka>=0) ka=memoryrec(sysname,source,ka,kb,length);
    kb=100;
    if (ka>=0) ka=memoryrec(finname,source,ka,kb,length);
    kb=100;
    if (ka>=0) ka=memoryrec(foutname,source,ka,kb,length);
    kb=100;
    if (ka>=0) ka=memoryrec(rmfile,source,ka,kb,length);
    kb=100*4;
    if (ka>=0) ka=memoryrec(psfile,source,ka,kb,length);
    kb=101*4*sizeof(int);
    if (ka>=0) ka=memoryrec((char *)lakatom,source,ka,kb,length);
    kb=101*4*sizeof(int);
    if (ka>=0) ka=memoryrec((char *)layatom,source,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)laycell,source,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)layorig,source,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)layfit,source,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)laxcell,source,ka,kb,length);
    kb=101*4*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)laxorig,source,ka,kb,length);
    kb=4*20*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)fitvars,source,ka,kb,length);
    kb=4*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)aweight,source,ka,kb,length);
    kb=4*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)magamp,source,ka,kb,length);
    kb=300*12*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)patom,source,ka,kb,length);

    kb=npoint*sizeof(int);
    if (ka>=0) ka=memoryrec((char *)pdnum,source,ka,kb,length);
    kb=natoms*natoms*natoms*sizeof(int);
    if (ka>=0) ka=memoryrec((char *)tevenadd,source,ka,kb,length);
    kb=natoms*natoms*natoms*sizeof(int);
    if (ka>=0) ka=memoryrec((char *)tevendim,source,ka,kb,length);
    kb=msorder*natoms*natoms*sizeof(int);
    if (ka>=0) ka=memoryrec((char *)tevencut,source,ka,kb,length);
    kb=natoms*natoms*sizeof(int);
    if (ka>=0) ka=memoryrec((char *)devenadd,source,ka,kb,length);
    kb=ntrieven*10*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)tevenpar,source,ka,kb,length);
    kb=ndbleven*7*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)devenpar,source,ka,kb,length);
    kb=natoms*natoms*natoms*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)talpha,source,ka,kb,length);
    kb=natoms*natoms*natoms*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)tgamma,source,ka,kb,length);
    kb=(trymax+50)*(nfit+10)*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)fithist,source,ka,kb,length);
    kb=(numpe+1)*10*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)tusage,source,ka,kb,length);
    kb=ntrielem*sizeof(Fcomplex);
    if (ka>=0) ka=memoryrec((char *)tevenelem,source,ka,kb,length);
    kb=ndbleven*radim*sizeof(Fcomplex);
    if (ka>=0) ka=memoryrec((char *)devenelem,source,ka,kb,length);
    kb=natoms*natoms*radim*sizeof(Fcomplex);
    if (ka>=0) ka=memoryrec((char *)devendetec,source,ka,kb,length);
    if (ka<0) error=901;

    if ((error==0)&&(phaseshift))
    { phaseshift=new Phaseshift [katoms];
      if (!phaseshift) error=102;
    }
    for (i=0;(error==0)&&(phaseshift)&&(i<katoms);++i)
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      if (ka+kb<=length) error=phaseshift[i].paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(radmatrix))
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      radmatrix=new Radialmatrix;
      if (!radmatrix) error=102;
      else if (ka+kb<=length) error=radmatrix->paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(meanpath))
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      meanpath=new Meanpath;
      if (!meanpath) error=102;
      else if (ka+kb<=length) error=meanpath->paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(vibrate))
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      vibrate=new Vibration;
      if (!vibrate) error=102;
      else if (ka+kb<=length) error=vibrate->paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(evenmat))
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      evenmat=new Rotamat;
      if (!evenmat) error=102;
      else if (ka+kb<=length) error=evenmat->paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(termmat))
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      termmat=new Rotamat;
      if (!termmat) error=102;
      else if (ka+kb<=length) error=termmat->paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(expix))
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      expix=new Expix;
      if (!expix) error=102;
      else if (ka+kb<=length) error=expix->paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(hanka))
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      hanka=new Hankel;
      if (!hanka) error=102;
      else if (ka+kb<=length) error=hanka->paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(hankb))
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      hankb=new Hankel;
      if (!hankb) error=102;
      else if (ka+kb<=length) error=hankb->paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&(pdintensity))
    { memorycopy((char *)&kb,source+ka,sizeof(int));
      pdintensity=new Pdintensity;
      if (!pdintensity) error=102;
      else if (ka+kb<=length) error=pdintensity->paimport(source+ka,kb);
      ka+=kb;
    }
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  mype=ma; numpe=mb; jobtime=(Jobtime *)temp;
  if (tusage)
  { for (i=0;i<(numpe+1)*10;++i) tusage[i]=0.0f;
    tusage[mype*10]=(float)mype;
  }
  return(error);
} //end of Mscdrun::paimport

int Mscdrun::sendjobs()
{ int i,sendsize;
  char *sendbuf;

  sendbuf=NULL;
  if ((numpe>1)&&(mype==0))
  { if (tusage) tusage[mype*10+1]+=jobtime->lastprocess();
    if (error==0)
    { sendsize=getlength(); sendbuf=new char [sendsize];
      if (!sendbuf) error=102;
    }
    else sendsize=0;
    if (error==0) error=paexport(sendbuf,sendsize);
    for (i=1;i<numpe;++i)
    { if (error==0)
        error=mpisend((char *)&sendsize,sizeof(int),i,3072+i);
      else mpisend((char *)&sendsize,sizeof(int),i,3072+i);
      if ((error==0)&&(sendsize>0))
        error=mpisend(sendbuf,sendsize,i,4096+i);
    }
    if (tusage) tusage[mype*10+2]+=jobtime->lastprocess();
  }
  if (sendbuf) delete [] sendbuf;

  return(error);
} //end of Mscdrun::sendjobs

int Mscdrun::receivejobs()
{ int recsize;
  char *recbuf;

  recbuf=NULL;
  if ((numpe>1)&&(mype>0))
  { if (tusage) tusage[mype*10+1]+=jobtime->lastprocess();
    if (error==103)
      error=mpireceive((char *)&recsize,sizeof(int),0,3072+mype);
    else mpireceive((char *)&recsize,sizeof(int),0,3072+mype);
    if ((error==0)&&(recsize>sizeof(int)))
    { recbuf=new char [recsize];
      if (!recbuf) error=102;
      if (error==0) error=mpireceive(recbuf,recsize,0,4096+mype);
    }
    if ((error==0)&&(recsize>sizeof(int)))
    { error=103; error=paimport(recbuf,recsize);
    }
    else if (error==0) error=753;
    if (error==0)
    { dispmode=1; displog=0; flogout=NULL;
    }
    if (tusage) tusage[mype*10+4]+=jobtime->lastprocess();
  }
  if (recbuf) delete [] recbuf;

  return(error);
} //end of Mscdrun::receivejobs

int Mscdrun::getsuplength(int order)
{ int ka;
  if ((error==0)&&(order==1))
  { ka=sizeof(int)+(3+4*20+300*12)*sizeof(float);
    if (pdnum) ka+=npoint*sizeof(int);
    if (tevenpar) ka+=ntrieven*10*sizeof(float);
    if (devenpar) ka+=ndbleven*7*sizeof(float);
    if (talpha) ka+=natoms*natoms*natoms*sizeof(float);
    if (tgamma) ka+=natoms*natoms*natoms*sizeof(float);
  }
  else if ((error==0)&&(order==2))
    ka=3*sizeof(int)+10*sizeof(float)+npoint*3*sizeof(float);
  else if ((error==0)&&(order==0)) ka=sizeof(int);
  else ka=0;
  return(ka);
} //end of Mscdrun::getsuplength

int Mscdrun::pasupexport(char *dest,int length,int order)
{ int i,j,ka,kb;
  float xa,xb,xc,xd;
  float tproc[10];

  if ((error==0)&&(dest))
  { ka=0; kb=getsuplength(order);
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    if (order==1)
    { kb=sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)&vinner,ka,kb,length);
      kb=sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)&tsample,ka,kb,length);
      kb=sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)&lattice,ka,kb,length);
      kb=4*20*sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)fitvars,ka,kb,length);
      kb=300*12*sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)patom,ka,kb,length);

      kb=npoint*sizeof(int);
      if (ka>=0) ka=memorysend(dest,(char *)pdnum,ka,kb,length);
      kb=ntrieven*10*sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)tevenpar,ka,kb,length);
      kb=ndbleven*7*sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)devenpar,ka,kb,length);
      kb=natoms*natoms*natoms*sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)talpha,ka,kb,length);
      kb=natoms*natoms*natoms*sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)tgamma,ka,kb,length);
    }
    else if (order==2)
    { for (i=0;i<10;++i)
      { if (tusage) tproc[i]=tusage[mype*10+i]; else tproc[i]=0.0f;
      }
      kb=10*sizeof(float);
      if (ka>=0) ka=memorysend(dest,(char *)&tproc,ka,kb,length);
      kb=sizeof(int);
      if (ka>=0) ka=memorysend(dest,(char *)&pdbeg,ka,kb,length);
      kb=sizeof(int);
      if (ka>=0) ka=memorysend(dest,(char *)&pdend,ka,kb,length);
      if (ka<0) error=901;
      for (i=0;i<npoint;++i)
      { j=pdnum[i];
        if (error==0) error=pdintensity->getpoint(j,&xd,&xd,&xd,&xd,
          &xd,&xa,&xb,&xc,&xd);
        if (error==0)
        { kb=sizeof(float);
          if (ka>=0) ka=memorysend(dest,(char *)&xa,ka,kb,length);
          kb=sizeof(float);
          if (ka>=0) ka=memorysend(dest,(char *)&xb,ka,kb,length);
          kb=sizeof(float);
          if (ka>=0) ka=memorysend(dest,(char *)&xc,ka,kb,length);
        }
      }
    }
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Mscdrun::pasupexport

int Mscdrun::pasupimport(char *source,int length,int order)
{ int i,j,ka,kb,ma,mb;
  float xa,xb,xc,xd,xe,xf,xg,xh,xi;
  float tproc[10];

  if ((error==0)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<sizeof(int))||(length!=kb)) ka=-1;
    if ((ka>=0)&&(order==1))
    { kb=sizeof(float);
      if (ka>=0) ka=memoryrec((char *)&vinner,source,ka,kb,length);
      kb=sizeof(float);
      if (ka>=0) ka=memoryrec((char *)&tsample,source,ka,kb,length);
      kb=sizeof(float);
      if (ka>=0) ka=memoryrec((char *)&lattice,source,ka,kb,length);
      kb=4*20*sizeof(float);
      if (ka>=0) ka=memoryrec((char *)fitvars,source,ka,kb,length);
      kb=300*12*sizeof(float);
      if (ka>=0) ka=memoryrec((char *)patom,source,ka,kb,length);

      kb=npoint*sizeof(int);
      if (ka>=0) ka=memoryrec((char *)pdnum,source,ka,kb,length);
      kb=ntrieven*10*sizeof(float);
      if (ka>=0) ka=memoryrec((char *)tevenpar,source,ka,kb,length);
      kb=ndbleven*7*sizeof(float);
      if (ka>=0) ka=memoryrec((char *)devenpar,source,ka,kb,length);
      kb=natoms*natoms*natoms*sizeof(float);
      if (ka>=0) ka=memoryrec((char *)talpha,source,ka,kb,length);
      kb=natoms*natoms*natoms*sizeof(float);
      if (ka>=0) ka=memoryrec((char *)tgamma,source,ka,kb,length);
    }
    else if ((ka>=0)&&(order==2))
    { kb=10*sizeof(float);
      if (ka>=0) ka=memoryrec((char *)tproc,source,ka,kb,length);
      kb=sizeof(int);
      if (ka>=0) ka=memoryrec((char *)&ma,source,ka,kb,length);
      kb=sizeof(int);
      if (ka>=0) ka=memoryrec((char *)&mb,source,ka,kb,length);
      if ((ka>=0)&&(ma>=0)) ka+=ma*3*sizeof(float);
      if (ka<0) error=901;
      for (i=ma;i<mb;++i)
      { j=pdnum[i];
        if (error==0) error=pdintensity->getpoint(j,&xa,&xb,&xc,&xd,
          &xe,&xf,&xg,&xh,&xi);
        kb=sizeof(float);
        if (ka>=0) ka=memoryrec((char *)&xf,source,ka,kb,length);
        kb=sizeof(float);
        if (ka>=0) ka=memoryrec((char *)&xg,source,ka,kb,length);
        kb=sizeof(float);
        if (ka>=0) ka=memoryrec((char *)&xh,source,ka,kb,length);
        if (ka<0) error=901;
        if (error==0) error=pdintensity->loadpoint(j,xa,xb,xc,xd,
          xe,xf,xg,xh,xi);
      }
      if (ka>=0) ka+=(npoint-mb)*3*sizeof(float);
      if ((error==0)&&(ka>=0)&&(tusage))
      { ma=(int)tproc[0];
        if ((ma>=0)&&(ma<numpe))
        { for (i=0;i<10;++i) tusage[ma*10+i]=tproc[i];
        }
        else error=901;
      }
    }
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Mscdrun::pasupimport

int Mscdrun::sendsup(int order)
{ int i,sendsize;
  char *sendbuf;

  sendbuf=NULL;
  if ((numpe>1)&&(((mype==0)&&(order<2))||((mype>0)&&(order==2))))
  { if (tusage) tusage[mype*10+1]+=jobtime->lastprocess();
    if (error==0)
    { sendsize=getsuplength(order); sendbuf=new char [sendsize];
      if (!sendbuf) error=102;
    }
    else sendsize=0;
    if (error==0) error=pasupexport(sendbuf,sendsize,order);
    if ((mype==0)&&(order<2))
    { for (i=1;i<numpe;++i)
      { if (error==0)
          error=mpisend((char *)&sendsize,sizeof(int),i,5120+i);
        else mpisend((char *)&sendsize,sizeof(int),i,5120+i);
        if ((error==0)&&(order==1))
          error=mpisend(sendbuf,sendsize,i,6144+i);
      }
    }
    else if ((mype>0)&&(order==2))
    { if (error==0)
        error=mpisend((char *)&sendsize,sizeof(int),0,7168+mype);
      else mpisend((char *)&sendsize,sizeof(int),0,7168+mype);
      if (error==0)
        error=mpisend(sendbuf,sendsize,0,8192+mype);
    }
    if (tusage) tusage[mype*10+2]+=jobtime->lastprocess();
  }
  if (sendbuf) delete [] sendbuf;

  return(error);
} //end of Mscdrun::sendsup

int Mscdrun::receivesup(int order)
{ int i,recsize,precsize;
  float xa,xb;
  char *recbuf;

  recbuf=NULL; recsize=precsize=0;
  if ((numpe>1)&&(((mype>0)&&(order==1))||((mype==0)&&(order==2))))
  { if (tusage) tusage[mype*10+1]+=jobtime->lastprocess();
    if (mype>0)
    { if (error==0)
        error=mpireceive((char *)&recsize,sizeof(int),0,5120+mype);
      else mpireceive((char *)&recsize,sizeof(int),0,5120+mype);
      if ((error==0)&&(recsize>sizeof(int)))
      { recbuf=new char [recsize];
        if (!recbuf) error=102;
        if (error==0) error=mpireceive(recbuf,recsize,0,6144+mype);
        if (error==0) error=pasupimport(recbuf,recsize,order);
      }
      else if (error==0) error=753;
    }
    else
    { precsize=recsize=0;
      for (i=1;i<numpe;++i)
      { if (error==0)
          error=mpireceive((char *)&recsize,sizeof(int),i,7168+i);
        else mpireceive((char *)&recsize,sizeof(int),i,7168+i);
        if ((error==0)&&(recsize>sizeof(int)))
        { if (i==1)
          { precsize=recsize; recbuf=new char [recsize];
          }
          else if (precsize<recsize)
          { if (recbuf) delete [] recbuf;
            precsize=recsize; recbuf=new char [recsize];
          }
          if (!recbuf) error=102;
          if (error==0) error=mpireceive(recbuf,recsize,i,8192+i);
          if (error==0) error=pasupimport(recbuf,recsize,order);
        }
        else if (error==0) error=753;
        xa=jobtime->process();
        xb=tusage[i*10+1]+tusage[i*10+2]+tusage[i*10+3]+tusage[i*10+4];
        if (xa>xb) tusage[i*10+2]+=(xa-xb);
      }
    }
    if (tusage) tusage[mype*10+3]+=jobtime->lastprocess();
  }
  if (recbuf) delete [] recbuf;

  return(error);
} //end of Mscdrun::receivesup

int Mscdrun::assistant(int order)
{ int i,j,k,m,afitmath;
  float xa,xb;
  float *afit,*xdata,*ydata,*ymod;

  afitmath=0; afit=xdata=ydata=ymod=NULL;
  if (error==0)
  { if ((npoint<2)||((numpe<2)&&(mype==0)))
    { pdbeg=0; pdend=npoint;
    }
    else if (numpe<2) pdbeg=pdend=npoint;
    else
    { k=npoint/numpe; m=npoint%numpe;
      if (mype<m)
      { pdbeg=(k+1)*mype; pdend=pdbeg+k+1;
      }
      else
      { pdbeg=(k+1)*m+k*(mype-m); pdend=pdbeg+k;
      }
    }
    if (pdbeg<0) pdbeg=0;
    if (pdend>npoint) pdend=npoint;
  }
  if ((error==0)&&(numpe>1)&&(mype==0)&&(order==0))
    error=sendsup(0);
  else if ((error==0)&&(numpe>1)&&(mype>0))
  { while (error==0) error=intensity(afitmath,afit,xdata,ydata,ymod);
    if (error==753) error=0;
  }
  if ((error==0)&&(mype==0)&&(order==0))
  { if (tusage) tusage[mype*10+1]+=jobtime->lastprocess();
    if (tusage)
    { xa=jobtime->process();
      if (xa>1.0e-3)
      { for (i=0;i<numpe;++i)
        { xb=0.0f; for (j=1;j<5;++j) xb+=tusage[i*10+j];
          if (xa>xb) tusage[i*10+4]+=(xa-xb);
        }
        for (i=0;i<numpe;++i)
        { xb=0.0f; for (j=1;j<5;++j) xb+=tusage[i*10+j];
          for (j=1;j<5;++j)
          { if (xb>1.0e-3) tusage[i*10+j+4]=tusage[i*10+j]/xb*100.0f;
            else tusage[i*10+j+4]=0.0f;
          }
        }
      }
      xa=0.0f;
      for (i=0;i<numpe;++i) for (j=1;j<5;++j) xa+=tusage[i*10+j];
      if (xa>1.0e-3)
      { tusage[numpe*10]=xa;
        for (j=1;j<5;++j)
        { xb=0.0f; for (i=0;i<numpe;++i) xb+=tusage[i*10+j];
          tusage[numpe*10+j]=xb; tusage[numpe*10+j+4]=xb/xa*100.0f;
        }
      }
    }
  }
  return(error);
} //end of Mscdrun::assistant

