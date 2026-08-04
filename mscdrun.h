/*
----------------------------------------------------------------------
  C++ classes designed for photoelectron diffraction software package
  Yufeng Chen, Michel A. Van Hove
  Lawrence Berkeley National Laboratory (LBNL), Berkeley, CA 94720
  Copyright (c) Van Hove Group 1997-1998. All rights reserved.
----------------------------------------------------------------------

  Yufeng Chen, LBNL MS 2-100, 1 Cyclotron Road, Berkeley, CA 94720
  Phone (510) 486-4581 Fax (510) 486-5530 Email. ychen@LBL.gov

  M. A. Van Hove, LBNL MS 66-100, 1 Cyclotron Road, Berkeley, CA 94720
  Phone (510) 486-6160 Fax (510) 486-4995 Email. MAVanhove@LBL.gov
----------------------------------------------------------------------

List of classes

1. Mscdrun
Description: execute a job request
Constructors and Destructor:
  Mscdrun(char *filein,int jobnum,int jobtotal,int dispmode,
    int displog) --- constructor, with input file name
      filein, for the number jobnum of the
      total number of jobtotal, with display mode dispmode, and if
      displog nonzero, save a log file
Member Functions:
  float getmemory() --- get memory it took in unit of bytes
  int readparameter() --- read input file
  int fitcheck() --- check fitting parameters and adjust scanning range
  int paradisp() --- display parameters
  int symtrivert() --- symmetry analysis for triple vertex scattering
    events
  int symdblvert() --- symmetry analysis for double vertex scattering
    events up to the detector


List of functions
1.  float fitpdintensity(int fitmath,int ndata,int nfit,int *yafit,
      float *afit,float *safit,float *xdata,float *ydata,float *ymod,
      float *dyda,Mscdrun *mscdrun) --- photoemission intensity function
        for fitting


----------------------------------------------------------------------
*/

#ifndef __MSCDRUN_H
#define __MSCDRUN_H

#include <fstream>

#include "userutil.h"
#include "fcomplex.h"
#include "phase.h"
#include "radmat.h"
#include "meanpath.h"
#include "vibrate.h"
#include "rotamat.h"
#include "msfuncs.h"
#include "jobtime.h"
#include "pdinten.h"

class Mscdrun
{ private:
    int datatype,scanmode,jobnum,jobtotal,dispmode,displog,linitial,
      natoms,katoms,eatoms,lnum,nfit,nfollow,ncurve,npoint,nkout,
      ndtheta,ndphi,ndangle,trynum,filenum,totlayer,
      msorder,raorder,finals,fitmode,fitmath,trymax,devstep,ntrieven,
      nscorse,ndbleven,nalpha,radim,ntrielem,nsymm,mscatter,
      beampol,mype,numpe,pdbeg,pdend,sizeint,error,ATA;
    float ftolerance,kmin,kmax,kstep,dtmin,dtmax,dtstep,dpmin,
      dpmax,dpstep,radius,depth,ltheta,lphi,mtheta,mphi,accepang,
      lattice,valence,bandgap,density,mweight,nearest,biggest,vinner,
      tdebye,tsample,pathcut,basemem,ATAweight1,ATAweight2;
    char *username,*sysname,*finname,*foutname,*rmfile,*psfile;
    int *lakatom,*layatom,*pdnum,*tevenadd,*tevendim,*tevencut,
      *devenadd;
    float *laycell,*layorig,*layfit,*laxcell,*laxorig,*fitvars,
      *aweight,*magamp,*patom,*tevenpar,*devenpar,*talpha,*tgamma,
      *fithist,*tusage;
    Textout *flogout;
    Fcomplex *tevenelem,*devenelem,*devendetec;
    Phaseshift *phaseshift;
    Radialmatrix *radmatrix;
    Meanpath *meanpath;
    Vibration *vibrate;
    Rotamat *evenmat,*termmat;
    Expix *expix;
    Hankel *hanka,*hankb;
    Pdintensity *pdintensity;
    Jobtime *jobtime;
  public:
    Mscdrun(int imype,int inumpe);
    ~Mscdrun();
    void init();
    int loadvars(char *ifinname,char *ifoutname,int ijobnum,
      int ijobtotal,int idispmode,int idisplog,Textout *iflogout);
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int sendjobs();
    int receivejobs();
    int getsuplength(int order);
    int pasupexport(char *dest,int length=0,int order=0);
    int pasupimport(char *source,int length=0,int order=0);
    int sendsup(int order);
    int receivesup(int order);
    int assistant(int order);
    int getnumpoint();
    int getnumcurve();
    int getnumfit();
    float getmemory();
    int readparameter();
    int paradisp();
    int paradisplog();
    int symtrivert();
    int symdblvert();
    int precutable();
    int intensity(int afitmath,float *afit,float *xdata,float *ydata,
      float *ymod);
    float fitphotoemission(int afitmath,int ndata,
      int nafit,int *yafit,float *afit,float *safit,float *xdata,
      float *ydata,float *ymod,float *dyda);
    int getfitdata(float *xdata,float *ydata,float *afit,float *safit);
    int fitarpefs(Mscdrun *mscdrun);
    int savecurve(char *filename,char *usermessage);
    int sendobjects();
    int receiveobjects();
    int sendpoints();
    int receivepoints();
  private:
    int fitcheck();
    Fcomplex evenelem(int akind,int ma,int na,int mb,int nb,
      float akin,float vka,float vkb,float beta);
    Fcomplex ATAevenelem(int akind,int ma,int na,int mb,int nb,
      float akin,float vka,float vkb,float beta, float ATAweight1,
      float ATAweight2);
    Fcomplex evenbelem(int alf,int ma,int na,int mb,int nb,
      float vkb,float beta);
    int alltrievent(int forcut,float akin);
    int onerotation(float *patoma,float *patomb,float *patomc,
      float *alpha,float *beta,float *gamma);
    int alldblevent(float akin,float *xdetec);
    int allevendetec(float akin,float *xdetec);
    int onevenemit(int ia,int ib,int alf,int am,float akin,
      float *xdetec,float *polaron,Fcomplex *aemitelem);
    Fcomplex onemidetec(float akin,int ie,int alf,int am,
      float *xdetec,float *polaron);
    Fcomplex matrixelement(int ali,int alf,int am,float akin);
    int allrotation();
    int summation(float akin,float *xdetec,float *polaron,
      float *suminten,float *bakinten,Fcomplex *asum,Fcomplex *bsum,
      Fcomplex *csum);
    int dispintensity(int messagenum,int pointnum=0,int afitmath=0,
      float areliable=0.0,float breliable=0.0,float akout=0.0,
      float adtheta=0.0,float adphi=0.0,float suminten=0.0);
    int makeatoms(float disturb);
    int maketripar();
    int makedblpar();
    float kinside(float akout,float avinner);
    float thetainside(int accepnum,float akin,float akout,
      float thetaout,float phiout,float altheta,float alphi,
      float *xdetec,float *polaron);
    int sendarray(char *sendbuf,char *source,int ka,int kb,
      int sendsize);
    int recarray(char *dest,char *recbuf,int ka,int kb,
      int recsize);
    int recmemnew(char *dest,char *recbuf,int ka,int kb,
      int recsize);
};

float fitpdintensity(int fitmath,int ndata,int nfit,int *yafit,
  float *afit,float *safit,float *xdata,float *ydata,float *ymod,
  float *dyda,Mscdrun *mscdrun);

#endif //__MSCDRUN_H

/*
----------------------------------------------------------------------
dimmensions
  psfile[100*4] lakatom[101*4] layatom[101*4] pdnum[npoint]
  tevenadd[natoms*natoms*natoms] --- address of triple vertex event
  tevendim[natoms*natoms*natoms] --- dimension of the triple vertex
    event
  tevencut[msorder*natoms*natoms] --- cut table for each order
  tevenelem[ntrielem] --- matrix element of triple vertex event
  talpha[natoms*natoms*natoms] tgamma[natoms*natoms*natoms]
  devenadd[natoms*natoms] --- address of unique double vertex events
  devenelem[ndbleven*radim] --- array element of double vertex event
  devendetec[natoms*natoms*radim] --- array element up to detector with
    consideration of phase correction with respect to the origin
  fitvars[4*20] laycell[101*4] layorig[101*4] layfit[101*4]
  laxcell[101*4] laxorig[101*4]
  aweight[4] patoms[300*12]
  tevenpar[ntrieven*10] devenpar[ndbleven*7]
  tusage[(numpe+1)*10] --- processor time distribution on computation,
    sending, receiving and idle
  phaseshift[katoms]


lakatom   0          1       2          3
          atom-kind  emiter  atomunits  layer-number

layfit    0             1             2              3
          lattice       ---           vinner         tdebye
          origin-vector ---           layer-spacing  scaling-factor

patom     0  1  2  3  4  5  6        7        8  9  10  11  12
          x  y  z  i  j  k  kindatom emitter

finals
    0     real calculation for both channels
    1     high channel (li+1) only
    2     low channel (li-1) only
    3     reference only
    4     scattering term only
    5     same as 0

layfit    0             1             2              3
          lattice       ---           vinner         tdebye
          origin-vector ---           layer-spacing  scaling-factor

layorig   0             1             2              3
          lattice       ---           vinner         tdebye
          origin-vector origin-angle  layer-spacing  scaling-factor

tevenpar   0        1            2             3       4
          lengtha  inverse_lena  inverse_lengb cosbeta atom_kind

          5        6       7     8     9
          eledim   memadd  ia    ib   ic

devenpar  0   1   2   3        4          5   6
          xa  ya  za  lengtha  atom_kind  ia  ib

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

tusage
    0     processor id
    1     computation time
    2     sending time
    3     receiving time
    4     idle time
    5     computation time percent
    6     sending time percent
    7     receiving time percent
    8     idle time percent
    9     reserved
----------------------------------------------------------------------
*/
