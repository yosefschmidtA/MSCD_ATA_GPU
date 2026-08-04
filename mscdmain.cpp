#include <iostream>
#include <fstream>

#include "userinfo.h"
#include "userutil.h"
#include "mscdjob.h"

int main(int argc,char **argv)
{ int i,displog,mype,numpe,machine,error;
  Textout *flogout;

  error=0; flogout=NULL;
  mpiinit(&argc,&argv);
  mype=mpigetmype(); numpe=mpigetnumpe();
  machine=computer();
  Userinfo userinfo(mype,numpe); userinfo.init();
  userinfo.loadparameter(1.37f,"MSCD",
    "Yufeng Chen and Michel A Van Hove",
    "Calculation of photoelectron diffraction and dichroism",machine);
  if (mype==0)
  { Textout conout;
    conout.string(userinfo.usermsg,0,1);
    if (userinfo.version==0.0) error=102;
  }
  Mscdjob mscdjob(mype,numpe);
  if (mype==0)
  { mscdjob.init();
    if ((error==0)&&(argc<2)) error=mscdjob.takejob("mscdin.txt");
    else
    { for (i=1;(error==0)&&(i<argc);++i) error=mscdjob.takejob(argv[i]);
    }
    displog=mscdjob.getdisplog();
    if ((error==0)&&(displog>0))
    { flogout=new Textout("mscdlist.txt");
      if (flogout) error=flogout->geterror();
      else error=102;
      if (error==0) flogout->string(userinfo.usermsg,0,1);
    }
    if (error==0) error=mscdjob.loadflog(flogout);
  }
  if ((error==0)&&(numpe>1)&&(mype==0)) error=mscdjob.sendjobs();
  else if ((numpe>1)&&(mype==0)) mscdjob.sendjobs();
  else if (numpe>1) error=mscdjob.receivejobs();
  if (error==0) error=mscdjob.execute(userinfo.usermsg);
  if (mype==0)
  { if (error==0) error=mscdjob.listjob("mscdout.txt",userinfo.usermsg);
    if (flogout) delete flogout;
    Textout conout;
    Errorinfo errorinfo(error);
    conout.string(errorinfo.errmsg,0,16+1);
  }
  mpiend();

  return(error);
} //end of main


/*
----------------------------------------------------------------------
MSCD version history


MSCD package began the programming on December 9, 1996, using ANSI C++
object oriented programming (OOP) technique and massage passing
interfacing (MPI) parallel communication.

Version 1.37, June 10, 1998
  Make the source code compatible to Macintosh CodeWarrier Compiler.
  Fixed a bug in net search mode.

Version 1.36, June 1, 1998
  Impoved emission angle search feature.

Version 1.35, April 21, 1998
  Number of emitters in different logical layers is considered.
  Total number of data points can be up to 30000.

Version 1.34, February 22, 1998
  Fixed a bug in calculation of large cluster.

Version 1.32, January 16, 1998
  Fixed a bug in calculation of reliability factors for two hologram
  data files.
  Calculate Van Hove and Fadley group's reliability factors in
  program caldif

Version 1.31, August 31, 1997
  Expand the Rehr-Albers approximation order from 3 (10x10 matrix) to
  4 (15x15 matrix). To choose 4th order of Rehr-Albers approximation,
  set raorder=4 in the input file.

Version 1.30, August 30, 1997
  Expand the Rehr-Albers approximation order from 2 (6x6 matrix) to
  3 (10x10 matrix). To choose 3rd order of R-A approximation,
  set raorder=3 in the input file.

Version 1.24, June 30, 1997
  1) Add usermac.cpp file. The program works on Macintosh computer.
  2) Add two functions in userutil.cpp, skipendline and getendline,
     making the program package supports PC, Unix and Macintosh
     ascii text format data file.
  3) Add usercomp.cpp file. The program works on COMPS network of
     workstation in parallel mode.
  4) Add cpu time analysis function for each processor in parallel
     mode. The processor time are divided into four categories,
     computation, sending, receiving and idle. To activate this
     feature, set display mode to 10. The total processor time
     and its distribution are shown in output file mscdlist.txt.

Version 1.23, June 12, 1997
  The maxmum number of quantum momenta lmax is expanded again to 60.
  To make it possible, the makecurve function in rotamat.cpp has been
  changed, using double precission floating point numbers.

Version 1.22, June 10, 1997
  The maximum energy and number of quantum momenta are expanded.
  The maximum kmax set to 25.0, corresponding energy 2.38 keV.
  The maximum lmax set to 30.

Version 1.21, June 3, 1997
  Removed a bug in fitting mode and a possible bug in test-running mode.

Version 1.20, May 28, 1997
  First parallel version, working on Cray T3E in MPP parallel mode.

Version 1.10, April 1, 1997
  Introduce the rotation transfer method to save computation time for
  circular polarized light and magnetic circular dichroism.

Version 1.00, March 1, 1997
  First version, working on PC, Sun workstation, Cray J90
  supercomputer, and Cray T3E supercomputer in sequential mode,
  covering all the features of the SCAT Version 3.65 package.
----------------------------------------------------------------------
*/

