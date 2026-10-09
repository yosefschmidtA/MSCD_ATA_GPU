// This file is for Sun workstation only

/*
----------------------------------------------------------------------
code table for different computers

    computer_name         code    utility_filename

    Sequential computers (1-99)
      IBM-PC              11      userpc.cpp
      Macintosh           12      usermac.cpp
      Sun_workstation     13      usersun.cpp
      Cray J90            14      userj90.cpp
      Cray C90            15      userc90.cpp

    Parallel computers   (101-199)
      Cray T3E           111      usert3e.cpp
      COMPS network      112      usercomp.cpp
----------------------------------------------------------------------
*/

#include <iostream>
#include <time.h>
#include <sys/times.h>

#include "userutil.h"

#ifndef CLK_TCK
#ifndef HZ
#include <sys/param.h>
#endif
#endif
#ifndef HZ
#define HZ 60
#endif
#ifndef CLK_TCK
#define CLK_TCK HZ
#endif

int computer(char *name,int size)
{ int computercode;

  computercode=14;
  if ((name!=NULL)&&(size>10))
  { if (size>25) size=25;
    stringcopy(name,"a Cray J90 cluster",size);
  }

  return(computercode);
} //end of computer

long timeprocessor()
{ long k,ta;
  struct tms timebuf;

  k=CLK_TCK; if (k<1) k=1;
  ta=times(&timebuf);
  if (ta<0) ta=-1;
  else ta=(timebuf.tms_utime+timebuf.tms_stime)/k;
  return(ta);
} //end of timeprocessor

int mpiinit(int *argc,char ***argv)
{ int error;
  error=(*argv)[*argc-1][0]; error=0;
  return(error);
} //end of mpiinit

int mpigetmype()
{ int myrank;
  myrank=0;
  return(myrank);
} //end of mpigetmype

int mpigetnumpe()
{ int size;
  size=1;
  return(size);
} //end of mpigetnumpe

int mpisend(char *buffer,int bufsize,int destpe,int tag)
{ int error;
  error=buffer[bufsize-1]; error=destpe; error=tag; error=0;
  return(error);
} //end of mpisend

int mpireceive(char *buffer,int bufsize,int sourcepe,int tag)
{ int error;
  error=buffer[bufsize-1]; error=sourcepe; error=tag; error=0;
  error=0;
  return(error);
} //end of mpireceive

int mpiend()
{ return(0);
} //end of mpiend

