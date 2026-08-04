// This file is for IBM PC and compatible only

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

#include "userutil.h"

int computer(char *name,int size)
{ int computercode;

  computercode=12;
  if ((name!=NULL)&&(size>10))
  { if (size>25) size=25;
    stringcopy(name,"a Macintosh computer",size);
  }

  return(computercode);
} //end of computer

long timeprocessor()
{ long ta;

  time_t timebuf;
  ta=time(&timebuf);
  if (ta<-10000)
  { ta<<=1; ta>>=1;
  }
  if (ta<0) ta=-1;
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

