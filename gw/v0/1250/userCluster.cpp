// This file is for Cray supercomputer only

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
      Linux Clusters     113      userCluster.cpp
----------------------------------------------------------------------
*/

#include <iostream>
#include <ctime>
#include <sys/times.h>
#include <mpi.h>
//#include <mpi++.h>

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
{ int ka,kb,kc,kd,numpe,computercode;

  computercode=113;
  if ((name!=NULL)&&(size>10))
  { numpe=mpigetnumpe();
    if (size>25) size=25;
    ka=numpe%10; kb=(numpe/10)%10; kc=(numpe/100)%10; kd=numpe/1000;
    if (kd>0) name[0]=kd+'0'; else name[0]=' ';
    if ((kd>0)||(kc>0)) name[1]=kc+'0'; else name[1]=' ';
    if ((kd>0)||(kc>0)||(kb>0)) name[2]=kb+'0'; else name[2]=' ';
    name[3]=ka+'0';
    if (numpe<=1) stringcopy(name+4," Linux Cluster",size-4);
    else stringcopy(name+4," Linux Cluster",size-4);
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
{ MPI_Init(argc,argv);
  return(0);
} //end of mpiinit

int mpigetmype()
{ int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  return(myrank);
} //end of mpigetmype

int mpigetnumpe()
{ int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  return(size);
} //end of mpigetnumpe

int mpisend(char *buffer,int bufsize,int destpe,int tag)
{ int error;
  if (buffer)
  { error=MPI_Send(buffer,bufsize,MPI_CHAR,destpe,tag,MPI_COMM_WORLD);
    if (error!=0) error=751;
  }
  else error=103;
  return(error);
} //end of mpisend

int mpireceive(char *buffer,int bufsize,int sourcepe,int tag)
{ int error;
  MPI_Status status;

  if (buffer)
  { error=MPI_Recv(buffer,bufsize,MPI_CHAR,sourcepe,tag,
      MPI_COMM_WORLD,&status);
    if (error!=0) error=752;
  }
  else error=103;
  return(error);
} //end of mpireceive

int mpiend()
{ MPI_Finalize();
  return(0);
} //end of mpiend

