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

1. Mscdjob
Description: management of jobs
Constructors and Destructor:
  Mscdjob() --- constructor.
Member Functions:
  float getmemory() --- get memory it took in unit of bytes
  int getdisplog() --- return non zero if log file need to be created
  int loadflog(Textout *iflogout) --- load log file stream
  int takejob(char *filename) --- add a job in job list, with the input
    file name as filename
  int listjob(char *filename,char *usermessage) --- save the job list
    with their input filename, output filename, or error message
  int execute() --- execute the list of jobs
----------------------------------------------------------------------
*/

#ifndef __MSCDJOB_H
#define __MSCDJOB_H

#include "userutil.h"

class Mscdjob
{ private:
    int datatype,jobtotal,dispmode,displog,mype,numpe,error;
    float basemem;
    int *errorcode;
    char *finname,*foutname;
    Textout *flogout;
  public:
    Mscdjob(int imype,int inumpe);
    ~Mscdjob();
    void init();
    int getlength();
    int paexport(char *dest,int length=0);
    int paimport(char *source,int length=0);
    int sendjobs();
    int receivejobs();
    float getmemory();
    int getdisplog();
    int loadflog(Textout *iflogout);
    int takejob(char *filename);
    int listjob(char *filename,char *usermessage);
    int execute(char *usermessage);
};

#endif //__MSCDJOB_H

