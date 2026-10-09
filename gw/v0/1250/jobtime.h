/*
----------------------------------------------------------------------
  C++ classes designed for relaxation energy software package
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

1. Jobtime
Description: calculation of elapse or processor time of a job
Constructors and Destructor:
  Jobtime() --- constructor.
Member Functions:
  void reset() --- reset job time count
  void stop() --- stop job time count
  void resume() --- resume job time count
  float lastelapse() --- return the elapse time from last visit
  float lastprocess() --- return the processor time from last visit
  float elapse() --- return the total elapse time from the beginning
  float process() --- return the total processor time from beginning
  time_t begintime() --- return the begining time
  time_t endtime() --- return the ending time
  int timestamp(Textout& textout) --- put a time stamp to the
      output stream textout
----------------------------------------------------------------------
*/

#ifndef __JOBTIME_H
#define __JOBTIME_H

#include <time.h>

#include "userutil.h"

class Jobtime
{ private:
    int lastvisit;
    long tcpubeg,tcpumid,tcpuend;
    time_t timebeg,timemid,timeend;
  public:
    Jobtime();
    void reset();
    void stop();
    void resume();
    float lastelapse();
    float lastprocess();
    float elapse();
    float process();
    time_t begintime(char *begdate=NULL,int size=0);
    time_t endtime(char *enddate=NULL,int size=0);
    int timestamp(Textout& textout);
};

#endif //__JOBTIME_H

