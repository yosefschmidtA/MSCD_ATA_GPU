#include <iostream>
#include <time.h>

#include "userutil.h"
#include "jobtime.h"

Jobtime::Jobtime()
{ reset();
} //end of Jobtime::Jobtime

void Jobtime::reset()
{ lastvisit=0;
  time(&timebeg); timemid=timeend=timebeg;
  tcpubeg=timeprocessor(); tcpumid=tcpuend=tcpubeg;
  return;
} //end of Jobtime::reset

void Jobtime::stop()
{ if (lastvisit!=255)
  { lastvisit=255;
    timemid=timeend; time(&timeend);
    tcpumid=tcpuend; tcpuend=timeprocessor();
  }
  return;
} //end of Jobtime::stop

void Jobtime::resume()
{ long ka;
  if (lastvisit==255)
  { lastvisit=0;
    ka=timeprocessor();
    ka-=tcpuend; tcpubeg+=ka; tcpumid+=ka; tcpuend+=ka;
  }
  return;
} //end of Jobtime::resume

float Jobtime::lastelapse()
{ float xa;
  if ((lastvisit==0)||(lastvisit==1))
  { lastvisit=1;
    timemid=timeend; time(&timeend);
    tcpumid=tcpuend; tcpuend=timeprocessor();
  }
  else if (lastvisit!=255) lastvisit=0;
  xa=(float)(timeend-timemid);
  return(xa);
} //end of Jobtime::lastelapse

float Jobtime::lastprocess()
{ float xa;
  if ((lastvisit==0)||(lastvisit==2))
  { lastvisit=2;
    timemid=timeend; time(&timeend);
    tcpumid=tcpuend; tcpuend=timeprocessor();
  }
  else if (lastvisit!=255) lastvisit=0;
  xa=(float)(tcpuend-tcpumid);
  return(xa);
} //end of Jobtime::lastprocess

float Jobtime::elapse()
{ float xa;
  if (lastvisit!=255)
  { timemid=timeend; time(&timeend);
    tcpumid=tcpuend; tcpuend=timeprocessor();
  }
  xa=(float)(timeend-timebeg);
  return(xa);
} //end of Jobtime::elapse

float Jobtime::process()
{ float xa;
  if (lastvisit!=255)
  { timemid=timeend; time(&timeend);
    tcpumid=tcpuend; tcpuend=timeprocessor();
  }
  xa=(float)(tcpuend-tcpubeg);
  return(xa);
} //end of Jobtime::process

time_t Jobtime::begintime(char *begdate,int size)
{
  if ((begdate)&&(size>12))
  { stringcopy(begdate,ctime(&timebeg)+4,7);
    begdate[6]=',';
    stringcopy(begdate+7,ctime(&timebeg)+19,6);
    begdate[12]='\0';
  }

  return(timebeg);
} //end of Jobtime::getbegintime

time_t Jobtime::endtime(char *enddate,int size)
{
  if ((enddate)&&(size>12))
  { stringcopy(enddate,ctime(&timeend)+4,7);
    enddate[6]=',';
    stringcopy(enddate+7,ctime(&timeend)+19,6);
    enddate[12]='\0';
  }

  return(timeend);
} //end of Jobtime::getendtime

int Jobtime::timestamp(Textout& textout)
{ int k,numpe,error;
  float timecpu;
  char namebuf[80];

  error=textout.geterror();
  if (error==0)
  { k=computer(namebuf,25); numpe=mpigetnumpe();
    timecpu=process();
    if ((timecpu>0.0)&&(numpe>1)) timecpu*=(float)numpe;
    textout.string("This calculation ",0,16);
    if ((timecpu<0.0)||(timecpu>2.0e7)) textout.string("was made on ");
    else if (timecpu<300.0)
    { k=(int)(timecpu+0.5); textout.string("took ");
      textout.integer(k,4); textout.string(" processor seconds on ");
    }
    else if (timecpu<36000.0)
    { k=(int)(timecpu/60.0+0.5); textout.string("took ");
      textout.integer(k,4); textout.string(" processor minutes on ");
    }
    else if (timecpu<864000.0)
    { k=(int)(timecpu/3600.0+0.5); textout.string("took ");
      textout.integer(k,4); textout.string(" processor hours on ");
    }
    else
    { k=(int)(timecpu/8.64e4+0.5);  textout.string("took ");
      textout.integer(k,4); textout.string(" processor days on ");
    }
    textout.string(namebuf,0,1);
    textout.string("     starting on ");
    textout.string(ctime(&timebeg));
    textout.string("   and ending on ");
    textout.string(ctime(&timeend),0,1);
    if (error==0) error=textout.geterror();
  }

  return(error);
} //end of Jobtime::timestamp

