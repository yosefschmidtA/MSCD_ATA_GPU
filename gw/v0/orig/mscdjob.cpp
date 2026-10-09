#include <iostream>
#include <fstream>
#include <iomanip>

#include "userutil.h"
#include "userinfo.h"
#include "mscdrun.h"
#include "mscdjob.h"

Mscdjob::Mscdjob(int imype,int inumpe)
{ mype=imype; numpe=inumpe;
  error=107; init();
} //end of Mscdjob::Mscdjob

Mscdjob::~Mscdjob()
{ if (finname) delete [] finname; if (foutname) delete [] foutname;
  if (errorcode) delete [] errorcode;
} //end of Mscdjob::~Mscdjob

void Mscdjob::init()
{ if ((error==107)||(error==103))
  { jobtotal=dispmode=displog=0; datatype=741; basemem=0.0f;
    flogout=NULL; finname=foutname=NULL; errorcode=NULL;
  }
  if (error==103)
  { errorcode=new int [100];
    finname=new char [100*100]; foutname=new char [100*100];
    if ((!finname)||(!foutname)||(!errorcode)) error=102;
  }
  else if (error==106)
  { if (errorcode)
    { errorcode=new int [100]; if (!errorcode) error=102;
    }
    if (finname)
    { finname=new char [100*100]; if (!finname) error=102;
    }
    if (foutname)
    { foutname=new char [100*100]; if (!foutname) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Mscdjob::init

int Mscdjob::getlength()
{ int ka;
  if (error==0) ka=sizeof(Mscdjob)+100*sizeof(int)+100*100*2+
    sizeof(int);
  else ka=0;
  return(ka);
} //end of Mscdjob::getlength

int Mscdjob::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Mscdjob);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=100*sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)errorcode,ka,kb,length);
    kb=100*100;
    if (ka>=0) ka=memorysend(dest,finname,ka,kb,length);
    kb=100*100;
    if (ka>=0) ka=memorysend(dest,foutname,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Mscdjob::paexport

int Mscdjob::paimport(char *source,int length)
{ int ka,kb,ma,mb;

  ma=mype; mb=numpe;
  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Mscdjob);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=100*sizeof(int);
    if (ka>=0) ka=memoryrec((char *)errorcode,source,ka,kb,length);
    kb=100*100;
    if (ka>=0) ka=memoryrec(finname,source,ka,kb,length);
    kb=100*100;
    if (ka>=0) ka=memoryrec(foutname,source,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  mype=ma; numpe=mb;
  return(error);
} //end of Mscdjob::paimport

int Mscdjob::sendjobs()
{ int i,sendsize;
  char *sendbuf;

  sendbuf=NULL;
  if ((numpe>1)&&(mype==0))
  { if (error==0)
    { sendsize=getlength(); sendbuf=new char [sendsize];
      if (!sendbuf) error=102;
    }
    else sendsize=0;
    if (error==0) error=paexport(sendbuf,sendsize);
    for (i=1;i<numpe;++i)
    { if (error==0)
        error=mpisend((char *)&sendsize,sizeof(int),i,1024+i);
      else mpisend((char *)&sendsize,sizeof(int),i,1024+i);
      if (error==0) error=mpisend(sendbuf,sendsize,i,2048+i);
    }
  }
  if (sendbuf) delete [] sendbuf;

  return(error);
} //end of Mscdjob::sendjobs

int Mscdjob::receivejobs()
{ int recsize;
  char *recbuf;

  recbuf=NULL;
  if ((numpe>1)&&(mype>0))
  { if (error==103)
      error=mpireceive((char *)&recsize,sizeof(int),0,1024+mype);
    else mpireceive((char *)&recsize,sizeof(int),0,1024+mype);
    if ((error==0)&&(recsize>sizeof(int)))
    { recbuf=new char [recsize];
      if (!recbuf) error=102;
      if (error==0) error=mpireceive(recbuf,recsize,0,2048+mype);
      if (error==0)
      { error=103; error=paimport(recbuf,recsize);
      }
    }
    else if (error==0) error=753;
    dispmode=1; displog=0; flogout=NULL;
  }
  if (recbuf) delete [] recbuf;

  return(error);
} //end of Mscdjob::receivejobs

float Mscdjob::getmemory()
{ float xa;
  if (error==0) xa=(float)(basemem+sizeof(Mscdjob)+20512*sizeof(char)+
    100*sizeof(int));
  else xa=-1.0f;
  return(xa);
} //end of Mscdjob::getmemory

int Mscdjob::takejob(char *filename)
{ int i,j,k,m,type,begrow,linenum,dispjob,memsize;
  float xa,xb,xc,tolerance;
  char *membuf,*filebuf;

  membuf=filebuf=NULL; memsize=256;
  if ((error==101)||(error==0))
  { membuf=new char [memsize];
    filebuf=new char [memsize];
    if ((!membuf)&&(!filebuf)) error=102;
    else error=0;
  }
  if ((error==0)&&(jobtotal<100-1))
  { spaceclear(filename);
    k=1;
    for (j=0;j<jobtotal;++j)
    { k=stringcomp(filename,finname+j*100,100);
      if (k==0) j=jobtotal;
    }
    if (k!=0)
    { std::ifstream fin(filename,std::ios::in);
      if (!fin) error=201+(datatype<<16);
      else
      { fin >> xa >> xb >> xc;
        type=(int)xa; begrow=(int)xb; linenum=(int)xc;
        if ((type==741)||(type==751))
          for (i=1;i<begrow;++i) skipendline(fin);
        if (!fin) error=207+(datatype<<16);
        else if (type==741)
        { for (i=0;i<linenum;++i)
          { getendline(fin,membuf,memsize);
            spacehide(membuf);
            j=nextnonwhite(membuf);
            k=nextwhite(membuf,j);
            stringcopy(filebuf,membuf+j,k-j+1);
            spaceclear(filebuf);
            if (((filebuf[0]=='p')||(filebuf[0]=='P'))&&
               ((filebuf[1]=='d')||(filebuf[1]=='D')||
               (filebuf[1]=='e')||(filebuf[1]=='E')))
            { j=nextnonwhite(membuf,k);
              k=nextwhite(membuf,j);
              stringcopy(filebuf,membuf+j,k-j+1);
              spaceback(filebuf);
              spaceclear(filebuf);
              for (j=0;j<jobtotal;++j)
              { k=stringcomp(filebuf,foutname+j*100,100);
                if (k==0) j=jobtotal;
              }
              if (k!=0)
              { fin >> xa >> xb >> xc;
                dispjob=(int)xb; tolerance=xc;
                if (!fin) error=207+(datatype<<16);
              }
              else
              { dispjob=0; tolerance=0.0f;
              }
              if ((error==0)&&(k!=0))
              { stringcopy(finname+jobtotal*100,filename,100);
                stringcopy(foutname+jobtotal*100,filebuf,100);
                ++jobtotal;

                if (jobtotal==1) displog=dispjob/10;
                else if (dispjob<10) displog=0;
                if (tolerance<1.0e-5) tolerance=0.0f;
                dispjob%=10; m=computer();
                if ((jobtotal==1)&&(m!=43)&&(m!=50)) dispmode=dispjob;
                else if ((jobtotal==1)&&(dispjob<4)) dispmode=dispjob;
                else if ((jobtotal==1)&&(dispjob==4)&&(tolerance>0.0))
                  dispmode=dispjob;
                else if (jobtotal==1) dispmode=0;
                else if (((dispmode==0)||(dispmode>3))&&
                  ((dispjob==0)||(dispjob>3)))
                  dispmode=0;
                else if (((dispmode==0)||(dispmode>3))&&
                  (dispjob==3)&&(tolerance>0.0))
                  dispmode=0;
                else if (dispmode==0) dispmode=dispjob;
                else if ((dispjob!=0)&&(dispmode>dispjob))
                  dispmode=dispjob;
              }
              i=linenum+1;
            }
          }
          if (i==linenum) error=207+(datatype<<16);
        }
        else if (type==751)
        { for (i=0;(error==0)&&(i<linenum);++i)
          { getendline(fin,membuf,memsize);
            spacehide(membuf);
            j=nextnonwhite(membuf);
            k=nextwhite(membuf,j);
            stringcopy(filebuf,membuf+j,k-j+1);
            spaceclear(filebuf);
            if ((filebuf[0]=='i')||(filebuf[0]=='I'))
            { j=nextnonwhite(membuf,k);
              k=nextwhite(membuf,j);
              stringcopy(filebuf,membuf+j,k-j+1);
              spaceback(filebuf);
              spaceclear(filebuf);
              k=stringcomp(filebuf,filename,100);
              if (k!=0) error=takejob(filebuf);
            }
          }
        }
        else error=207+(datatype<<16);
      }
    }
  }
  if (membuf) delete [] membuf; if (filebuf) delete [] filebuf;
  return(error);
} //end of Mscdjob:takejob

int Mscdjob::getdisplog()
{ if (error!=0) displog=0;
  return(displog);
} //end of Mscdjob::getdisplog

int Mscdjob::loadflog(Textout *iflogout)
{ if ((error==0)&&(displog>0)&&(iflogout)) flogout=iflogout;
  else if ((error==0)&&(displog>0))
  { displog=0; flogout=0; error=901;
  }
  else
  { displog=0; flogout=0;
  }
  return(error);
} //end of Mscdjob::loadlogfd

int Mscdjob::listjob(char *filename,char *usermessage)
{ int i,begrow,linenum;

  linenum=jobtotal*3+1;
  if ((error==0)&&((dispmode==0)||(dispmode>2)))
  { Textout conout;
    conout.charfill('_',37,16+1);
    conout.string(" Job report",0,1);
    for (i=0;i<jobtotal;++i)
    { conout.integer(i+1,4,16); conout.space(3);
      conout.string((finname+i*100),0,1);
      if (errorcode[i]==0)
      { conout.string("       result saved in file ");
        conout.string((foutname+i*100),0,1);
      }
      else
      { Errorinfo errorinfo(errorcode[i]);
        conout.space(7); conout.string(errorinfo.errmsg,0,1);
      }
    }
    conout.charfill('_',37,16+1);
    conout.string(" Job report saved in file mscdout.txt",0,1);
    if (displog>0)
      conout.string(" Messages saved in file mscdlist.txt",0,1);
  }

  if ((error==0)&&(displog>0)&&(flogout))
  { flogout->charfill('_',37,16+1);
    flogout->string(" Job report",0,1);
    for (i=0;i<jobtotal;++i)
    { flogout->integer(i+1,4,16); flogout->space(3);
      flogout->string((finname+i*100),0,1);
      if (errorcode[i]==0)
      { flogout->string("       result saved in file ");
        flogout->string((foutname+i*100),0,1);
      }
      else
      { Errorinfo errorinfo(errorcode[i]);
        flogout->space(7); flogout->string(errorinfo.errmsg,0,1);
      }
    }
    flogout->charfill('_',37,16+1);
    flogout->string(" Job report saved in file mscdout.txt",0,1);
    if (displog>0)
      flogout->string(" Messages saved in file mscdlist.txt",0,1);
  }

  if (error==0)
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { i=0; begrow=5;
      while (usermessage[i]!='\0')
      { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
        ++i;
      }
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(linenum,6); fileout.space(5);
      fileout.string("datakind beginning-row linenumbers",0,1);
      fileout.string(usermessage,0,1);

      fileout.string(" Job report",0,1);
      for (i=0;i<jobtotal;++i)
      { fileout.integer(i+1,4,16);
        fileout.space(3);
        fileout.string((finname+i*100),0,1);
        if (errorcode[i]==0)
        { fileout.string("       result saved in file ");
          fileout.string((foutname+i*100),0,1);
        }
        else
        { Errorinfo errorinfo(errorcode[i]);
          fileout.string(errorinfo.errmsg,0,1);
        }
      }
      if (error==0) error=fileout.geterror();
    }
    else error+=(datatype<<16);
  }
  return(error);
} //end of Mscdjob::listjob

int Mscdjob::execute(char *usermessage)
{ int jobnum;
  float xa;

  Textout conout;
  for (jobnum=0;(error==0)&&(jobnum<jobtotal);++jobnum)
  { if (mype==0)
    { if (dispmode>4) waitenter();
      if ((dispmode==0)||(dispmode>1))
      { conout.charfill('_',37,16+1);
        conout.string("Job "); conout.integer(jobnum+1,3);
        conout.string(" of "); conout.integer(jobtotal,3,1);
      }
      if (flogout)
      { flogout->charfill('_',37,16+1);
        flogout->string("Job "); flogout->integer(jobnum+1,3);
        flogout->string(" of "); flogout->integer(jobtotal,3,1);
      }
    }
    if (error==0)
    { Mscdrun *mscdrun=new Mscdrun(mype,numpe);
      if (!mscdrun) error=102;
      if (mype==0)
      { mscdrun->init();
        error=mscdrun->loadvars((finname+jobnum*100),
          (foutname+jobnum*100),jobnum,jobtotal,dispmode,displog,
          flogout);
        if (error==0)
          error=mscdrun->readparameter();
        if (error==0)
          error=mscdrun->paradisp();
        if (error==0)
          error=mscdrun->paradisplog();
        if (error==0)
          error=mscdrun->symtrivert();
        if (error==0)
          error=mscdrun->symdblvert();
        if (error==0)
          error=mscdrun->precutable();
        if (error==0)
        { xa=mscdrun->getmemory();
          if (basemem<xa) basemem=xa;
        }
      }
      if ((error==0)&&(numpe>1)&&(mype==0))
        error=mscdrun->sendjobs();
      else if ((numpe>1)&&(mype==0))
        mscdrun->sendjobs();
      else if (numpe>1)
        error=mscdrun->receivejobs();
      if (error==0)
        error=mscdrun->assistant(1);
      if ((error==0)&&(mype==0))
        error=mscdrun->fitarpefs(mscdrun);
      if ((error==0)&&(mype==0))
        error=mscdrun->assistant(0);
      if ((error==0)&&(mype==0))
        error=mscdrun->savecurve((foutname+jobnum*100),
          usermessage);
      errorcode[jobnum]=error; error=0;
      if (mscdrun) delete mscdrun;
    }
  }

  return(error);
} //end of Mscdjob::execute

