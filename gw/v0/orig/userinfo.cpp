#include <iostream>

#include "userutil.h"
#include "userinfo.h"

Userinfo::Userinfo(int imype,int inumpe)
{ mype=imype; numpe=inumpe; error=107; init();
} //end of Userinfo::Userinfo

Userinfo::~Userinfo()
{ if (author) delete [] author; if (usermsg) delete [] usermsg;
} //end of Userinfo::~Userinfo

void Userinfo::init()
{ if ((error==107)||(error==103))
  { author=usermsg=NULL;
  }
  if (error==103)
  { author=new char [80]; usermsg=new char [1024];
    if ((!author)||(!usermsg)) error=102;
  }
  else if (error==106)
  { if (author)
    { author=new char [80]; if (!author) error=102;
    }
    if (usermsg)
    { usermsg=new char [1024]; if (!usermsg) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Userinfo::init

int Userinfo::loadparameter(float iversion,char *filename,
  char *iauthor,char *title,int machine)
{ int i,j,ka,kb,kc;
  char dashline[80];

  if ((error==101)||(error==0))
  { error=0;
    if (iversion<1.00) version=1.0f;
    else if (iversion>10.0) version=10.0f;
    else version=iversion;
    version+=1.0e-3f;
    ka=0; usermsg[0]='\0';
    for (i=0;i<52;++i)
    { if (iauthor[i]!='\0') author[i]=iauthor[i];
      else break;
    }
    author[i]='\0';
    for (i=0;i<65;++i) dashline[i]='-';
    dashline[i]='\0';
    ka+=stringlength(dashline)+stringlength(title)+
      stringlength(filename)+12;
    if (ka<1000)
      ka=stringappend(usermsg,dashline,0,"\n",title,"\n",filename);

    kc=stringlength(author)+stringlength(filename);
    if (error==0) kb=sizeof(int); else kb=sizeof(int);
    if (machine>100)
    { ka+=10; kc+=10;
      if (ka<1000) ka=stringappend(usermsg," Parallel ");
    }
    else if (kb<4)
    { ka+=8; kc+=8;
      if (ka<1000) ka=stringappend(usermsg," 16 Bit ");
    }
    if ((ka<1000)&&(kc<52)) ka=stringappend(usermsg," Version ");
    else if (ka<1000) ka=stringappend(usermsg," ");
    i=(int)version; j=(int)(version*10.0)-i*10;
    kb=(int)(version*100.0)-i*100-j*10;
    if (ka<1000)
    { usermsg[ka-1]=(char)(i+'0');
      usermsg[ka]='.';
      usermsg[ka+1]=(char)(j+'0');
      usermsg[ka+2]=(char)(kb+'0');
      usermsg[ka+3]=' ';
      usermsg[ka+4]='\0';
    }
    ka+=5+stringlength(author)+200;
    if (ka<1000) ka=stringappend(usermsg,author,0,"\n");
    ka=stringappend(usermsg,"Lawrence Berkeley National Laboratory",0,
      " (LBNL), Berkeley, CA 94720\nCopyright (c) Van Hove Group",
      " 1997-1998. All rights reserved\n",dashline,"\n");
  }
  else
  { version=0.0f;
  }
  return(error);
} //end of Userinfo::Userinfo

float Userinfo::getmemory()
{ float xa;
  xa=(float)(sizeof(Userinfo)+(1024+64)*sizeof(char));
  return(xa);
} //end of Userinfo::getmemory

Errorinfo::Errorinfo(int errorcode)
{ int ka,filecode,subcode,length;

  error=errorcode; errmsg=NULL;
  length=100; errmsg=new char [length];
  if ((errmsg)&&(error!=0))
  { ka=0; errmsg[0]='\0';
    subcode=(errorcode&0xffff);
    filecode=(errorcode>>16);
    if ((subcode>=1000)||(subcode==0)) subcode=901;
    ka=stringappend(errmsg,"Error ");
    errmsg[ka-1]=(char)(subcode/100+'0');
    errmsg[ka]=(char)(subcode/10%10+'0');
    errmsg[ka+1]=(char)(subcode%10+'0');
    errmsg[ka+2]='\0';
    ka+=3;
    if ((subcode/10 != 20) && (subcode/10 != 70) &&
      (subcode/10 != 71)) filecode=0;
    else if ((filecode > 100) && (filecode < 300))
      ka=stringappend(errmsg,", photoemission");
    else if ((filecode >= 300) && (filecode < 500))
      ka=stringappend(errmsg,", photoemission");
    else if ((filecode >= 710) && (filecode < 720))
      ka=stringappend(errmsg,", phase shift");
    else if ((filecode >= 720) && (filecode < 730))
      ka=stringappend(errmsg,", radial matrix");
    else if ((filecode >= 730) && (filecode < 740))
      ka=stringappend(errmsg,", report");
    else if ((filecode >= 740) && (filecode < 750))
      ka=stringappend(errmsg,", mscd input");
    else if ((filecode >= 750) && (filecode < 760))
      ka=stringappend(errmsg,", batch input");
    else if ((filecode >= 810) && (filecode < 820))
      ka=stringappend(errmsg,", potential");
    else if ((filecode >= 820) && (filecode < 830))
      ka=stringappend(errmsg,", psrm input");
    else if ((filecode >= 830) && (filecode < 840))
      ka=stringappend(errmsg,", eigen wave");
    else if ((filecode >= 840) && (filecode < 850))
      ka=stringappend(errmsg,", spectrum");
    else if ((filecode >= 850) && (filecode < 860))
      ka=stringappend(errmsg,", peak intensity");
    else if ((filecode >= 910) && (filecode < 920))
      ka=stringappend(errmsg,", hologram");
    else if ((filecode >= 920) && (filecode < 930))
      ka=stringappend(errmsg,", holo input");
    else if ((filecode >= 930) && (filecode < 940))
      ka=stringappend(errmsg,", scattering factor");
    else if ((filecode >= 940) && (filecode < 950))
      ka=stringappend(errmsg,", mean free path");
    else if ((filecode >= 950) && (filecode < 960))
      ka=stringappend(errmsg,", thermal vibration");
    else if (filecode==961)
      ka=stringappend(errmsg,", expectation");
    else if (filecode==962)
      ka=stringappend(errmsg,", gaussain output");
    else if (filecode==963)
      ka=stringappend(errmsg,", relaxation");
    else if (filecode==964)
      ka=stringappend(errmsg,", rys roots");
    else
      ka=stringappend(errmsg,",");

    if (subcode == 101)
      ka=stringappend(errmsg," data not ready");
    else if (subcode == 102)
      ka=stringappend(errmsg," out of memory");
    else if (subcode == 103)
      ka=stringappend(errmsg," memory not ready");
    else if (subcode == 104)
      ka=stringappend(errmsg," out of allocated memory");
    else if (subcode == 105)
      ka=stringappend(errmsg," 32 bit compiling required");
    else if (subcode == 106)
      ka=stringappend(errmsg," memory not paimported");
    else if (subcode == 107)
      ka=stringappend(errmsg," memory not flushed");
    else if (subcode == 108)
      ka=stringappend(errmsg," floating number size mismatch");
    else if (subcode == 201)
      ka=stringappend(errmsg," file not found");
    else if (subcode == 202)
      ka=stringappend(errmsg," file open error");
    else if (subcode == 203)
      ka=stringappend(errmsg," file creation error");
    else if (subcode == 204)
      ka=stringappend(errmsg," data reading error");
    else if (subcode == 205)
      ka=stringappend(errmsg," data writing error");
    else if (subcode == 206)
      ka=stringappend(errmsg," end-of-file reading error");
    else if (subcode == 207)
      ka=stringappend(errmsg," data format error");
    else if (subcode == 208)
      ka=stringappend(errmsg," no data available");
    else if (subcode == 209)
      ka=stringappend(errmsg," data type error");
    else if (subcode == 211)
      ka=stringappend(errmsg," singular matrix error");
    else if (subcode == 221)
      ka=stringappend(errmsg," incorrect application type");
    else if (subcode == 222)
      ka=stringappend(errmsg," form contain no data");
    else if (subcode == 223)
      ka=stringappend(errmsg," form submitted with unknown method");
    else if (subcode == 224)
      ka=stringappend(errmsg," required field not submitted");
    else if (subcode == 225)
      ka=stringappend(errmsg," access denied from your web server");
    else if (subcode == 226)
      ka=stringappend(errmsg," invalid email address");
    else if (subcode == 227)
      ka=stringappend(errmsg," invalid url address");
    else if (subcode == 231)
      ka=stringappend(errmsg," invalid input format");
    else if (subcode == 232)
      ka=stringappend(errmsg," invalid output format");
    else if (subcode == 233)
      ka=stringappend(errmsg," no input error");
    else if (subcode == 234)
      ka=stringappend(errmsg," no output error");
    else if (subcode == 235)
      ka=stringappend(errmsg," memory not enough for output");
    else if (subcode == 241)
      ka=stringappend(errmsg," password too simple");
    else if (subcode == 601)
      ka=stringappend(errmsg," too few atoms (no or one)");
    else if (subcode == 602)
      ka=stringappend(errmsg," too many atoms");
    else if (subcode == 603)
      ka=stringappend(errmsg," too many kinds of atoms");
    else if (subcode == 604)
      ka=stringappend(errmsg," emitter not found");
    else if (subcode == 605)
      ka=stringappend(errmsg," too small distance between atoms");
    else if (subcode == 611)
      ka=stringappend(errmsg," no available experimental data");
    else if (subcode == 612)
      ka=stringappend(errmsg," experimental chi data too small");
    else if (subcode == 613)
      ka=stringappend(errmsg," data files not match");
    else if (subcode == 614)
      ka=stringappend(errmsg," too few input data points");
    else if (subcode == 615)
      ka=stringappend(errmsg," too many input data points");
    else if (subcode == 616)
      ka=stringappend(errmsg," too few output data points");
    else if (subcode == 617)
      ka=stringappend(errmsg," too many output data points");
    else if (subcode == 618)
      ka=stringappend(errmsg," no available hologram data");
    else if (subcode == 621)
      ka=stringappend(errmsg," symmetry search error");
    else if (subcode == 622)
      ka=stringappend(errmsg," symmetry order error");
    else if (subcode == 623)
      ka=stringappend(errmsg," subshell name error");
    else if (subcode == 624)
      ka=stringappend(errmsg," binding energy adjustment error");
    else if (subcode == 625)
      ka=stringappend(errmsg," binding energy may be too small");
    else if (subcode == 626)
      ka=stringappend(errmsg," binding energy may be too large");
    else if (subcode == 627)
      ka=stringappend(errmsg," binding energy be too small");
    else if (subcode == 628)
      ka=stringappend(errmsg," binding energy be too large");
    else if (subcode == 631)
      ka=stringappend(errmsg," output file name not found");
    else if (subcode == 632)
      ka=stringappend(errmsg," duplicate output file name");
    else if (subcode == 633)
      ka=stringappend(errmsg," list file open error");
    else if (subcode == 641)
      ka=stringappend(errmsg," too few fitting parameters");
    else if (subcode == 642)
      ka=stringappend(errmsg," too many fitting parameters");
    else if (subcode == 651)
      ka=stringappend(errmsg," atom not found");
    else if (subcode == 701)
      ka=stringappend(errmsg," too few data points");
    else if (subcode == 702)
      ka=stringappend(errmsg," too many data points");
    else if (subcode == 703)
      ka=stringappend(errmsg," too few curves");
    else if (subcode == 704)
      ka=stringappend(errmsg," too many curves");
    else if (subcode == 705)
      ka=stringappend(errmsg," no available data");
    else if (subcode == 706)
      ka=stringappend(errmsg," data too small");
    else if (subcode == 707)
      ka=stringappend(errmsg," data order error");
    else if (subcode == 711)
      ka=stringappend(errmsg," energy order error");
    else if (subcode == 712)
      ka=stringappend(errmsg," energy step too small");
    else if (subcode == 713)
      ka=stringappend(errmsg," too few energy points");
    else if (subcode == 714)
      ka=stringappend(errmsg," too many energy points");
    else if (subcode == 715)
      ka=stringappend(errmsg," too few angle points");
    else if (subcode == 716)
      ka=stringappend(errmsg," too many angle points");
    else if (subcode == 717)
      ka=stringappend(errmsg," zero or negative energy");
    else if (subcode == 751)
      ka=stringappend(errmsg," message sending error");
    else if (subcode == 752)
      ka=stringappend(errmsg," message receiving error");
    else if (subcode == 753)
      ka=stringappend(errmsg," received empty message");
    else if (subcode == 754)
      ka=stringappend(errmsg," received error message");
    else if (subcode == 761)
      ka=stringappend(errmsg," too many processors");
    else if (subcode == 901)
      ka=stringappend(errmsg," programming error occured");
    else if (subcode == 999)
      ka=stringappend(errmsg," unknown error occured");
    else
      ka=stringappend(errmsg," unknown error occured");

    ka=stringappend(errmsg,", program terminated");
  }
  else if (errmsg)
  { ka=stringcopy(errmsg,"Program terminated normally");
  }
  else if (error==0) error=102;
} //end of Errorinfo::Errorinfo

Errorinfo::~Errorinfo()
{ delete [] errmsg;
} //end of Errorinfo::~Errorinfo

float Errorinfo::getmemory()
{ float xa;
  xa=(float)(sizeof(Errorinfo)+100*sizeof(char));
  return(xa);
} //end of Userinfo::getmemory


