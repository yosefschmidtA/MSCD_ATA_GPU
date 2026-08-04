#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "userutil.h"
using namespace std;


Textout::Textout()
{ error=0; fdout=0;
} //end of Textout::Textout

Textout::Textout(const char *filename, std::ios_base::openmode mode)
{ fout.open(filename,mode);
  if (!fout)
  { error=203; fdout=0;
  }
  else
  { error=0; fdout=10;
  }
} //end of Textout::Textout

int Textout::open(const char *filename, std::ios_base::openmode mode)
{ if (fdout<1) fout.open(filename,mode);
  if (!fout)
  { error=203; fdout=0;
  }
  else
  { error=0; fdout=10;
  }
  return(error);
} //end of Textout::Textout

int Textout::geterror()
{ return(error);
} //end of Textout::geterror

int Textout::newline()
{ if ((error==0)&&(fdout<1)) cout << endl;
  else if (error==0)
  { fout << endl;
    if (!fout) error=205;
  }
  return(error);
} //end of Textout::newline

int Textout::newline(int linenum)
{ int i;
  if ((error==0)&&(fdout<1))
  { for (i=0;i<linenum;++i) cout << endl;
  }
  else if (error==0)
  { for (i=0;i<linenum;++i) fout << endl;
    if (!fout) error=205;
  }
  return(error);
} //end of Textout::newline

int Textout::space(int width,int linenum)
{ int i,begline,endline;
  begline=(linenum>>4)&0xf; endline=linenum&0xf;
  if ((error==0)&&(fdout<1))
  { for (i=0;i<begline;++i) cout << endl;
    for (i=0;i<width;++i) cout << ' ';
    for (i=0;i<endline;++i) cout << endl;
  }
  else if (error==0)
  { for (i=0;i<begline;++i) fout << endl;
    for (i=0;i<width;++i) fout << ' ';
    for (i=0;i<endline;++i) fout << endl;
    if (!fout) error=205;
  }
  return(error);
} //end of Textout::space

int Textout::charfill(char ch,int width,int linenum)
{ int i,begline,endline;
  begline=(linenum>>4)&0xf; endline=linenum&0xf;
  if ((error==0)&&(fdout<1))
  { for (i=0;i<begline;++i) cout << endl;
    for (i=0;i<width;++i) cout << ch;
    for (i=0;i<endline;++i) cout << endl;
  }
  else if (error==0)
  { for (i=0;i<begline;++i) fout << endl;
    for (i=0;i<width;++i) fout << ch;
    for (i=0;i<endline;++i) fout << endl;
    if (!fout) error=205;
  }
  return(error);
} //end of Textout::charfill

int Textout::string(char *string,int width,int linenum)
{ int i,begline,endline;
  begline=(linenum>>4)&0xf; endline=linenum&0xf;
  if ((error==0)&&(fdout<1))
  { for (i=0;i<begline;++i) cout << endl;
    if (width>0) cout << setw(width);
    cout << string;
    for (i=0;i<endline;++i) cout << endl;
  }
  else if (error==0)
  { for (i=0;i<begline;++i) fout << endl;
    if (width>0) fout << setw(width);
    fout << string;
    for (i=0;i<endline;++i) fout << endl;
    if (!fout) error=205;
  }
  return(error);
} //end of Textout::string

int Textout::integer(int number,int width,int linenum)
{ int i,begline,endline;
  begline=(linenum>>4)&0xf; endline=linenum&0xf;
  if ((error==0)&&(fdout<1))
  { for (i=0;i<begline;++i) cout << endl;
    if (width>0) cout << setw(width);
    cout << number;
    for (i=0;i<endline;++i) cout << endl;
  }
  else if (error==0)
  { for (i=0;i<begline;++i) fout << endl;
    if (width>0) fout << setw(width);
    fout << number;
    for (i=0;i<endline;++i) fout << endl;
    if (!fout) error=205;
  }
  return(error);
} //end of Textout::integer

int Textout::floating(float realnumber,float widthprecision,
  int linenum)
{ int i,k,begline,endline,width,precision,format;

  begline=(linenum>>4)&0xf; endline=linenum&0xf;
  width=(int)widthprecision;
  precision=(int)((widthprecision-width)*10+0.1);
  if ((precision>width-2)&&(width>2)) precision=width-2;
  else if (precision>width-2) precision=0;
  if (precision>20) precision=20;
  format=(linenum>>8)&3;
  if (realnumber>0.1) k=(int)log10((double)realnumber);
  else if (realnumber<-0.1) k=(int)log10(-(double)realnumber);
  else k=0;
  if ((format==3)&&(k>precision)) format=2;
  else if ((format==3)&&(k>0)) precision-=k;
  if ((error==0)&&(fdout<1))
  { for (i=0;i<begline;++i) cout << endl;
    if (width>0) cout << setw(width);
    if (precision>0) cout << setprecision(precision);
    if ((precision>0)&&((format==1)||(format==3)))
      cout << setiosflags(ios::fixed) << resetiosflags(ios::scientific);
    else if ((precision>0)&&(format==2))
      cout << resetiosflags(ios::fixed) << setiosflags(ios::scientific);
    else
      cout << resetiosflags(ios::fixed|ios::scientific);
    if ((width>0)&&(precision<1)) cout << (int)realnumber;
    else cout << realnumber;
    for (i=0;i<endline;++i) cout << endl;
  }
  else if (error==0)
  { for (i=0;i<begline;++i) fout << endl;
    if (width>0) fout << setw(width);
    if (precision>0) fout << setprecision(precision);
    if ((precision>0)&&(format==1))
      fout << setiosflags(ios::fixed) << resetiosflags(ios::scientific);
    else if ((precision>0)&&(format==2))
      fout << resetiosflags(ios::fixed) << setiosflags(ios::scientific);
    else
      fout << resetiosflags(ios::fixed|ios::scientific);
    if ((width>0)&&(precision<1)) fout << (int)realnumber;
    else fout << realnumber;
    for (i=0;i<endline;++i) fout << endl;
    if (!fout) error=205;
  }
  return(error==0);
} //end of Textout::floating

int Textout::flush()
{ 
  if ((error==0)&&(fdout<1)) cout.flush();
  else if (error==0)
  { fout.flush();
    if (!fout) error=205;
  }
  return(error);
} //end of Textout::flush


int iabs(int ix)
{ return ((ix>=0)?ix:-ix);
} //end of iabs

char confine(char x, char min, char max)
{
  if (x < min) x=min;
  else if (x > max) x=max;
  return(x);
} //end of confine

int confine(int x, int min, int max)
{
  if (x < min) x=min;
  else if (x > max) x=max;
  return(x);
} //end of confine

float confine(float x, float min, float max)
{
  if (x < min) x=min;
  else if (x > max) x=max;
  return(x);
} //end of confine

// ya=decimal*(xa/decimal+0.5)
float round(float xa,float decimal)
{ int k;
  float xb;
  if (decimal<0.0) decimal=-decimal;
  xb=confine(decimal,1.0e-10f,1.0e10f);
  if (xa>=0.0) k=(int)(xa/xb+0.5);
  else k=(int)(xa/xb-0.5);
  xa=k*xb;
  return(xa);
} //end of round

int spacehide(char *membuf)
{ int i,j,k,m,ch;
  i=k=m=0;
  while ((membuf)&&(i<256)&&(membuf[i]!='\0'))
  { ch=membuf[i];
    if ((k==0)&&((ch=='\'')||(ch == '\"')))
    { k=ch; m=i;
    }
    else if ((k==0)&&((ch==',')||(ch=='\t')||(ch=='\n')||(ch=='\r')))
      membuf[i]=' ';
    else if ((k!=0)&&(ch==k))
    { for (j=m;j<=i;++j)
      { ch=membuf[j];
        if (ch==k) membuf[j]=' ';
        else if (ch==' ') membuf[j]='\"';
      }
      k=0;
    }
    ++i;
  }
  if ((membuf)&&(k!=0)&&(membuf[m]=='\"')) membuf[m]=' ';
  
  return(i+1);
} //end of spacehide

int spaceback(char *membuf)
{ int i;
  i=0;
  while ((membuf)&&(i<256)&&(membuf[i]!='\0'))
  { if (membuf[i]=='\"') membuf[i]=' ';
    ++i;
  }
  return(i+1);
} //end of spaceback

int spaceclear(char *membuf)
{ int i,j,k;
  k=0; while ((membuf)&&(k<256)&&(membuf[k]!='\0')) ++k;
  ++k;
  if (k>1)
  { for (i=0;i<256;++i)
      if ((membuf[i]!=' ')&&(membuf[i]!='\t')&&(membuf[i]!='\n')&&
        (membuf[i]!='\r')&&(membuf[i]!='\'')&&(membuf[i]!='\"'))
        break;
    for (j=k-2;j>=0;--j)
      if ((membuf[j]!=' ')&&(membuf[j]!='\t')&&(membuf[j]!='\n')&&
        (membuf[i]!='\r')&&(membuf[j]!='\'')&&(membuf[j]!='\"'))
        break;
    if (i>j) k=1;
    else k=j-i+2;
    if (i>0) for (j=0;j<k-1;++j) membuf[j]=membuf[j+i];
    membuf[k-1]='\0';
  }
  else
  { k=1; membuf[k-1]='\0';
  }

  return(k);
} //end of spaceclear

void skipdelim(istream& is)
{ char ch;
  ch=' ';
  while ((is)&&((ch==' ')||(ch=='\t')||(ch=='\n')||(ch=='\r')||
    (ch=='(')||(ch==')')||(ch==','))) is>>ch;
  if (is) is.putback(ch);
  return;
} //end of skipdelim

int skipendline(istream& is)
{ char ch;
  ch=0;
  while ((is)&&(ch!='\n')&&(ch!='\r')) is.get(ch);
  return(ch);
} //end of skipendline

int getendline(istream& is,char *membuf,int memsize)
{ int i;
  char ch;

  i=0; ch=' ';
  while ((is)&&(membuf)&&(ch!='\n')&&(ch!='\r')&&
    ((memsize==0)||(i<memsize-1)))
  { is.get(ch); membuf[i]=ch;
    if ((ch!='\n')&&(ch!='\r')) ++i;
  }
  membuf[i]='\0';
  if (is) ch=0; else ch=-1;
  return(ch);
} //end of skipendline

int skipbytes(istream& is,int nbyte)
{ int i,k;
  char ch;

  if ((nbyte<1)||(nbyte>256)) k=256; else k=nbyte;
  ch=0;
  if (is)
  { for (i=0;i<k;++i)
    { if ((ch!='\n')&&(ch!='\r')) is.get(ch);
      else i=k;
    }
  }

  return(ch);
} //end of skipbytes

int stringlength(const char *source)
{ int i;
  i=0; while ((source)&&(source[i]!='\0')) ++i;
  return(i+1);
} //end of stringlength

int stringcopy(char *dest,const char *source,int nbyte)
{ int k;
  k=0;
  if ((dest)&&(source))
  { while (((nbyte==0)||(k<nbyte-1))&&(source[k]!='\0'))
    { dest[k]=source[k]; ++k;
    }
    dest[k]='\0'; ++k;
  }
  return(k);
} //end of stringcopy

int stringappend(char *dest,const char *source,int nbyte,
  const char *sourcea,const char *sourceb,const char *sourcec,
  const char *sourced)
{ int ka,kb,kc;

  if ((dest)&&(source))
  { ka=stringlength(dest); kb=stringlength(source);
    if ((nbyte!=0)&&(kb>nbyte)) kb=nbyte;
    stringcopy(dest+ka-1,source,kb);
    ka+=kb-1;
    if (sourcea)
    { kc=stringlength(sourcea);
      if ((nbyte!=0)&&(kc>nbyte-kb)) kc=nbyte-kb;
      stringcopy(dest+ka-1,sourcea,kc);
      ka+=kc-1; kb+=kc-1;
    }
    if (sourceb)
    { kc=stringlength(sourceb);
      if ((nbyte!=0)&&(kc>nbyte-kb)) kc=nbyte-kb;
      stringcopy(dest+ka-1,sourceb,kc);
      ka+=kc-1; kb+=kc-1;
    }
    if (sourcec)
    { kc=stringlength(sourcec);
      if ((nbyte!=0)&&(kc>nbyte-kb)) kc=nbyte-kb;
      stringcopy(dest+ka-1,sourcec,kc);
      ka+=kc-1; kb+=kc-1;
    }
    if (sourced)
    { kc=stringlength(sourced);
      if ((nbyte!=0)&&(kc>nbyte-kb)) kc=nbyte-kb;
      stringcopy(dest+ka-1,sourced,kc);
      ka+=kc-1; kb+=kc-1;
    }
  }
  else ka=0;
  return(ka);
} //end of stringappend

int stringcomp(const char *source1,const char *source2,int nbyte)
{ int i,k,cha,chb;
  i=k=0;
  if ((!source1)&&(!source2)) k=0;
  else if (!source1) k=-1;
  else if (!source2) k=1;
  else
  { while (k==0)
    { cha=source1[i]; chb=source2[i]; ++i;
      if ((nbyte!=0)&&(i>nbyte-1)) break;
      else if ((cha=='\0')&&(chb=='\0')) break;
      else if ((cha=='\0')&&(chb!='\0')) k=-1;
      else if ((cha!='\0')&&(chb=='\0')) k=1;
      else
      { if ((cha>='A')&&(cha<='Z')) cha+='a'-'A';
        if ((chb>='A')&&(chb<='Z')) chb+='a'-'A';
        if (cha>chb) k=1;
        else if (cha<chb) k=-1;
      }
    }
  }
  return(k);
} //end of strigcomp

int stringclear(char *membuf)
{ char cha,chb;
  int i,k;
  i=k=0;
  while ((membuf)&&(i<256)&&(membuf[i]!='\0'))
  { if (k>0) cha=membuf[k-1]; else cha=' ';
    chb=membuf[i++];
    if ((chb==' ')||(chb=='\t')||(chb=='\n')||(chb=='\r')||
      (chb==',')||(chb=='=')||(chb=='\\')||(chb=='/')||(chb==':')||
      (chb==';')||(chb=='*')||(chb=='-')||(chb=='_'))
    { if (cha!=' ') membuf[k++]=' ';
    }
    else membuf[k++]=chb;
  }
  if ((k>0)&&(membuf[k-1]==' ')) --k;
  membuf[k++]='\0';
  return(k);
} //end of stringclear

int stringtoint(char *source)
{ int i,j,k,base,ma;

  if (source) k=stringlength(source);
  else k=0;
  i=0; base=10; ma=0;
  while ((i<k)&&((source[i]==' ')||(source[i]=='\t')||
    (source[i]=='\r')||(source[i]=='\n')||(source[i]==','))) ++i;
  if ((i<k-2)&&(source[i]=='0')&&((source[i+1]=='x')||
    (source[i+1]=='X')))
  { base=16; i+=2;
  }
  while (i<k)
  { if ((source[i]>='0')&&(source[i]<='9'))
      j=(int)(source[i]-'0');
    else if ((source[i]>='a')&&(source[i]<='f'))
      j=10+(int)(source[i]-'a');
    else if ((source[i]>='A')&&(source[i]<='F'))
      j=10+(int)(source[i]-'A');
    else if (source[i]==',') j=-1;
    else j=base+100;
    if ((j>=0)&&(j<base)) ma=ma*base+j;
    ++i;
  }
  return(ma);
} //end of stringtoint

float stringtofloat(char *source)
{ float xa;

  if (source) xa=(float)strtod(source,NULL);
  else xa=0.0f;
  return(xa);
} //end of stringtofloat

int stringmatch(char *source1,char *source2)
{ int i,ja,jb,jc,jd,ka,kb,kc,match;

  match=0;
  if (source1) ka=stringlength(source1);
  else ka=0;
  if (source2) kb=stringlength(source2);
  else kb=0;
  if ((ka>0)&&(kb>0))
  { ja=nextnonwhite(source1,0,ka);
    jb=nextnonwhite(source1,ka-1,0,1)+1;
    jc=nextnonwhite(source2,0,kb);
    jd=nextnonwhite(source2,kb-1,0,1)+1;
    ka=jb-ja; kb=jd-jc;
    if (ka>kb) kc=ka-kb;
    else kc=kb-ka;

    if (ka>kb)
    { for (i=ja;i<ja+kc;++i)
      { if (stringcomp(source1+i,source2+jc,kb)==0)
        { match=1; i=ja+kc;
        }
      }
    }
    else
    { for (i=jc;i<jc+kc;++i)
      { if (stringcomp(source1+ja,source2+i,ka)==0)
        { match=1; i=jc+kc;
        }
      }
    }
  }

  return(match);
} //end of stringmatch

int nextchar(char *source,char chx,int begadd,int endadd,
  int direction)
{ int i,ka,cha,chb;

  cha=chb=chx;
  if ((cha>='A')&&(cha<='Z')) cha=cha-'A'+'a';
  if ((chb>='a')&&(chb<='z')) chb=chb-'a'+'A';
  if (source) ka=begadd+stringlength(source+begadd);
  else ka=0;
  if (begadd<0) begadd=0;
  else if (begadd>ka) begadd=ka;
  if (endadd<0) endadd=0;
  else if (endadd>ka) endadd=ka;
  if ((direction==0)&&((endadd<=0)||(endadd>ka))) endadd=ka;
  else if ((direction!=0)&&(endadd<0)) endadd=0;
  if ((direction==0)&&(begadd<endadd))
  { i=begadd;
    while ((i<endadd)&&(source[i]!=cha)&&(source[i]!=chb)) ++i;
  }
  else if ((direction!=0)&&(begadd>endadd))
  { i=begadd;
    while ((i>=endadd)&&(source[i]!=cha)&&(source[i]!=chb)) --i;
  }
  else i=endadd;
  return(i);
} //end of nextchar

int nextnonchar(char *source,char chx,int begadd,int endadd,
  int direction)
{ int i,ka,cha,chb;

  cha=chb=chx;
  if ((cha>='A')&&(cha<='Z')) cha=cha-'A'+'a';
  if ((chb>='a')&&(chb<='z')) chb=chb-'a'+'A';
  if (source) ka=begadd+stringlength(source+begadd);
  else ka=0;
  if (begadd<0) begadd=0;
  else if (begadd>ka) begadd=ka;
  if (endadd<0) endadd=0;
  else if (endadd>ka) endadd=ka;
  if ((direction==0)&&((endadd<=0)||(endadd>ka))) endadd=ka;
  else if ((direction!=0)&&(endadd<0)) endadd=0;
  if ((direction==0)&&(begadd<endadd))
  { i=begadd;
    while ((i<endadd)&&((source[i]==cha)||(source[i]==chb))) ++i;
  }
  else if ((direction!=0)&&(begadd>endadd))
  { i=begadd;
    while ((i>=endadd)&&((source[i]==cha)||(source[i]==chb))) --i;
  }
  else i=endadd;
  return(i);
} //end of nextnonchar

int nextwhite(char *source,int begadd,int endadd,int direction)
{ int i,ka;
  if (source) ka=begadd+stringlength(source+begadd);
  else ka=0;
  if (begadd<0) begadd=0;
  else if (begadd>ka) begadd=ka;
  if (endadd<0) endadd=0;
  else if (endadd>ka) endadd=ka;
  if ((direction==0)&&((endadd<=0)||(endadd>ka))) endadd=ka;
  else if ((direction!=0)&&(endadd<0)) endadd=0;
  if ((direction==0)&&(begadd<endadd))
  { i=begadd;
    while ((i<endadd)&&(source[i]!=' ')&&(source[i]!='\t')&&
      (source[i]!='\r')&&(source[i]!='\n')) ++i;
  }
  else if ((direction!=0)&&(begadd>endadd))
  { i=begadd;
    while ((i>=endadd)&&(source[i]!=' ')&&(source[i]!='\t')&&
      (source[i]!='\r')&&(source[i]!='\n')) --i;
  }
  else i=endadd;
  if (i<0) i=0;
  else if (i>ka-1) i=ka-1;

  return(i);
} //end of nextwhit

int nextnonwhite(char *source,int begadd,int endadd,int direction)
{ int i,ka;
  if (source) ka=begadd+stringlength(source+begadd);
  else ka=0;
  if (begadd<0) begadd=0;
  else if (begadd>ka) begadd=ka;
  if (endadd<0) endadd=0;
  else if (endadd>ka) endadd=ka;
  if ((direction==0)&&((endadd<=0)||(endadd>ka))) endadd=ka;
  else if ((direction!=0)&&(endadd<0)) endadd=0;
  if ((direction==0)&&(begadd<endadd))
  { i=begadd;
    while ((i<endadd)&&((source[i]==' ')||(source[i]=='\t')||
      (source[i]=='\r')||(source[i]=='\n'))) ++i;
  }
  else if ((direction!=0)&&(begadd>endadd))
  { i=begadd;
    while ((i>=endadd)&&((source[i]==' ')||(source[i]=='\t')||
      (source[i]=='\r')||(source[i]=='\n'))) --i;
  }
  else i=endadd;
  if (i<0) i=0;
  else if (i>ka-1) i=ka-1;
  
  return(i);
} //end of nextnonwhite

int memorycopy(char *dest,const char *source,int nbyte)
{ int i,k;
  if ((dest)&&(source))
  { k=nbyte;
    for (i=0;i<k;++i) dest[i]=source[i];
  }
  else k=0;
  return(k);
} //end of memeorycopy

int memorysend(char *sendbuf,char *source,int address,int nbyte,
  int destsize)
{ if ((sendbuf)&&(source)&&(address+nbyte<=destsize))
  { memorycopy(sendbuf+address,source,nbyte);
    address+=nbyte;
  }
  else if ((sendbuf)&&(source)) address=-1;
  return(address);
} //end of memorysend

int memoryrec(char *dest,char *recbuf,int address,int nbyte,
  int destsize)
{ if ((dest)&&(recbuf)&&(address+nbyte<=destsize))
  { memorycopy(dest,recbuf+address,nbyte);
    address+=nbyte;
  }
  else if ((dest)&&(recbuf)) address=-1;
  return(address);
} //end of memoryrec

void waitenter()
{ int k;

  cout << endl << "Pause --- strike Enter to continue ";
  k=0;
  while ((k!='\n')&&(k!='\r')) k=cin.get();
  cout << endl;
  return;
} //end of waitenter

