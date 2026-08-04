#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "userutil.h"
#include "polation.h"
#include "xpspec.h"

Spectrum::Spectrum()
{ error=107; init();
} //end of Spectrum::Spectrum

Spectrum::~Spectrum()
{ if (ebind) delete [] ebind; if (xpscount) delete [] xpscount;
  if (curvepar) delete [] curvepar; if (comment) delete [] comment;
  if (emission) delete [] emission; if (bestfit) delete [] bestfit;
} //end of Spectrum::~Spectrum

void Spectrum::init()
{ if ((error==107)||(error==103))
  { datatype=0; ncurve=npoint=npeak=numint=commentsize=0;
    comment=NULL;
    ebind=xpscount=curvepar=emission=bestfit=NULL;
    if (error==0) mpoint=sizeof(int);
    else mpoint=sizeof(int);
    if (mpoint<4) mpoint=5000;
    else mpoint=25000;
  }
  if (error==103)
  { ebind=new float [mpoint]; xpscount=new float [mpoint];
    curvepar=new float [256*8];
    if ((!ebind)||(!xpscount)||(!curvepar)) error=102;
  }
  else if (error==106)
  { if (ebind)
    { ebind=new float [mpoint];
      if (!ebind) error=102;
    }
    if (xpscount)
    { xpscount=new float [mpoint];
      if (!xpscount) error=102;
    }
    if (curvepar)
    { curvepar=new float [256*8];
      if (!curvepar) error=102;
    }
    if (comment)
    { comment=new char [commentsize];
      if (!comment) error=102;
    }
    if (emission)
    { emission=new float [ncurve*numint];
      if (!emission) error=102;
    }
    if (bestfit)
    { bestfit=new float [ncurve*nfit];
      if (!bestfit) error=102;
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} //end of Spectrum::init

int Spectrum::getlength()
{ int ka;
  if (error==0) ka=sizeof(int)+sizeof(Spectrum)+
    commentsize*sizeof(char)+
    (2*mpoint+ncurve*npeak*8+ncurve*nfit+256*8)*sizeof(float);
  else ka=0;
  return(ka);
} //end of Spectrum::getlength

int Spectrum::paexport(char *dest,int length)
{ int ka,kb;
  if ((error==0)&&(dest))
  { ka=0; kb=getlength();
    if (length==0) length=kb;
    else if (length!=kb) ka=-1;
    kb=sizeof(int);
    if (ka>=0) ka=memorysend(dest,(char *)&length,ka,kb,length);
    kb=sizeof(Spectrum);
    if (ka>=0) ka=memorysend(dest,(char *)this,ka,kb,length);
    kb=mpoint*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)ebind,ka,kb,length);
    kb=mpoint*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)xpscount,ka,kb,length);
    kb=256*8*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)curvepar,ka,kb,length);
    kb=commentsize;
    if (ka>=0) ka=memorysend(dest,(char *)comment,ka,kb,length);
    kb=ncurve*numint*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)emission,ka,kb,length);
    kb=ncurve*nfit*sizeof(float);
    if (ka>=0) ka=memorysend(dest,(char *)bestfit,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Spectrum::paexport

int Spectrum::paimport(char *source,int length)
{ int ka,kb;

  if ((error==103)&&(source))
  { memorycopy((char *)&kb,source,sizeof(int));
    ka=sizeof(int);
    if ((length==0)&&(kb>0)) length=kb;
    else if ((length<0)||(kb<0)||(length!=kb)) ka=-1;
    kb=sizeof(Spectrum);
    if (ka>=0) ka=memoryrec((char *)this,source,ka,kb,length);
    if (ka<0) error=901;
    else if (error==0)
    { error=106; init();
    }
    else error=754;
  }
  else ka=-1;
  if (error==0)
  { kb=mpoint*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)ebind,source,ka,kb,length);
    kb=mpoint*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)xpscount,source,ka,kb,length);
    kb=256*8*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)curvepar,source,ka,kb,length);
    kb=commentsize;
    if (ka>=0) ka=memoryrec((char *)comment,source,ka,kb,length);
    kb=ncurve*numint*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)emission,source,ka,kb,length);
    kb=ncurve*nfit*sizeof(float);
    if (ka>=0) ka=memoryrec((char *)bestfit,source,ka,kb,length);
    if ((error==0)&&((ka<0)||(ka!=length))) error=901;
  }
  return(error);
} //end of Spectrum::paimport

int Spectrum::getnumpoint()
{ int k;
  if (error==0) k=npoint;
  else k=0;
  return(k);
} //end of Spectrum::getnumpoint

int Spectrum::getnumcurve()
{ int k;
  if (error==0) k=ncurve;
  else k=0;
  return(k);
} //end of Spectrum::getnumcurve

float Spectrum::getmemory()
{ float xa;
  if (error==0) xa=(float)(sizeof(Spectrum)+commentsize*sizeof(char)+
    (mpoint*2+ncurve*npeak*8+ncurve*nfit+256*8)*sizeof(float));
  else xa=-1.0f;
  return(xa);
} //end of Spectrum::getmemory

/*
----------------------------------------------------------------------
curvepar  0        1       2      3        4
          ephoton  theta   phi    memadd   memend+1

          5              6-7
          reliability    reserved
----------------------------------------------------------------------
*/
int Spectrum::loadcurve(char *filename)
{ int i,j,ka,kb,kc,kd,ke,nregion,ncycle,nrepeat,begrow,ndata,
    linenum,memsize;
  float xa,xb,estart,estop,estep,xpstime,ephoton,ework,atheta,aphi;
  double dya;
  float *tep,*txc;
  char *commentbuf,*membuf;
  FILE *fhandle;

  commentbuf=membuf=NULL; tep=txc=NULL; memsize=256; ndata=0;
  if ((error==101)||(error==0))
  { commentbuf=new char [100*80]; membuf=new char [memsize];
    tep=new float [1024]; txc=new float [1024];
    if ((!commentbuf)||(!membuf)||(!tep)||(!txc)) error=102;
    else error=0;
  }

  if ((error==0)&&(ncurve<256))
  { datatype=841;
    i=begrow=ndata=nregion=ncycle=nrepeat=xpsmode=0;
    ephoton=ework=0.0f;
    std::ifstream fin(filename,std::ios::in);
    if (!fin) error=201+(datatype<<16);
    else
    { fin >> membuf[0];
      if ((fin)&&((membuf[0]==' ')||(membuf[0]=='\t')||
        ((membuf[0]>='0')&&(membuf[0]<='9'))))
      { fin >> xa >> xb;
        if ((membuf[0]>='0')&&(membuf[0]<='9'))
          xa+=(float)((membuf[0]-'0')*100+0.001);
        if (fin) ka=(int)xa;
        else ka=0;
        if (ka==datatype)  //read MSCD data format
        { begrow=(int)xb; ephoton=atheta=aphi=0.0f;
          fin >> xa; skipendline(fin);
          linenum=(int)xa;
          if (linenum<0) linenum=0;
          else if (linenum>256) linenum=256;
          kb=0; commentbuf[0]='\0';
          for (i=2;i<begrow;++i)
          { getendline(fin,membuf,memsize);
            if ((i<50)||(i>begrow-45))
            { kb+=stringlength(membuf)+3;
              if (kb<100*80) kb=stringappend(commentbuf,membuf,0,"\n");
            }
            else if (i==50)
            { kb+=40;
              if (kb<100*80) kb=stringappend(commentbuf,
                "\n  ... more ... more ... more ...\n\n");
            }
          }
          commentbuf[100*80-1]='\0';
          if (!fin) error=207+(datatype<<16);
          else if ((npoint==0)&&(linenum>0))
            error=loadcomment(commentbuf,100*80);
          if ((error==0)&&(linenum>0))
          { ka=1; 
            ndata=linenum; kb=0;
            kb=nextnonwhite(membuf,kb);
            kb=nextwhite(membuf,kb);
            kb=nextnonwhite(membuf,kb);

            kb=nextwhite(membuf,kb);
            kb=nextnonwhite(membuf,kb);
            ephoton=stringtofloat(membuf+kb);

            kb=nextwhite(membuf,kb);
            kb=nextnonwhite(membuf,kb);
            atheta=stringtofloat(membuf+kb);

            kb=nextwhite(membuf,kb);
            kb=nextnonwhite(membuf,kb);
            aphi=stringtofloat(membuf+kb);

            if ((ephoton<1.0)||(ephoton>10000.0))
              ephoton=atheta=aphi=0.0f;
          }
          else if (error==0)
          { fin >> xa; ka=(int)xa;
            skipendline(fin);
          }
          else ka=0;
          if (ncurve+ka>256) ka=256-ncurve;
          for (i=0;(error==0)&&(i<ka);++i)
          { if (linenum==0)
            { fin >> xa >> xb >> ephoton >> atheta >> aphi;
              ndata=(int)xb;
              skipendline(fin);
            }
            if (npoint+ndata>mpoint) ndata=mpoint-npoint;
            curvepar[(ncurve+i)*8]=ephoton;
            curvepar[(ncurve+i)*8+1]=atheta;
            curvepar[(ncurve+i)*8+2]=aphi;
            curvepar[(ncurve+i)*8+3]=(float)npoint;
            curvepar[(ncurve+i)*8+4]=(float)(npoint+ndata);
            for (j=0;(fin)&&(j<ndata);++j)
            { fin >> xa >> xb;
              ebind[npoint+j]=xa;
              xpscount[npoint+j]=xb;
              skipendline(fin);
              if ((j>0)&&(ebind[npoint+j]-ebind[npoint+j-1]<1.0e-5))
                error=207+(datatype<<16);
            }
            npoint+=ndata;
            if (!fin) error=204+(datatype<<16);
            else if (npoint>=mpoint) i=ka;
          }
          ncurve+=ka;
        }
      }

      else if ((fin)&&(membuf[0]==';'))  //read Eli's data format
      { ke=1; ncycle=nregion=0; ephoton=ework=0.0f;
        while ((fin)&&(ke!=0))
        { getendline(fin,membuf,memsize);
          if (fin)
          { ka=stringcomp(membuf,";  # cycles",12);
            kb=stringcomp(membuf,";  # regions",13);
            kc=stringcomp(membuf,";  initial photon energy",25);
            kd=stringcomp(membuf,";  work function",17);
            ke=stringcomp(membuf,"; region",9);
          }
          else ka=kb=kc=kd=ke=-1;
          if (ka==0) j=12;
          else if (kb==0) j=13;
          else if (kc==0) j=25;
          else if (kd==0) j=17;
          else j=0;
          if ((ka==0)||(kb==0)||(kc==0)||(kd==0))
          { xa=stringtofloat(membuf+j);
            if (ka==0) ncycle=(int)xa;
            else if (kb==0) nregion=(int)xa;
            else if (kc==0) ephoton=xa;
            else ework=xa;
          }
        }
        ndata=0;
        for (i=0;i<nregion;++i)
        { ke=1; kc=0;
          while ((fin)&&(ke!=0))
          { getendline(fin,membuf,memsize);
            if (fin)
            { ka=stringcomp(membuf,";  numpts",10);
              ke=stringcomp(membuf,"x",2);
            }
            else ka=ke=-1;
            if (ka==0) j=10;
            else j=0;
            if (ka==0)
            { xa=stringtofloat(membuf+j);
              kc=(int)xa;
            }
          }
          if (ndata+kc>1024) kc=1024-ndata;
          if ((ncycle>=1)&&(ncycle<=100)) xb=1.0f/(float)ncycle;
          else xb=1.0f;
          for (j=0;j<kc;++j)
          { fin >> xa;
            tep[ndata+j]=xa;
          }
          skipendline(fin);
          skipendline(fin);
          for (j=0;j<kc;++j)
          { fin >> xa;
            txc[ndata+j]=xa*xb;
          }
          ndata+=kc;
        }
        if (!fin) error=207+(datatype<<16);
        else
        { for (j=0;j<ndata/2;++j)
          { xa=tep[j]; tep[j]=tep[ndata-j-1]; tep[ndata-j-1]=xa;
            xa=txc[j]; txc[j]=txc[ndata-j-1]; txc[ndata-j-1]=xa;
          }
          for (j=0;j<ndata-1;++j)
          { if (tep[j+1]-tep[j]<1.0e-5)
            { for (kb=j+1;kb<ndata-1;++kb)
              { tep[kb]=tep[kb+1]; txc[kb]=txc[kb+1];
              }
              --j; --ndata;
            }
          }
          if (ndata<5) error=614+(datatype<<16);
        }

        if ((error==0)&&(ndata>250))
        { Polation aspline; aspline.init(ndata);
          if (error==0) error=aspline.loadcurve(tep,txc);
          if (error==0)
          { ka=ndata; ndata=250;
            xb=(tep[ka-1]-tep[0])/float(ndata-1);
            for (j=0;j<ndata;++j)
            { xa=tep[0]+(float)j*xb;
              tep[j]=xa;
              txc[j]=aspline.fspline(xa);
            }
          }
        }
        if (error==0)
        { if (npoint+ndata>mpoint) ndata=mpoint-npoint;
          for (j=0;j<ndata;++j)
          { ebind[npoint+j]=tep[j]; xpscount[npoint+j]=txc[j];
          }
          curvepar[ncurve*8]=ephoton;
          curvepar[ncurve*8+1]=0.0f;
          curvepar[ncurve*8+2]=0.0f;
          curvepar[ncurve*8+3]=(float)npoint;
          curvepar[ncurve*8+4]=(float)(npoint+ndata);
          ++ncurve; npoint+=ndata;
        }
      }

      else if ((fin)&&(membuf[0]=='F'))  //read Phi's ascii format
      { ke=1; ncycle=nregion=0; ephoton=ework=0.0f;
        while ((fin)&&(ke!=0))
        { getendline(fin,membuf,memsize);
          if (fin)
          { ka=stringcomp(membuf,"anodee =",9);
            kb=stringcomp(membuf,"nocycles",9);
            kc=stringcomp(membuf,"noregions",10);
            ke=stringcomp(membuf,"regionacqtimes",15);
          }
          else ka=kb=kc=ke=-1;
          if (ka==0) j=9;
          else if (kb==0) j=10;
          else if (kc==0) j=11;
          else j=0;
          if ((ka==0)||(kb==0)||(kc==0))
          { xa=stringtofloat(membuf+j);
            if (ka==0) ephoton=xa;
            else if (kb==0) ncycle=(int)xa;
            else nregion=(int)xa;
          }
        }
        ndata=0;
        for (i=0;i<nregion;++i)
        { ke=1; kc=0;
          fin >> xpstime;
          while ((fin)&&(ke!=0))
          { getendline(fin,membuf,memsize);
            if (fin) ke=stringcomp(membuf,"xr",3);
          }
          if ((fin)&&(ke==0))
          { j=0;
            while (membuf[j]!='\0')
            { if ((membuf[j]=='[')||(membuf[j]==']'))
                membuf[j]=' ';
              ++j;
            }
            j=4;
            xa=stringtofloat(membuf+j);
            kc=(int)xa;
          }
          if (ndata+kc>1024) kc=1024-ndata;
          if ((ncycle>=1)&&(ncycle<=100)&&(xpstime>1.0e-5))
            xb=0.001f/(float)ncycle/xpstime;
          else xb=1.0f;
          for (j=0;j<kc;++j)
          { fin >> xa;
            tep[ndata+j]=xa;
          }
          skipendline(fin);
          skipendline(fin);
          for (j=0;j<kc;++j)
          { fin >> xa;
            txc[ndata+j]=xa*xb;
          }
          ndata+=kc;
        }
        if (!fin) error=207+(datatype<<16);
        else
        { for (j=0;j<ndata/2;++j)
          { xa=tep[j]; tep[j]=tep[ndata-j-1]; tep[ndata-j-1]=xa;
            xa=txc[j]; txc[j]=txc[ndata-j-1]; txc[ndata-j-1]=xa;
          }
          for (j=0;j<ndata-1;++j)
          { if (tep[j+1]-tep[j]<1.0e-5)
            { for (kb=j+1;kb<ndata-1;++kb)
              { tep[kb]=tep[kb+1]; txc[kb]=txc[kb+1];
              }
              --j; --ndata;
            }
          }
          if (ndata<5) error=614+(datatype<<16);
        }

        if ((error==0)&&(ndata>250))
        { Polation aspline; aspline.init(ndata);
          if (error==0) error=aspline.loadcurve(tep,txc);
          if (error==0)
          { ka=ndata; ndata=250;
            xb=(tep[ka-1]-tep[0])/float(ndata-1);
            for (j=0;j<ndata;++j)
            { xa=tep[0]+(float)j*xb;
              tep[j]=xa;
              txc[j]=aspline.fspline(xa);
            }
          }
        }
        if (error==0)
        { if (npoint+ndata>mpoint) ndata=mpoint-npoint;
          for (j=0;j<ndata;++j)
          { ebind[npoint+j]=tep[j]; xpscount[npoint+j]=txc[j];
          }
          curvepar[ncurve*8]=ephoton;
          curvepar[ncurve*8+1]=0.0f;
          curvepar[ncurve*8+2]=0.0f;
          curvepar[ncurve*8+3]=(float)npoint;
          curvepar[ncurve*8+4]=(float)(npoint+ndata);
          ++ncurve; npoint+=ndata;
        }
      }

      else if (fin)     //read Phi's binary format
      { fin.close();
        fhandle=fopen(filename,"rb");
        if (error==0) ka=sizeof(int);
        else ka=sizeof(int);
        if (error==0) kb=sizeof(double);
        else kb=sizeof(double);
        if (!fhandle) error=201+(datatype<<16);
        else if (ka<4) error=105;
        else if (kb!=8) error=108;
        else
        { ka=kb=kc=0;
          fread(&ka,4,1,fhandle);
          if (ka!=0) error=207+(datatype<<16);
          else
          { fread(&ka,4,1,fhandle);
            xpsmode=ka+1;
          }
          if ((error==0)&&((xpsmode<1)||(xpsmode>2)))
            error=207+(datatype<<16);
          else if (error==0)
          { fread(membuf,1,76,fhandle);
            fread(&kb,4,1,fhandle);
            nregion=kb;
            fread(&kb,4,1,fhandle);
            fread(&kb,4,1,fhandle);
            ncycle=kb;
            fread(membuf,1,56,fhandle);
            fread(&dya,8,1,fhandle);
            ework=(float)dya;
            fread(membuf,1,40,fhandle);
            fread(&dya,8,1,fhandle);
            ephoton=(float)dya;
            fread(membuf,1,103,fhandle);
          }
          ndata=0;
          for (i=0;(error==0)&&(fhandle)&&(i<kb);++i)
          { dya=0.0; kc=0;
            fread(&kc,4,1,fhandle);
            fread(&dya,8,1,fhandle);
            estart=(float)dya;
            fread(&dya,8,1,fhandle);
            estop=estart+(float)dya;
            fread(membuf,1,16,fhandle);
            fread(&dya,8,1,fhandle);
            estep=(float)dya;
            fread(&dya,8,1,fhandle);
            xpstime=(float)dya;
            fread(&kc,4,1,fhandle);
            nrepeat=kc;
            xb=(float)(ncycle*xpstime*nrepeat);
            if (xb>1.0e-5) xb=1.0f/xb;
            else xb=1.0f;
            fread(membuf,1,20,fhandle);
            fread(&kc,4,1,fhandle);
            if (ndata+kc>1024) kc=1024-ndata;
            for (j=0;j<kc;++j)
            { fread(&dya,8,1,fhandle);
              xa=estart-estep*(float)j;
              tep[ndata+j]=xa;
              txc[ndata+j]=(float)dya*xb;
            }
            ndata+=kc;
          }
          if (!fhandle) error=207+(datatype<<16);
          else
          { for (j=0;j<ndata/2;++j)
            { xa=tep[j]; tep[j]=tep[ndata-j-1]; tep[ndata-j-1]=xa;
              xa=txc[j]; txc[j]=txc[ndata-j-1]; txc[ndata-j-1]=xa;
            }
            for (j=0;j<ndata-1;++j)
            { if (tep[j+1]-tep[j]<1.0e-5)
              { for (kb=j+1;kb<ndata-1;++kb)
                { tep[kb]=tep[kb+1]; txc[kb]=txc[kb+1];
                }
                --j; --ndata;
              }
            }
            if (ndata<5) error=614+(datatype<<16);
          }
          if ((error==0)&&(ndata>250))
          { Polation aspline; aspline.init(ndata);
            if (error==0) error=aspline.loadcurve(tep,txc);
            if (error==0)
            { ka=ndata; ndata=250;
              xb=(tep[ka-1]-tep[0])/float(ndata-1);
              for (j=0;j<ndata;++j)
              { xa=tep[0]+(float)j*xb;
                tep[j]=xa;
                txc[j]=aspline.fspline(xa);
              }
            }
          }
          if (error==0)
          { if (npoint+ndata>mpoint) ndata=mpoint-npoint;
            for (j=0;j<ndata;++j)
            { ebind[npoint+j]=tep[j]; xpscount[npoint+j]=txc[j];
            }
            curvepar[ncurve*8]=ephoton;
            curvepar[ncurve*8+1]=0.0f;
            curvepar[ncurve*8+2]=0.0f;
            curvepar[ncurve*8+3]=(float)npoint;
            curvepar[ncurve*8+4]=(float)(npoint+ndata);
            ++ncurve; npoint+=ndata;
          }
        }
        fclose(fhandle);
      }
    }
  }

  if (commentbuf) delete [] commentbuf; if (membuf) delete [] membuf;
  if (tep) delete [] tep; if (txc) delete [] txc;
  return(error);
} //end of Spectrum::loadcurve

int Spectrum::clean(int replace)
{ int i,j,k;
  float xa,xb,xc;

  if (error==0)
    for (i=0;i<ncurve;++i)
      for (k=5;k<8;++k)
        curvepar[i*8+k]=0.0f;

  if ((error==0)&&(replace!=0))
  { for (i=0;i<ncurve-1;++i)
    { for (j=i+1;j<ncurve;++j)
      { xa=(float)(fabs(curvepar[i*8]-curvepar[j*8])+
          fabs(curvepar[i*8+1]-curvepar[j*8+1])+
          fabs(curvepar[i*8+2]-curvepar[j*8+2]));
        if (xa<1.0e-5) j=ncurve+10;
      }
      if (j>ncurve)
      { for (j=i;j<ncurve-1;++j)
        { for (k=0;k<8;++k)
          curvepar[j*8+k]=curvepar[(j+1)*8+k];
        }
        --ncurve; --i;
      }
    }
  }
  if (error==0)
  { for (i=0;i<ncurve-1;++i)
    { for (j=i+1;j<ncurve;++j)
      { xa=curvepar[i*8]; xb=curvepar[j*8];
        if (xa==xb)
        { xa=curvepar[i*8+1]; xb=curvepar[j*8+1];
        }
        if (xa==xb)
        { xa=curvepar[i*8+2]; xb=curvepar[j*8+2];
        }
        if (xa>xb)
        { for (k=0;k<8;++k)
          { xc=curvepar[i*8+k];
            curvepar[i*8+k]=curvepar[j*8+k];
            curvepar[j*8+k]=xc;
          }
        }
      }
    }
  }
  if (error==0)
  { npoint=0;
    for (i=0;i<ncurve;++i)
    { k=(int)(curvepar[i*8+4]-curvepar[i*8+3]+0.001);
      npoint+=k;
    }
  }

  return(error);
} //end of Spectrum::clean

int Spectrum::loadcomment(char *commentbuf,int bufsize)
{ int i;

  if ((error==0)&&(ncurve==0))
  { commentsize=stringlength(commentbuf);
    if (commentsize>bufsize) commentsize=bufsize;
    if (comment) delete [] comment;
    comment=new char [commentsize];
    if (!comment) error=102;
    else
    { for (i=0;i<commentsize-1;++i) comment[i]=commentbuf[i];
      if ((comment[commentsize-2]=='\n')||
        (comment[commentsize-2]=='\r')) comment[commentsize-2]='\0';
      comment[commentsize-1]='\0';
    }
  }
  if ((error==0)&&((ncurve>0)||(commentsize<5)))
  { if (comment) delete [] comment;
    commentsize=0; comment=NULL;
  }

  return(error);
} //end of Spectrum::loadcomment

int Spectrum::savecurve(char *filename,char *usermessage)
{ int i,j,k,ndata,begrow,linenum;
  float ephoton,atheta,aphi;

  if (error==0)
  { Textout fileout(filename);
    error=fileout.geterror();
    if (error==0)
    { datatype=841;
      i=0; begrow=0;
      if (ncurve<2) linenum=npoint;
      else linenum=0;
      if ((linenum>0)&&(comment))
      { while (comment[i]!='\0')
        { if ((comment[i]=='\n')||(comment[i]=='\r')) ++begrow;
          ++i;
        }
      }
      else
      { while (usermessage[i]!='\0')
        { if ((usermessage[i]=='\n')||(usermessage[i]=='\r')) ++begrow;
          ++i;
        }
      }
      if ((linenum>0)&&(comment)) begrow+=3;
      else if (linenum>0) begrow+=9;
      else begrow+=7;
      fileout.integer(datatype,6); fileout.integer(begrow,6);
      fileout.integer(linenum,6); fileout.space(5);
      if (linenum>0)
        fileout.string("datakind beginning-row linenumbers",0,1);
      else fileout.string("datakind beginning-row multi-curves",0,1);
      if ((linenum>0)&&(comment)) fileout.string(comment,0,1);
      else
      { fileout.string(usermessage,0,1);
        fileout.string("   X-ray photoemission spectrum",0,1);
        fileout.string("   counts versus photoelectron kinetic energy");
        fileout.string(" (from Fermi level)",0,1);
        fileout.string(
          "     parameters: curve point photon-energy(eV) theta phi",
          0,1);
        fileout.string(
          "     columns: binding-energy (eV) counts (K/sec)",0,1);
        fileout.integer(ncurve,6);
        fileout.integer(npoint,6);
        fileout.string("  ncurve npoint",0,1);
      }
      for (i=0;i<ncurve;++i)
      { ndata=(int)(curvepar[i*8+4]-curvepar[i*8+3]+0.001);
        ephoton=curvepar[i*8];
        atheta=curvepar[i*8+1];
        aphi=curvepar[i*8+2];
        k=(int)curvepar[i*8+3];
        if ((linenum==0)||(!comment))
        { fileout.integer(i+1,4);
          fileout.integer(ndata,6);
          fileout.floating(ephoton,10.2f,256);
          fileout.floating(atheta,8.1f,256);
          fileout.floating(aphi,8.1f,256);
          fileout.string(
            "     ------------------------",0,1);
        }
        for (j=0;j<ndata;++j)
        { fileout.floating(ebind[k+j],10.2f,256);
          fileout.floating(xpscount[k+j],11.3f,256+1);
        }
      }
      if (error==0) error=fileout.geterror();
    }
    else error+=(datatype<<16);
  }

  return(error);
} //end of Spectrum::savecurve

