#include <math.h>

#include "mscdrun.h"
#include "pdchifit.h"

Pdchifit::Pdchifit()
{ error=107; init();
} // end of Pdchifit::Pdchifit

Pdchifit::~Pdchifit()
{ if (xdata) delete [] xdata; if (ydata) delete [] ydata;
  if (invsig) delete [] invsig; if (ymod) delete [] ymod;
  if (afit) delete [] afit; if (safit) delete [] safit;
  if (yafit) delete [] yafit; if (dafit) delete [] dafit;
  if (atry) delete [] atry; if (beta) delete [] beta;
  if (dyda) delete [] dyda;
  if (covar) delete [] covar; if (alpha) delete [] alpha;
} // end of Pdchifit::~Pdchifit

void Pdchifit::init(int indata,int infit)
{ if ((error==107)||(error==103))
  { xdata=ydata=invsig=ymod=afit=safit=NULL;
    yafit=NULL; dafit=atry=beta=dyda=covar=alpha=NULL;
  }
  if (error==103)
  { if (indata<1) error=701;
    else if (indata>4096) error=702;
    else if (infit<0) error=641;
    else if (infit>20) error=642;
    else
    { ndata=indata; nfit=infit; if (ndata<1) nfit=0;
      xdata=new float [ndata]; ydata=new float [ndata];
      invsig=new float [ndata]; ymod=new float [ndata];
      afit=new float [nfit]; safit=new float [nfit];
      yafit=new int [nfit]; dafit=new float [nfit+1];
      atry=new float [nfit]; beta=new float [nfit];
      dyda=new float [nfit*ndata];
      covar=new float [(nfit+1)*nfit]; alpha=new float [nfit*nfit];
      if ((!xdata)||(!ydata)||(!invsig)||(!ymod)||(!afit)||(!safit)||
        (!yafit)||(!dafit)||(!atry)||(!beta)||(!dyda)||(!covar)||
        (!alpha))
      { error=102; ndata=nfit=0;
      }
    }
  }
  if (error==106) error=0;
  else if (error==103) error=101;
  else if (error==107) error=103;

  return;
} // end of Pdchifit::init

float Pdchifit::getmemory()
{ float xa;
  if (error==0)
    xa=(float)(sizeof(Pdchifit)+nfit*4*sizeof(int)+
      (ndata*4+nfit*7+nfit*nfit+nfit*ndata+1)*sizeof(float));
  else xa=-1.0f;
  return(xa);
} //end of Pdchifit::getmemory

int Pdchifit::loadcurve(int ifitmath,int itrymax,float itolerance,
  float *ixdata,float *iydata,float *iafit,float *isafit,
  Mscdrun *imscdrun,
  float (*ifitfunc)(int fitmath,int ndata,int nfit,int *yafit,
    float *afit,float *safit,float *xdata,float *ydata,float *ymod,
    float *dyda,Mscdrun *mscdrun))
{ int i,k;
  float xa;

  if ((error==101)||(error==0))
  { error=0; trynum=0;
    if ((ifitmath<0)||(ifitmath>5)) fitmath=2;
    else fitmath=ifitmath;
    if ((fitmath>0)&&(ndata<5)) error=701;
  }
  if (error==0)
  { if (itrymax>10000) trymax=10000;
    else if (itrymax>0) trymax=itrymax;
    else trymax=0;
    if (itolerance>1.0e-5) tolerance=itolerance;
    else if (itolerance<-1.0e-5) tolerance=-itolerance;
    else tolerance=1.0e-5f;
    xa=0.0f;
    for (i=0;i<ndata;++i)
    { xdata[i]=ixdata[i]; ydata[i]=iydata[i];
      xa+=2.0f*ydata[i]*ydata[i];
    }
    if (xa < 1.0e-30) xa=1.0f;
    for (i=0;i<ndata;++i) invsig[i]=1.0f/xa;
    mfit=0;
    for (i=0;i<nfit;++i)
    { afit[i]=iafit[i]; k=(int)(-isafit[i])-1;
      if (isafit[i]>1.0e-5)
      { yafit[i]=mfit; ++mfit;
        if ((afit[i]>1.0)||(afit[i]<-1.0)) safit[i]=isafit[i]*afit[i];
        else safit[i]=isafit[i];
      }
      else if ((k<0)||(k>=i))
      { yafit[i]=-(nfit+1); safit[i]=0.0f;
      }
      else
      { yafit[i]=-(k+1);
        if ((afit[k]>1.0e-5)||(afit[k]<-1.0e-5))
          safit[i]=afit[i]/afit[k];
        else safit[i]=1.0f;
      }
    }
    mscdrun=imscdrun; fitfunc=ifitfunc;
  }
  return(error);
} // end of Pdchifit::loadcurve

int Pdchifit::dofit(float *ybest,float *abest,float *reliability,
  int *lastrynum)
{ int i;

  if ((error==0)&&(fitmath==5)) error=pdangsearch();
  if ((error==0)&&(fitmath>2)&&(fitmath<5)) error=netsearch();
  if ((error==0)&&(fitmath>1)&&(fitmath<4)) error=downsimplex();
  if ((error==0)&&(fitmath>0)&&(fitmath<4)) error=levenmarq();
  if (error==0)
  { for (i=0;i<nfit;++i) atry[i]=afit[i];
    reliable=(*fitfunc)(0,ndata,nfit,yafit,atry,safit,xdata,ydata,
      ymod,dyda,mscdrun);
    ++trynum; 
    if (reliable>100.0) error=(int)reliable;
    else
    { for (i=0;i<ndata;++i) ybest[i]=ymod[i];
      for (i=0;i<nfit;++i) abest[i]=afit[i];
      *reliability=reliable; *lastrynum=trynum;
    }
  }

  return(error);
} // end of Pdchifit::dofit

int Pdchifit::pdangsearch()
{ int i,ia,ib,ka,kb,kc,devstep;
  float relibest;

  relibest=20.0f; ka=devstep=0;
  if ((error==0)&&(fitmath==5)&&(nfit==2))
  { ka=(int)(sqrt((trymax+2.0)/2.0)+1.0);
    devstep=(int)safit[0]; kc=(int)safit[1];
    if ((devstep>1)&&(devstep<9)&&(ka*devstep>kc)) ka=kc/devstep;
    else if (ka>kc) ka=kc;
  }
  if ((error==0)&&(fitmath==5)&&(nfit==2)&&(ka>1))
  { kb=(trymax-2)/(ka-1);
    if (kb<2) kb=2; else if (kb>20) kb=20;
    kc=(180/(kb-1)+4)/5*5; kb=180/kc+1;
    if (kb<2) kb=2;
    for (ia=0;ia<ka;++ia) for (ib=0;ib<kb;++ib)
    { if ((devstep<9)&&(ia==0)&&(ib>0)) continue;
      if ((devstep>0)&&(devstep<9)) atry[0]=(float)ia*devstep;
      else if (devstep>=9) atry[0]=(float)(ia-(ka-1)/2.0);
      else atry[0]=(float)ia;
      atry[1]=(float)(ib*180.0/(kb-1.0));
      reliable=(*fitfunc)(3,ndata,nfit,yafit,atry,safit,xdata,ydata,
        ymod,dyda,mscdrun);
      ++trynum; if (reliable>100.0) error=(int)reliable;
      if (reliable<relibest)
      { relibest=reliable;
        for (i=0;i<nfit;++i) beta[i]=atry[i];
      }
    }
  }
  if (error==0) for (i=0;i<nfit;++i) afit[i]=beta[i];
  return(error);
} // end of Pdchifit::pdangsearch

int Pdchifit::netsearch()
{ int i,j,k,m,ia,ib,ic,ja,jb,jc,ka,kb,kc;
  float relibest;

  if (fitmath==4) m=trymax;
  else if (fitmath==3) m=trymax*3/4;
  else m=0;

  ja=jb=jc=-1; ka=kb=kc=1;
  for (k=0;(error==0)&&(k<nfit);++k)
  { if ((j=yafit[k])>=0) atry[k]=afit[k];
    else if (j>=-k) atry[k]=afit[-j-1]*safit[k];
    else atry[k]=afit[k];
    beta[k]=atry[k];
    if ((j>=0)&&(ja<0))
    { ja=k; ka=m;
    }
    else if ((j>=0)&&(jb<0))
    { jb=k; ka=kb=(int)sqrt((double)m);
    }
    else if ((j>=0)&&(jc<0))
    { jc=k; ka=kb=kc=(int)pow((double)m,1.0/3.0);
    }
  }
  if ((error==0)&&(jb<0)) jb=ja;
  if ((error==0)&&(jc<0)) jc=ja;

  relibest=20.0f;
  for (ia=0;(error==0)&&(ia<ka);++ia) for (ib=0;ib<kb;++ib)
    for (ic=0;ic<kc;++ic)
  { if (ja>=0) atry[ja]=afit[ja]+(ia-(ka-1)/2)*safit[ja];
    if (jb>ja) atry[jb]=afit[jb]+(ib-(kb-1)/2)*safit[jb];
    if (jc>ja) atry[jc]=afit[jc]+(ic-(kc-1)/2)*safit[jc];
    for (k=0;k<nfit;++k)
    { j=yafit[k];
      if ((ja>=0)&&(ja==-j-1)) atry[k]=atry[ja]*safit[k];
      else if ((jb>ja)&&(jb==-j-1)) atry[k]=atry[jb]*safit[k];
      else if ((jc>ja)&&(jc==-j-1)) atry[k]=atry[jc]*safit[k];
    }
    if ((atry[ja]*afit[ja]>=0.0)&&(atry[jb]*afit[jb]>=0.0)&&
      (atry[jc]*afit[jc]>=0.0))
    { reliable=(*fitfunc)(3,ndata,nfit,yafit,atry,safit,xdata,ydata,
        ymod,dyda,mscdrun);
      ++trynum; if (reliable>100.0) error=(int)reliable;
    }
    else reliable=relibest+1.0f;
    if (reliable<relibest)
    { relibest=reliable;
      for (i=0;i<nfit;++i) beta[i]=atry[i];
    }
  }
  if (error==0) for (i=0;i<nfit;++i) afit[i]=beta[i];
  return(error);
} // end of Pdchifit::netsearch

int Pdchifit::downsimplex()
{ int i,j,k,m,ilow,ihigh,inhigh,itst;
  float rtol,ftol,sum,temp;

  ftol=0.0f;
  if (error==0)
  { ftol=(float)sqrt(fabs(tolerance));
    for (m=0;m<nfit+1;++m)
    { if (m==0) i=0;
      else i=yafit[m-1]+1;
      if ((m==0)||(i>0))
      { for (k=0;k<nfit;++k)
        { j=yafit[k];
          if ((j>=0)&&(m==0)) covar[i*mfit+j]=afit[k];
          else if ((j>=0)&&(j==i-1)) covar[i*mfit+j]=afit[k]-safit[k];
          else if ((j>=0)&&(i>0)) covar[i*mfit+j]=afit[k]+safit[k];
        }
      }
    }
  }
  for (j=0;(error==0)&&(j<mfit);++j)
  { sum=0.0f;
    for (i=0;i<mfit+1;++i) sum+=covar[i*mfit+j];
    beta[j]=sum;
  }
  for (i=0;(error==0)&&(i<mfit+1);++i)
  { for (k=0;k<nfit;++k)
    { if ((j=yafit[k])>=0) atry[k]=covar[i*mfit+j];
      else if (j>=-k) atry[k]=atry[-j-1]*safit[k];
      else atry[k]=afit[k];
    }
    reliable=(*fitfunc)(2,ndata,nfit,yafit,atry,safit,xdata,ydata,
      ymod,dyda,mscdrun);
    dafit[i]=reliable; ++trynum;
    if (reliable>100.0) error=(int)reliable;
  }

  itst=0;
  while ((error==0)&&(itst<mfit+3)&&(trynum<trymax-mfit-2))
  { ilow=ihigh=inhigh=0;
    for (i=0;i<mfit+1;++i)
    { if (dafit[i]<dafit[ilow]) ilow=i;
      else if (dafit[i]>dafit[ihigh])
      { inhigh=ihigh; ihigh=i;
      }
      else if (dafit[i]>dafit[inhigh]) inhigh=i;
    }
    rtol=(float)fabs(dafit[ihigh]-dafit[ilow]);
    if (rtol>1.0e-5)
      rtol=rtol/(float)(fabs(dafit[ihigh])+fabs(dafit[ilow]));
    if (dafit[ilow]<ftol) itst=mfit+3;
    else if (rtol>0.01) itst=0;
    else ++itst;
    if (itst<mfit+3)
    { reliable=downhill(ihigh,-1.0f);
      if (reliable>100.0) error=(int)reliable;
      else if (reliable<=dafit[ilow])
      { reliable=downhill(ihigh,2.0f);
        if (reliable>100.0) error=(int)reliable;
      }
      else if (reliable>=dafit[inhigh])
      { temp=dafit[ihigh];
        reliable=downhill(ihigh,0.5f);
        if (reliable>100.0) error=(int)reliable;
        else if (reliable>=temp)
        { for (i=0;i<mfit+1;++i)
          { if (i!=ilow)
            { for (j=0;j<mfit;++j)
                covar[i*mfit+j]=
                  0.5f*(covar[i*mfit+j]+covar[ilow*mfit+j]);
              for (k=0;k<nfit;++k)
              { if ((j=yafit[k])>=0) atry[k]=covar[i*mfit+j];
                else if (j>=-k) atry[k]=atry[-j-1]*safit[k];
                else atry[k]=afit[k];
              }
              reliable=(*fitfunc)(2,ndata,nfit,yafit,atry,safit,xdata,
                ydata,ymod,dyda,mscdrun);
              dafit[i]=reliable; ++trynum;
              if (reliable>100.0) error=(int)reliable;
            }
          }
          for (j=0;j<mfit;++j)
          { sum=0.0f;
            for (i=0;i<mfit+1;++i) sum+=covar[i*mfit+j];
            beta[j]=sum;
          }
        }
      }
      if (reliable>100.0) error=(int)reliable;
    }
  }
  ilow=0;
  for (i=0;i<mfit+1;++i) if (dafit[i]<dafit[ilow]) ilow=i;
  for (k=0;(error==0)&&(k<nfit);++k)
  { if ((j=yafit[k])>=0) afit[k]=covar[ilow*mfit+j];
    else if (j>=-k) afit[k]=afit[-j-1]*safit[k];
  }
  return(error);
} // end of Pdchifit::downsimplex

int Pdchifit::levenmarq()
{ int itst;
  float pchisq;

  alamda=-1.0f;
  if (error==0) error=marqmin();
  itst=0;
  while ((error==0)&&(itst<mfit+3)&&(trynum<trymax-mfit-1))
  { pchisq=chisq;
    error=marqmin();
    if (pchisq-chisq>0.01*(chisq+pchisq)) itst=0;
    else ++itst;
    if (reliable<tolerance) itst=mfit+3;
  }

  return(error);
} // end of Pdchifit::levenmarq

float Pdchifit::downhill(int ihigh,float factor)
{ int j,k;
  float fac1,fac2;

  if (nfit>0) fac1=(float)((1.0-factor)/nfit);
  else fac1=0.0f;
  fac2=fac1-factor;
  for (k=0;(error==0)&&(k<nfit);++k)
  { if ((j=yafit[k])>=0) atry[k]=beta[j]*fac1-covar[ihigh*mfit+j]*fac2;
    else if (j>=-k) atry[k]=atry[-j-1]*safit[k];
    else atry[k]=afit[k];
  }
  if (error==0)
  { reliable=(*fitfunc)(2,ndata,nfit,yafit,atry,safit,xdata,ydata,
      ymod,dyda,mscdrun);
    ++trynum;
    if (reliable>100.0) error=(int)reliable;
    else if (reliable<dafit[ihigh])
    { dafit[ihigh]=reliable;
      for (k=0;k<nfit;++k)
      { if ((j=yafit[k])>=0)
        { beta[j]+=atry[k]-covar[ihigh*mfit+j];
          covar[ihigh*mfit+j]=atry[k];
        }
      }
    }
  }
  else reliable=(float)error;
  return(reliable);
} // end of Pdchifit::downhill

int Pdchifit::marqmin()
{ int i,j,k,p,q;

  if ((error==0)&&(alamda < 0.0))
  { alamda=0.01f;
    for (i=0;i<nfit;++i)
    { if (((j=yafit[i])>=-i)&&(j<0)) atry[i]=afit[-j-1]*safit[i];
      else atry[i]=afit[i];
    }
    error=marqcof();
    ochisq=chisq+1.0f;
  }
  else if (error==0)
  { for (p=0;p<nfit;++p)
    { if ((j=yafit[p])>=0)
      { for (q=0;q<nfit;++q)
        { if ((k=yafit[q])>=0) covar[j*mfit+k]=alpha[j*mfit+k];
        }
        covar[j*mfit+j]=alpha[j*mfit+j]*(1.0f+alamda);
        dafit[j]=beta[j];
      }
    }
    error=gausjord(covar,dafit,mfit,1);
    if ((error==0) && (alamda!=0.0))
    { for (p=0;p<nfit;++p)
      { if ((j=yafit[p])>=0) atry[p]=afit[p]+dafit[j];
        else if (j>=-p) atry[p]=atry[-j-1]*safit[p];
      }
      error=marqcof();
    }
  }
  if ((error==0) && (alamda>0.0) && (chisq<ochisq))
  { ochisq=chisq;
    if (alamda>1.0e-5) alamda*=0.1f;
    for (p=0;p<nfit;++p)
    { if ((j=yafit[p])>=0)
      { for (q=0;q<nfit;++q)
          if ((k=yafit[q])>=0) alpha[j*mfit+k]=covar[j*mfit+k];
        beta[j]=dafit[j];
        afit[p]=atry[p];
      }
      else if (j>=-p) afit[p]=atry[p];
    }
  }
  else if ((error==0) && (alamda>0.0))
  { chisq=ochisq;
    if (alamda<100.0) alamda*=10.0f;
  }

  return(error);
} // end of Pdchifit::marqmin

int Pdchifit::marqcof()
{ int i,j,k,m,p;
  float wt,dy;

  if (error==0)
  { for (j=0;j<nfit*nfit;++j) covar[j]=0.0f;
    for (j=0;j<nfit;++j) dafit[j]=0.0f;
    chisq=0.0f;
    reliable=(*fitfunc)(1,ndata,nfit,yafit,atry,safit,xdata,ydata,
      ymod,dyda,mscdrun);
    trynum+=mfit+1;
    if (reliable>100.0) error=(int)reliable;
  }
  for (i=0;(error==0)&&(i<ndata);++i)
  { dy=ydata[i]-ymod[i];
    for (p=0;p<nfit;++p)
    { if ((j=yafit[p])>=0)
      { wt=dyda[p*ndata+i]*invsig[i];
        dafit[j]+=dy*wt;
        for (m=0;m<=p;++m)
          if ((k=yafit[m])>=0) covar[j*mfit+k]+=wt*dyda[m*ndata+i];
      }
    }
    chisq+=dy*dy*invsig[i];
  }
  for (j=0;(error==0)&&(j<mfit);++j) for (k=j+1;k<mfit;++k)
    covar[j*mfit+k]=covar[k*mfit+j];
  return(error);
} // end of Pdchifit::marqcof

int Pdchifit::gausjord(float *amat,float *bmat,int na,int nb)
{ int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,p,q;
  float big,dum,pivinv,temp;

  indxc=indxr=ipiv=NULL; icol=irow=0;
  if (error==0)
  { indxc=new int [na]; indxr=new int [na]; ipiv=new int [na];
  }
  if ((error==0) && ((!indxc)||(!indxr)||(!ipiv))) error=102;
  if (error==0)
  { for (j=0;j<na;++j) ipiv[j]=0;
    for (i=0;i<na;++i)
    { big=0.0f;
      for (j=0;j<na;++j)
      { if (ipiv[j]!=1)
        { for (k=0;k<na;++k)
          { if (ipiv[k]==0)
            { if (big<=amat[j*na+k])
              { big=amat[j*na+k]; irow=j; icol=k;
              }
              else if (big<=-amat[j*na+k])
              { big=-amat[j*na+k]; irow=j; icol=k;
              }
            }
            else if (ipiv[k]>1)
            { error=211; i=j=k=na;
            }
          }
        }
      }
      ++ipiv[icol];
      if (irow != icol)
      { for (p=0;p<na;++p)
        { temp=amat[irow*na+p]; amat[irow*na+p]=amat[icol*na+p];
          amat[icol*na+p]=temp;
        }
        for (p=0;p<nb;++p)
        { temp=bmat[irow*nb+p]; bmat[irow*nb+p]=bmat[icol*nb+p];
          bmat[icol*nb+p]=temp;
        }
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if (amat[icol*na+icol]==0.0) error=211;
      else
      { pivinv=1.0f/amat[icol*na+icol];
        amat[icol*na+icol]=1.0f;
        for (p=0;p<na;++p) amat[icol*na+p]*=pivinv;
        for (p=0;p<nb;++p) bmat[icol*nb+p]*=pivinv;
        for (q=0;q<na;++q)
          if (q!=icol)
          { dum=amat[q*na+icol]; amat[q*na+icol]=0.0f;
            for (p=0;p<na;++p) amat[q*na+p]-=amat[icol*na+p]*dum;
            for (p=0;p<nb;++p) bmat[q*nb+p]-=bmat[icol*nb+p]*dum;
          }
      }
    }
    for (p=na-1;p>=0;--p)
    { if (indxr[p]!=indxc[p])
        for (k=0;k<na;++k)
        { temp=amat[k*na+indxr[p]];
          amat[k*na+indxr[p]]=amat[k*na+indxc[p]];
          amat[k*na+indxc[p]]=temp;
        }
    }
  }
  if (indxc) delete[] indxc; if (indxr) delete[] indxr;
  if (ipiv) delete[] ipiv;
  return(error);
} // end of Pdchifit::gausjord

