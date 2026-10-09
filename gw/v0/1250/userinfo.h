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

1. Userinfo
Description: provides version, author and copyright information.
Constructors and Destructors:
  Userinfo(float version,char *filename,char *author,char *title)
Member function
  float getmemory() --- get memory it took in unit of bytes

2. Errorinfo
Description: provides error message corresponding to the input error
  code.
Constructors and Destructors:
  Errorinfo(float errorcode)
Member function
  float getmemory() --- get memory it took in bytes

----------------------------------------------------------------------
*/

#ifndef __USERINFO_H
#define __USERINFO_H

#ifndef NULL
#define NULL 0
#endif

class Userinfo
{
  public:
    int mype,numpe,error;
    float version;
    char *author,*usermsg;
  public:
    Userinfo(int mype=0,int numpe=1);
    ~Userinfo();
    void init();
    int loadparameter(float version,char *filename,char *author,
      char *title,int machine=0);
    float getmemory();
};

class Errorinfo
{
  public:
    int error;
    char *errmsg;
  public:
    Errorinfo(int errorcode);
    ~Errorinfo();
    float getmemory();
};

#endif //__USERINFO_H

/*
List of Error Codes
  0       no error occured

  101     data not ready
  102     out of memory
  103     memory not ready
  104     out of allocated memory
  105     32 bit compiling required
  106     memory not paimported
  107     memory not flushed
  108     floating number size mismatch

  201     File not found
  202     File open error
  203     File creation error
  204     Data reading error
  205     Data writing error
  206     End-of-file reading error
  207     Data format error
  208     No data available
  209     Data type error

  211     Singular matrix error

  221     incorrect application type
  222     form contain no data
  223     form submitted with unknown method
  224     required field not submitted
  225     access denied from your web server
  226     invalid email address
  227     invalid url address

  231     invalid input format
  232     invalid output format
  233     no input error
  234     no output error
  235     memory not enough for output

  241     password too simple

  601     too few atoms (no or one)
  602     too many atoms
  603     too many kinds of atoms
  604     emitter not found
  605     too small distance between atoms

  611     no available experimental data
  612     experimental chi data too small
  613     data files not match
  614     too few input data points
  615     too many input data points
  616     too few output data points
  617     too many output data points
  618     no available hologram data

  621     symmetry search error
  622     symmetry order error
  623     subshell name error
  624     binding energy adjustment error
  625     binding energy may be too small
  626     binding energy may be too large
  627     binding energy be too small
  628     binding energy be too large

  631     output file name not found
  632     duplicate output file name
  633     list file open error

  641     too few fitting parameters
  642     too many fitting parameters

  651     atom not found

  701     too few data points
  702     too many data points
  703     too few curves
  704     too many curves
  705     no available data
  706     data too small
  707     data order error

  711     energy order error
  712     energy step too small
  713     too few energy points
  714     too many energy points
  715     too few angle points
  716     too many angle points
  717     zero or negative energy

  751     message sending error
  752     message receiving error
  753     received empty message
  754     received error message

  761     too many processors

  901     programming error occured
  999     unknown error occured


Symbol definition

   symbol       definition

        k       photoemission wave vector in unit of (1/angstrom)
       kx       photoemission wave vector component in x direction
       ky       photoemission wave vector component in y direction
    theta       photoemission detector polar angle in unit of degree
      phi       photoemission detector azimuthal angle in unit of degree
intensity       photoemission intensity
      chi       photoemission chi=(intensity-bakground)/bakground
        l       angular momentum index
       li       initial angular momentum
       lf       final angular momentum ( lf = li-1 and li+1 )
 phase(l)       phase shift of a scattering event for l
    I(lf)       integral of the radial component of the final state (lf)
phase(lf)       phase shift of the final state (lf)


Data file format

  first line    data-kind  begining-row  linenumbers

   data-kind    kind of data file
begining-row    number of first row of the data body
 linenumbers    number of lines or points of the data body


  datakind is a three-digit code, the first digit from right is the
scanning mode, the second is the rotation mechanism, the third is
the chi calculation algorithm.

first digit from right
        1       energy dependent curves or hologram in (k intensity chi)
        2       polar angle dependent curves in (theta intensity chi)
        3       azimuthal angle dependent curves in (phi intensity chi)
        4       solid angle dependent hologram in (theta phi intensity
                  chi)
        5-8     same as 1-4, but normalize chi

second digit from right
        1       theta and phi rotations operate on analyzer
        2       theta and phi rotations operate on sample
        3       theta rotation operate on analyzer, phi rotation on
                  sample

third digit from right
        1       calculated photoelectron diffraction pattern, chi
                  calculated by theory
                  (chi=(intensity-reference)/reference)
        2       calculated photoelectron diffraction pattern, chi
                  calculated by theory
                  (chi=(intensity-background)/background)
        3       experimental photoelectron diffraction pattern, chi
                  calculated by experiment
                  (chi=(intensity-background)/background)

        4       same as 3, but no intensity data

        5-6     reserved

other data types

   711  phase shift data in (k phases(l=0 1 2 ...))

   721  radial matrix data in
         (k R(li+1) phase(li+1) R(li-1) phase(li-1))

   731  calculation report

   741  input data file for photoemission calculation

   751  a batch of input files (filename should be scatin.txt)

   811  potential data in (r (angstrom) rV (eV-angstrom))

   821  input file for phase shift or radial matrix calculation
        (psrmin.txt)
   831  subshell eigen wave function data in (r (angstrom) u (arb unit))

   841  XPS spectrum data file in (energy(eV) count)
   842  XPS spectrum peakfit positions and intensities in
        (Eb (eV) counts background peak1, peak2, peak3)

   851  XPS spectrum peak intensities versus kinetic energies for
        all peaks in (k intensity ... reliability)
   852  XPS spectrum peak intensities versus polar angles for
        all peaks in (angle intensity ... reliability)
   853  XPS spectrum peak intensities versus azimuthal angles for
        all peaks in (angle intensity ... reliability)

   911  real space hologram image file in (x y z (angs) intensity
        (arb unit))
   921  input file for hologram transformation

   931  effective scattering factor as function of energy and scattering
        angle
   941  electron inelastic mean free path data as function of energy

   951  thermal vibrational mean square relative displacement data

   961  expectation energy data in (eigenvalue expectation (eV))
   962  gaussian output
   963  relaxation energy data in (eigenvalue relaxation (eV))
   964  roots and weight factors of rys polynomial

*/
