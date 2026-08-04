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

1. Textout
Description: output number or text into file stream
Constructor and Destructor:
  Textout() --- constructor
  Textout(const char *filename, int mode=ios::out|ios::trunc)
    --- constructor, makes an ofstream,
      opens a file for writing, and connects to it
Member Functions:
  int open(const char *filename,int mode=ios::out|ios::trunc)
    --- open a file for a specific object
  int geterror() --- output error code
  int newline() --- output into the stream for a linefeed
  int newline(int linenum) --- output into the stream for a number
    of linefeeds
  int space(int width,int linenum=0) --- output
    into the stream for a space string with width, and linenum
    if the rightmost bit of linenum is not zero
  int charfill(char ch,int width=0,int linenum=0) --- output
    into the stream for a number of character ch with width, and linenum
    if the rightmost bit of linenum is not zero
  int string(char *string,int width=0,int linenum=0) --- output
    into the stream for a string with width, and linenum if the
    rightmost bit of linenum is not zero
  int integer(int number,int width=0,int linenum=0) --- output
    into the stream for an integer number with width, and
    linenum if the rightmost bit of linenum is not zero
  int floating(float realnumber,float widthprecision=0,int linenum=0)
    --- output into the stream for a floating point number with width
    (integer part of the widthprecision) and precision (fractional
    part of the widthprecision), and linenum if the rightmost bit of
    linenum is not zero.
  int flush() --- flush the output stream

  The bit3-0 of linenum is the number of new lines after the text.
  The bit7-4 of linenum is the number of new lines before the text.
  The bit9-8 of linenum indicate default format if 00, fixed format if
    01, scientific format if 10, fixed format and fixed digits if 11.


List of functions

1.  int iabs(int ix) --- return the absolute integer value
2.  (char, int, float) confine(x,min,max) --- return the confined
      value within the specified minimum and maximum values.
3.  float round(float xa,float decimal) --- return
      ya=decimal*(int)(xa/decimal+0.5)
4.  int spacehide(char *membuf) --- convert spaces within quatation
      marks to double quatation and return the number of bytes of the
      string.
5.  int spaceback(char *membuf) --- convert doulbe quatation marks to
      spaces and return the number of bytes of the string.
6.  int spaceclear(char *membuf) --- delete starting and ending spaces
      and return the number of bytes of the string
7.  void skipdelim(istream& is) --- skip deliminators
8.  int skipendline(istream& is) --- skip input stream until the
      end of line (auto-detect PC, Unix or Macintosh text format)
9.  int getendline(istream& is,char *membuf,int memsize) ---
      read input stream until the end of line and saved into
      memory buffer (membuf) (auto-detect PC, Unix or Macintosh
      text format)
10. int skipbytes(istream& is,int nbyte) --- skil input stream number
      of bytes (nbyte) before the end of line
11. int stringlength(const char *source) --- return the length of the
      string
12. int stringcopy(char *dest,const char *source,int nbyte=0)
      --- copy first n bytes of the source string or full string
      if the length less than nbyte to the dest string
13. int stringappend(char *dest,const char *source,int nbyte=0,
      const char *sourcea=NULL,const char *sourceb=NULL,
      const char *sourcec=NULL,const char *sourced=NULL)
      --- append first nbyte bytes of the source string
      to the destination string
14. int stringcomp(const char *source1,const char *source2,int nbyte)
      --- compare first n bytes of the source1 and source2 strings, or
      one full string if its length less than nbyte, without case
      sensitivity, return 1 if source1>source2, 0 if source1=source2,
      -1 if source1<source2
15. int stringclear(char *membuf) --- clean the memory buffer membuf
      to allow only one space character 0x20 between other ascii
      characters, deliminaters are translated as space character
16. int stringtoint(char *string) --- convert a string into a
      integer number and return it
17. int stringtofloat(char *string) --- convert a string into a
      floating point number and return it
18. int stringmatch(char *source1,char *source2) ---
      return 0 if one string is included in the other
19. int nextchar(char *source,char chx,int begadd=0,int endadd=0,
      int direction=0)
      --- return the address of next character chx
      in forward direction if direction = 0 or backward if
      direction = 1, from begadd to endadd
20. int nextnonchar(char *source,char chx,int begadd=0,int endadd,
      int direction=0)
      --- return the address of next character non-chx
      in forward direction if direction = 0 or backward if
      direction = 1, from begadd to endadd
21. int nextwhite(char *string,int begadd=0,int endadd,
      int direction=0)
      --- return the address of next white character in
      forward direction if direction = 0 or backward if direction = 1,
      from begadd to endadd
22. int nextnonwhite(char *string,int begadd=0,int endadd,
      int direction=0)
      --- return the address of next non-white character in
      forward direction if direction = 0 or backward if direction = 1,
      from begadd to endadd
23. int memorycopy(char *dest,const char *source,int nbyte) --- copy
      nbyte bytes of the source memeory block to the dest memory block
24. int memorysend(char *dest,char *source,int address,int nbyte,
      int destsize) --- copy nbyte bytes of source memory block to
      the destination memory buffer from address (address)
25. int memoryrec(char *dest,char *source,int address,int nbyte,
      int destsize) --- copy nbyte bytes source memory block from
      address (address) to the destination memory buffer
26. void waitenter() --- pause for striking Enter key

Machine dependent functions

27. int computer(char *name=NULL,int size=0) --- copy name of computer
      into name with length limit of size, and return code of the
      computer

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

28. long timeprocessor() --- return current CPU time in seconds
29. int mpiinit(int *argc,char ***argv) --- initialization of
      message passing interface for parallelization
30. int mpigetmype() --- get id or rank of the current processor
      element
31. int mpigetnumpe() --- get total number of processor elements
32. int mpisend(char *buffer,int bufsize,int destpe,int tag) ---
      send out a memory buffer with size of bufsize to the destination
      processor destpe with message tag (tag)
33. int mpireceive(char *buffer,int bufsize,int sourcepe,int tag)
      --- receive from source processor element sourcepe a memory
      buffer with size of bufsize and put in buffer
34. int mpiend() --- finalize the message passing interface
35. int sendemail(char *mailpath,char *message,char *stamp,
      char *fromaddress=NULL,char *toaddress=NULL,char *subject=NULL,
      char *ccaddress=NULL,char *bccaddress=NULL,int stamp=1)
      --- send email from fromaddress to toaddress, cc to ccaddress,
      bcc to bccaddres, with subject, message and stamp if stamp != 0
----------------------------------------------------------------------
*/

#ifndef __USERUTIL_H
#define __USERUTIL_H

#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;

#ifndef NULL
#define NULL 0
#endif

class Textout
{ private:
    int error,fdout;
    ofstream fout;
  public:
    Textout();
    Textout(const char *filename, std::ios_base::openmode mode=std::ios::out | std::ios::trunc);
    int open(const char *filename, std::ios_base::openmode mode=std::ios::out | std::ios::trunc);
    int geterror();
    int newline();
    int newline(int linenum);
    int space(int width,int linenum=0);
    int charfill(char ch,int width=0,int linenum=0);
    int string(char *string,int width=0,int linenum=0);
    int integer(int number,int width=0,int linenum=0);
    int floating(float realnumber,float widthprecision=0,
      int linenum=0);
    int flush();
};

int iabs(int ix);
char confine(char x,char xmin,char xmax);
int confine(int x,int xmin,int xmax);
float confine(float x,float xmin,float xmax);
float round(float xa,float decimal);
int spacehide(char *membuf);
int spaceback(char *membuf);
int spaceclear(char *membuf);
void skipdelim(istream& is);
int skipendline(istream& is);
int getendline(istream& is,char *membuf,int memsize);
int skipbytes(istream& is,int nbyte);
int stringlength(const char *source);
int stringcopy(char *dest,const char *source,int nbyte=0);
int stringappend(char *dest,const char *source,int nbyte=0,
  const char *sourcea=NULL,const char *sourceb=NULL,
  const char *sourcec=NULL,const char *sourced=NULL);
int stringcomp(const char *source1,const char *source2,int nbyte=0);
int stringclear(char *membuf);
int stringtoint(char *source);
float stringtofloat(char *source);
int stringmatch(char *source1,char *source2);
int nextchar(char *source,char chx,int begadd=0,int endadd=0,
  int direction=0);
int nextnonchar(char *source,char chx,int begadd=0,int endadd=0,
  int direction=0);
int nextwhite(char *source,int begadd=0,int endadd=0,
  int direction=0);
int nextnonwhite(char *source,int begadd=0,int endadd=0,
  int direction=0);
int memorycopy(char *dest,const char *source,int nbyte);
int memorysend(char *dest,char *source,int address,int nbyte,
  int destsize);
int memoryrec(char *dest,char *source,int address,int nbyte,
  int destsize);
void waitenter();

//machine dependent functions
int computer(char *name=NULL,int size=0);
long timeprocessor();
int mpiinit(int *argc,char ***argv);
int mpigetmype();
int mpigetnumpe();
int mpisend(char *buffer,int bufsize,int destpe,int tag);
int mpireceive(char *buffer,int bufsize,int sourcepe,int tag);
int mpiend();
int sendemail(char *mailpath,char *message,
  char *fromaddress=NULL,char *toaddress=NULL,char *subject=NULL,
  char *ccaddress=NULL,char *bccaddress=NULL,int stamp=1);

#endif //__USERUTIL_H

