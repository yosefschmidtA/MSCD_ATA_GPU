SUFFIXES = .cpp
.SUFFIXES: $(SUFFIXES)
MPIFLAGS =
MPILIBS  = -lmpi -llam
COMPILE.cc = mpic++ -c $(CCFLAGS) $(CPPFLAGS) $(MPIFLAGS)
LINK.cc = mpic++ $(CCFLAGS) $(CPPFLAGS)
.cpp:
	$(LINK.cc)  -o $@ $< $(LDLIBS)
.cpp.o:
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
        
CPPFLAGS = -O3 -std=c++98 -Wall
mscdobj   = userinfo.o userutil.o cartesia.o polation.o \
            curvefit.o pdinten.o pdintena.o fcomplex.o \
            msfuncs.o pdchifit.o vibrate.o meanpath.o \
            phase.o radmat.o rotamat.o mscdrun.o \
            mscdruna.o mscdrunb_not_reanalize.o mscdrunc.o mscdrund.o \
            mscdrune.o mscdjob.o mscdmain.o \
            jobtime.o userCluster.o
mscdexe   = randmscd_parallel

calchiobj = userinfo.o userutil.o cartesia.o polation.o \
            curvefit.o pdinten.o pdintena.o \
            calchi.o
calchiexe = calchi

calnoxobj = userinfo.o userutil.o cartesia.o polation.o \
            curvefit.o pdinten.o pdintena.o \
            calnox.o
calnoxexe = calnox

caldifobj = userinfo.o userutil.o cartesia.o polation.o \
            curvefit.o pdinten.o pdintena.o \
            caldif.o
caldifexe = caldif

calfacobj = userinfo.o userutil.o polation.o phase.o \
            msfuncs.o rotamat.o fcomplex.o scatter.o \
            calfac.o
calfacexe = calfac

poconvobj = userinfo.o userutil.o potentia.o poconv.o
poconvexe = poconv

psconvobj = userinfo.o userutil.o fcomplex.o phase.o psconv.o
psconvexe = psconv

rmconvobj = userinfo.o userutil.o radmat.o rmconv.o
rmconvexe = rmconv

calmfpobj = userinfo.o userutil.o meanpath.o calmfp.o
calmfpexe = calmfp

calvibobj = userinfo.o userutil.o vibrate.o calvib.o
calvibexe = calvib

spconvobj = userinfo.o userutil.o polation.o xpspec.o \
            spconv.o
spconvexe = spconv

xpspeakobj = userinfo.o userutil.o polation.o xpspec.o \
            xpspeca.o curvefit.o xpspeak.o
xpspeakexe = xpspeak

mscdall : $(mscdexe) $(calchiexe) $(calnoxexe) $(caldifexe) \
          $(calfacexe) $(poconvexe) $(psconvexe) $(rmconvexe) \
          $(calmfpexe) $(calvibexe) $(spconvexe) $(xpspeakexe)

$(mscdexe) : $(mscdobj)
	$(LINK.cc) $(mscdobj) -o $(mscdexe)

$(calchiexe) : $(calchiobj)
	$(LINK.cc) $(calchiobj) -o $(calchiexe)

$(calnoxexe) : $(calnoxobj)
	$(LINK.cc) $(calnoxobj) -o $(calnoxexe)

$(caldifexe) : $(caldifobj)
	$(LINK.cc) $(caldifobj) -o $(caldifexe)

$(calfacexe) : $(calfacobj)
	$(LINK.cc) $(calfacobj) -o $(calfacexe)

$(poconvexe) : $(poconvobj)
	$(LINK.cc) $(poconvobj) -o $(poconvexe)

$(psconvexe) : $(psconvobj)
	$(LINK.cc) $(psconvobj) -o $(psconvexe)

$(rmconvexe) : $(rmconvobj)
	$(LINK.cc) $(rmconvobj) -o $(rmconvexe)

$(calmfpexe) : $(calmfpobj)
	$(LINK.cc) $(calmfpobj) -o $(calmfpexe)

$(calvibexe) : $(calvibobj)
	$(LINK.cc) $(calvibobj) -o $(calvibexe)

$(spconvexe) : $(spconvobj)
	$(LINK.cc) $(spconvobj) -o $(spconvexe)

$(xpspeakexe) : $(xpspeakobj)
	$(LINK.cc) $(xpspeakobj) -o $(xpspeakexe)

