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

# ---- port de GPU, Fase 1 (PLANO_CUDA.md) ----
# O .cu e' compilado sozinho pelo nvcc: nenhum cabecalho do programa entra
# nele, porque fcomplex.h:46 (o friend com argumento padrao que obriga o
# -fpermissive) nao passa pelo front-end de device. Os .o do C++ tem de ser
# construidos com -DMSCDGPU, entao o alvo exige rm -f *.o antes:
#   rm -f *.o && make randmscd_gpu \
#     CPPFLAGS="-O3 -std=c++98 -w -fpermissive -fopenmp -DMSCDGPU"
#   MSCD_GPU=validate mpirun --use-hwthread-cpus -np 1 randmscd_gpu Cov0.txt
NVCC     = nvcc
# -fmad=false NAO e' opcional. O nvcc contrai a*b+c em FMA por padrao; o g++
# aqui compila para x86-64 base, sem FMA, e nao contrai. A diferenca e' de 1
# ulp, mas beta (indice de rotmata) e o k de fexpix sao INDICES INTEIROS
# tirados de um float: 1 ulp perto da fronteira troca a entrada da tabela.
# Medido: com FMA, 80.209 elementos acima de 1e-3 e maxrel 6,2.
NVFLAGS  = -O3 -arch=sm_89 -lineinfo -fmad=false
CUDALIBS = -lcudart
gpuexe   = randmscd_gpu

mscdgpu.o : mscdgpu.cu mscdgpu.h
	$(NVCC) $(NVFLAGS) -c mscdgpu.cu -o mscdgpu.o

$(gpuexe) : $(mscdobj) mscdgpu.o
	$(LINK.cc) $(mscdobj) mscdgpu.o -o $(gpuexe) $(CUDALIBS)

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

