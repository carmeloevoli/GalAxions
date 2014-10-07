ARCH = $(shell uname -s)
CC = icc # g++
CFLAGS = -O3 -Wall -fno-common -fopenmp  
HEALPIXDIR = #/afs/desy.de/user/c/cevoli/x3dir/libs/Healpix_3.11

INCDIR=. -I$(CFITSIO_DIR)/include -I$(GSL_DIR)/include -I$(HEALPIXDIR)/include # -I$(ROOTSYS)/include 

LIBS =-lm -L/phenod/data/cevoli/libs -lgfortran
LIBS+=-L$(GSL_DIR)/lib -lgsl -lgslcblas
LIBS+=-L$(CFITSIO_DIR)/lib -lcfitsio
#LIBS+=-L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -ldl
#LIBS+=$(HEALPIXDIR)/lib/libchealpix.a 

EXEC  = galAxions.exe	
OBJS  = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
INCS  = $(wildcard *.h)	

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) 

%.o: %.cpp $(INCS)
	$(CC) $(CFLAGS) -I$(INCDIR) -c -o $@ $<

clean:
	rm -f $(EXEC)
	rm  -f *.o
	rm  -f *~	

tar: clean
	tar cfz galaxions.tar.gz *.cpp *.h Makefile --exclude=output/* 

#NEDIR=/afs/desy.de/user/c/cevoli/WORK/DM/AXION/NE2001/src.NE2001
#NEDIR=.
#FITSIOLIB=$(FITSIODIR)/lib
#HEALPIXLIB=$(HEALPIXDIR)/lib
#LIBS=-lm ./libNE2001.a -lgfortran -L/sw/lib -lgsl -lgslcblas -L$(FITSIODIR)/lib -lcfitsio -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,//Users/maccione/software/root/lib -lm -ldl
#$(CC) $(CFLAGS) -o $@ $^ $(HEALPIXDIR)/lib/libchealpix.a $(NEDIR)/libNE2001.a $(LIBS) 
