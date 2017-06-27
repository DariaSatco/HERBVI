# -*- Makefile -*- (for emacs)

#
# This Makefile is intended for compiling Herwig++ plugins
# You can find plugins here: INSERT URL
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#		

# location of include files
RIVETPATH=/afs/cern.ch/user/d/dsatco/cern_project/local
HEPMCPATH=/afs/cern.ch/sw/lcg/releases/LCG_85/HepMC/2.06.09/x86_64-slc6-gcc49-opt
LHAPDFPATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_85/MCGenerators/lhapdf/6.1.6/x86_64-slc6-gcc49-opt
FASTJETPATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_85/fastjet/3.2.0/x86_64-slc6-gcc49-opt
BOOSTPATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_85/Boost/1.61.0/x86_64-slc6-gcc49-opt
F77=gfortran
nproc=8

INCDIR=`pwd`/include

RIVETINCLUDE= -I$(RIVETPATH)/include -I$(HEPMCPATH)/include -I$(LHAPDFPATH)/include -I$(FASTJETPATH)/include  -I$()/include  -I$(BOOSTPATH)/include 
RIVETLIB= -L$(RIVETPATH)/lib  -L$(HEPMCPATH)/lib  -L$(LHAPDFPATH)/lib  -L$(FASTJETPATH)/lib  -L$(BOOSTPATH)/lib  -lHepMC -lHepMCfio -lfastjet -lsiscone_spherical -lsiscone -lRivet -lLHAPDF
# LIBDIR=`pwd`/../lib
#RIVETINCLUDE  = -I$(INCDIR) -I${boost}/include

#$(HEPMCINCLUDE) $(HEPMCINCLUDE) 

# INCLUDE = $(RIVETINCLUDE) 
# RIVETLIB = -L$(LIBDIR) -lRivet -lHepMC -lHepMCfio
#
# C++ flags
# 
CXX=g++
DEBUG=-g
CXXFLAGS= -O2 -Wall -std=c++11 $(DEBUG)

ALLCCFILES=$(shell echo *.cc)

default: herwig65.exe

#herwig6521.inc:
#	wget http://www.hep.phy.cam.ac.uk/theory/webber/Herwig/herwig6521.inc

#HERWIG65.INC:
#	wget http://www.hep.phy.cam.ac.uk/theory/webber/Herwig/HERWIG65.INC

#herwig6521.f:
#	wget http://www.hep.phy.cam.ac.uk/theory/webber/Herwig/herwig6521.f

include: HERWIG65.INC herwig6521.f herwig6521.inc
	mkdir include
	mkdir include/herwig6510;
	cp HERWIG65.INC herwig6521.inc  include/herwig6510;
	cp HERWIG65.INC  include/herwig6510/herwig65.inc

lib/libjimmy.so: include
	wget www.hepforge.org/archive/jimmy/jimmy-4.31.tar.gz;
	tar xzf jimmy-4.31.tar.gz;
	cd  jimmy-4.31; ./configure --prefix=`pwd`/.. --with-herwig6521=`pwd`/..; make -j ${nproc} install
herwig: jimmy


herwig6521.o: herwig6521.f herwig6521.inc include 
	gfortran -c -g herwig6521.f

herwigmain.o: herwigmain.f herwig6521.inc include lib/libjimmy.so
	gfortran -c -g herwigmain.f -I$(INCDIR)/jimmy

	gfortran -c -g -fno-automatic herbvi.f

herwig65.exe: driver.cc herwigmain.o herwig6521.o lib/libjimmy.so
	g++ $(CXXFLAGS) driver.cc herwig6521.o herbvi.o herwigmain.o -lgfortran -lgfortranbegin $(RIVETINCLUDE) $(RIVETLIB) -Llib -ljimmy -o herwig65.exe

clean:
	rm -rf lib include herwig65.exe jimmy* *.o herwig.yoda

