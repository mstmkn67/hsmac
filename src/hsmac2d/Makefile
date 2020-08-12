CC = g++
CFLAG = -c -O3 -I $(PF_FILES)/include
CCFLAG = -o hsmac2d  -L $(PF_FILES)/lib/linux_64 -lplatform
OBJECT=Timer.o gourmain.o cmdlineopts.o \
	BoundaryCondition.o StaggeredLattice2d.o Hsmac2d.o hsmac2d_driver.o

all :$(OBJECT)
	$(CC) $(OBJECT) $(CCFLAG)

Timer.o : Timer.cpp Timer.h
	$(CC) Timer.cpp $(CFLAG)

gourmain.o : udf/gourmain.cpp udf/gourmain.h
	$(CC) udf/gourmain.cpp $(CFLAG)

cmdlineopts.o : udf/cmdlineopts.cpp udf/cmdlineopts.h
	$(CC) udf/cmdlineopts.cpp $(CFLAG)

StaggeredLattice2d.o : StaggeredLattice2d.cpp StaggeredLattice2d.h
	$(CC) StaggeredLattice2d.cpp $(CFLAG)

BoundaryCondition.o : BoundaryCondition.cpp BoundaryCondition.h
	$(CC) BoundaryCondition.cpp $(CFLAG)

Hsmac2d.o : Hsmac2d.cpp Hsmac2d.h
	$(CC) Hsmac2d.cpp $(CFLAG)

hsmac2d_driver.o : hsmac2d_driver.cpp
	$(CC) hsmac2d_driver.cpp $(CFLAG)

clean : 
	rm -f *.o hsmac2d

install :
	strip hsmac2d
	cp hsmac2d ../../bin
