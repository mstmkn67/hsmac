CC = g++
CFLAG = -c -O3 -I $(PF_FILES)/include
CCFLAG = -o hsmac3d  -L $(PF_FILES)/lib/linux_64 -lplatform
OBJECT=Timer.o gourmain.o cmdlineopts.o \
	Hsmac3dBC.o StaggeredLattice3d.o Hsmac3d.o \
	Hsmac3dSimulator.o hsmac3d_driver.o

all :$(OBJECT)
	$(CC) $(OBJECT) $(CCFLAG)

Timer.o : Timer.cpp Timer.h
	$(CC) Timer.cpp $(CFLAG)

gourmain.o : udf/gourmain.cpp udf/gourmain.h
	$(CC) udf/gourmain.cpp $(CFLAG)

cmdlineopts.o : udf/cmdlineopts.cpp udf/cmdlineopts.h
	$(CC) udf/cmdlineopts.cpp $(CFLAG)

StaggeredLattice3d.o : flow/StaggeredLattice3d.cpp flow/StaggeredLattice3d.h
	$(CC) flow/StaggeredLattice3d.cpp $(CFLAG)

Hsmac3dBC.o : flow/Hsmac3dBC.cpp flow/Hsmac3dBC.h
	$(CC) flow/Hsmac3dBC.cpp $(CFLAG)

Hsmac3d.o : flow/Hsmac3d.cpp flow/Hsmac3d.h
	$(CC) flow/Hsmac3d.cpp $(CFLAG)

Hsmac3dSimulator.o : Hsmac3dSimulator.cpp Hsmac3dSimulator.h
	$(CC) Hsmac3dSimulator.cpp $(CFLAG)

hsmac3d_driver.o : hsmac3d_driver.cpp
	$(CC) hsmac3d_driver.cpp $(CFLAG)

clean : 
	rm -f *.o hsmac3d

install :
	strip hsmac3d
	cp hsmac3d ../../bin
