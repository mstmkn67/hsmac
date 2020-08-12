#ifndef _HSMAC3D_BC_H_
#define _HSMAC3D_BC_H_

#include "StaggeredLattice3d.h"
#include <cmath>
#include <iostream>
using namespace std;

class Hsmac3dBC{
public:
	Hsmac3dBC(StaggeredLattice3d* lattice);
	virtual ~Hsmac3dBC();
	virtual void update()=0;
protected:
	StaggeredLattice3d* lattice;
};

class Hsmac3dWall:public Hsmac3dBC{
public:
	Hsmac3dWall(StaggeredLattice3d* lattice,
              double u0,double v0,double w0,
	            double u1,double v1,double w1,
	            double u2,double v2,double w2,
	            double u3,double v3,double w3,
	            double u4,double v4,double w4,
	            double u5,double v5,double w5);
	virtual ~Hsmac3dWall();
	virtual void update();
private:
	double u0,v0,w0,u1,v1,w1,u2,v2,w2,u3,v3,w3,u4,v4,w4,u5,v5,w5;
};

class Hsmac3dChannel:public Hsmac3dBC{
public:
	Hsmac3dChannel(StaggeredLattice3d* lattice,
	               double u0,double v0,double w0,
	               double u2,double v2,double w2,
	               double u4,double v4,double w4,
	               double u5,double v5,double w5,
	               double delp);
	virtual ~Hsmac3dChannel();
	virtual void update();
private:
	double u0,v0,w0,u2,v2,w2,u4,v4,w4,u5,v5,w5;
	double delp;
};

class Hsmac3dLeesEdwards:public Hsmac3dBC{
public:
	Hsmac3dLeesEdwards(StaggeredLattice3d* lattice,
	                   double shear_rate,double dt,int* iteration);
	virtual ~Hsmac3dLeesEdwards();
	virtual void update();
private:
	int* iteration;
	double shear_rate;
	double dt;
};

#endif // _HSMAC3D_BC_H_
