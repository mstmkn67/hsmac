#ifndef _BOUNDARY_CONDITION_H_
#define _BOUNDARY_CONDITION_H_

#include "StaggeredLattice2d.h"
#include <cmath>
#include <iostream>
using namespace std;

class BoundaryCondition{
public:
	BoundaryCondition(StaggeredLattice2d* lattice);
	virtual ~BoundaryCondition();
	virtual void update()=0;
protected:
	StaggeredLattice2d* lattice;
};

class Wall:public BoundaryCondition{
public:
	Wall(StaggeredLattice2d* lattice,
       double u0,double v0,
	     double u1,double v1,
	     double u2,double v2,
	     double u3,double v3);
	virtual ~Wall();
	virtual void update();
private:
	double u0,v0,u1,v1,u2,v2,u3,v3;
};

class Channel:public BoundaryCondition{
public:
	Channel(StaggeredLattice2d* lattice,
	        double u0,double v0,double u2,double v2,double delp);
	virtual ~Channel();
	virtual void update();
private:
	double u0,v0,u2,v2;
	double delp;
};

class Oscillation:public BoundaryCondition{
public:
	Oscillation(StaggeredLattice2d* lattice,
	            double u2,double omega2,double dt,int* iteration);
	virtual ~Oscillation();
	virtual void update();
private:
	double u2,omega2;
	int* iteration;
	double dt;
};

class LeesEdwards:public BoundaryCondition{
public:
	LeesEdwards(StaggeredLattice2d* lattice,
	            double shear_rate,double dt,int* iteration);
	virtual ~LeesEdwards();
	virtual void update();
private:
	int* iteration;
	double shear_rate;
	double dt;
};

#endif // _BOUDNARY_CONDITION_H_
