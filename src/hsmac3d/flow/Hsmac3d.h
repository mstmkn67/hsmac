#ifndef _HSMAC3D_H_
#define _HSMAC3D_H_

#include "StaggeredLattice3d.h"
#include "Hsmac3dBC.h"
#include <cmath>
using namespace std;

class Hsmac3d{
public:
	Hsmac3d(StaggeredLattice3d* lattice,Hsmac3dBC* bc,
	        double viscosity,double density,double dt,
	        int pressure_iteration,double convergence,double relaxation_coef);
	virtual ~Hsmac3d();
	virtual void initial();
	virtual void update();
	//
	StaggeredLattice3d* lattice;
	Hsmac3dBC* bc;
protected:
	virtual void advance();
	virtual void calc_velocity();
	virtual void calc_pressure();
private:
	//
	double viscosity;
	double density;
	double dt;
	//
	int pressure_iteration;
	double convergence;
	double relaxation_coef;
};

#endif // _HSMAC3D_H_
