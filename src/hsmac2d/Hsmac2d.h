#ifndef _HSMAC2D_H_
#define _HSMAC2D_H_

#include "udfmanager.h"
#include "StaggeredLattice2d.h"
#include "BoundaryCondition.h"
#include <cmath>
using namespace std;

class Hsmac2d{
public:
	Hsmac2d(UDFManager* in,UDFManager* out);
	virtual ~Hsmac2d();
	virtual void update();
protected:
	virtual void input();
	virtual void initial();
	virtual void advance();
	virtual void calc_velocity();
	virtual void calc_pressure();
	virtual void output();
private:
	UDFManager *in,*out;
	StaggeredLattice2d* lattice;
	BoundaryCondition* bc;
	int iteration;
	//
	double viscosity;
	double density;
	double dt;
	//
	int pressure_iteration;
	double convergence;
	double relaxation_coef;
};

#endif // _HSMAC2D_H_
