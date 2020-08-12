#ifndef _STAGGERED_LATTICE2D_H_
#define _STAGGERED_LATTICE2D_H_

#include "SimpleArray.h"
#include <iostream>
using namespace std;

class StaggeredLattice2d{
public:
	StaggeredLattice2d(int nx,int ny,double lx,double ly,
	                   int nux=-1,int nuy=-1,int nvx=-1,int nvy=-1);
	virtual ~StaggeredLattice2d();
	int nx,ny;
	int nux,nuy;
	int nvx,nvy;
	double lx,ly;
	double dx,dy;
	Array2d<double> u0,v0;
	Array2d<double> u,v,p;
private:
};

#endif // _STAGGERED_LATTICE2D_H_
