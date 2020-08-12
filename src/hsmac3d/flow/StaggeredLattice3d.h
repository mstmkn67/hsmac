#ifndef _STAGGERED_LATTICE3D_H_
#define _STAGGERED_LATTICE3D_H_

#include "../SimpleArray.h"
#include <iostream>
using namespace std;

class StaggeredLattice3d{
public:
	StaggeredLattice3d(int nx,int ny,int nz,double lx,double ly,double lz,
	                   int nux=-1,int nuy=-1,int nuz=-1,
	                   int nvx=-1,int nvy=-1,int nvz=-1,
	                   int nwx=-1,int nwy=-1,int nwz=-1);
	virtual ~StaggeredLattice3d();
	int nx,ny,nz;
	int nux,nuy,nuz;
	int nvx,nvy,nvz;
	int nwx,nwy,nwz;
	double lx,ly,lz;
	double dx,dy,dz;
	Array3d<double> u0,v0,w0;
	Array3d<double> u,v,w,p;
	Array3d<double> fx,fy,fz;
private:
};

#endif // _STAGGERED_LATTICE3D_H_
