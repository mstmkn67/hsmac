#include "StaggeredLattice2d.h"

StaggeredLattice2d::StaggeredLattice2d(int _nx,int _ny,double _lx,double _ly,
                                       int _nux,int _nuy,int _nvx,int _nvy)
:nx(_nx),ny(_ny),lx(_lx),ly(_ly),nux(_nux),nuy(_nuy),nvx(_nvx),nvy(_nvy),
 u0(0,nx,0,ny+1),v0(0,nx+1,0,ny),
 u(0,nx,0,ny+1),v(0,nx+1,0,ny),p(0,nx+1,0,ny+1){
	dx=lx/nx;dy=ly/ny;
	if(nux>0){
		u.setBounds(0,nux+1,0,nuy+1);
		u0.setBounds(0,nux+1,0,nuy+1);
	}else{
		nux=nx-1;nuy=ny;
	}
	if(nvx>0){
		v.setBounds(0,nvx+1,0,nvy+1);
		v0.setBounds(0,nvx+1,0,nvy+1);
	}else{
		nvx=nx;nvy=ny-1;
	}
}

StaggeredLattice2d::~StaggeredLattice2d(){
}
