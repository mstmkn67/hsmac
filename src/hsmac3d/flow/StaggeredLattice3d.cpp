#include "StaggeredLattice3d.h"

StaggeredLattice3d::StaggeredLattice3d(int _nx,int _ny,int _nz,
                                       double _lx,double _ly,double _lz,
                                       int _nux,int _nuy,int _nuz,
                                       int _nvx,int _nvy,int _nvz,
                                       int _nwx,int _nwy,int _nwz)
:nx(_nx),ny(_ny),nz(_nz),lx(_lx),ly(_ly),lz(_lz),
 nux(_nux),nuy(_nuy),nuz(_nuz),
 nvx(_nvx),nvy(_nvy),nvz(_nvz),
 nwx(_nwx),nwy(_nwy),nwz(_nwz),
 u0(0,nx,0,ny+1,0,nz+1),v0(0,nx+1,0,ny,0,nz+1),w0(0,nx+1,0,ny+1,0,nz),
 u(0,nx,0,ny+1,0,nz+1), v(0,nx+1,0,ny,0,nz+1), w(0,nx+1,0,ny+1,0,nz),
 p(0,nx+1,0,ny+1,0,nz+1),
 fx(0,nx,0,ny+1,0,nz+1),fy(0,nx+1,0,ny,0,nz+1),fz(0,nx+1,0,ny+1,0,nz){
	dx=lx/nx;dy=ly/ny;dz=lz/nz;
	if(nux>0){
		u.setBounds(0,nux+1,0,nuy+1,0,nuz+1);
		u0.setBounds(0,nux+1,0,nuy+1,0,nuz+1);
		fx.setBounds(0,nux+1,0,nuy+1,0,nuz+1);
	}else{
		nux=nx-1;nuy=ny;nuz=nz;
	}
	if(nvx>0){
		v.setBounds(0,nvx+1,0,nvy+1,0,nvz+1);
		v0.setBounds(0,nvx+1,0,nvy+1,0,nvz+1);
		fy.setBounds(0,nvx+1,0,nvy+1,0,nvz+1);
	}else{
		nvx=nx;nvy=ny-1;nvz=nz;
	}
	if(nwx>0){
		w.setBounds(0,nwx+1,0,nwy+1,0,nwz+1);
		w0.setBounds(0,nwx+1,0,nwy+1,0,nwz+1);
		fz.setBounds(0,nwx+1,0,nwy+1,0,nwz+1);
	}else{
		nwx=nx;nwy=ny;nwz=nz-1;
	}
}

StaggeredLattice3d::~StaggeredLattice3d(){
}
