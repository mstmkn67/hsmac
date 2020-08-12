#include "Hsmac3dBC.h"

Hsmac3dBC::Hsmac3dBC(StaggeredLattice3d* l)
:lattice(l){
}

Hsmac3dBC::~Hsmac3dBC(){
}

/////////////////////////////////////////////

Hsmac3dWall::Hsmac3dWall(StaggeredLattice3d* lattice,
           double _u0,double _v0,double _w0,double _u1,double _v1,double _w1,
           double _u2,double _v2,double _w2,double _u3,double _v3,double _w3,
           double _u4,double _v4,double _w4,double _u5,double _v5,double _w5)
:Hsmac3dBC(lattice),u0(_u0),v0(_v0),w0(_w0),u1(_u1),v1(_v1),w1(_w1),
                    u2(_u2),v2(_v2),w2(_w2),u3(_u3),v3(_v3),w3(_w3),
                    u4(_u4),v4(_v4),w4(_w4),u5(_u5),v5(_v5),w5(_w5){
}

Hsmac3dWall::~Hsmac3dWall(){
}

void Hsmac3dWall::update(){
	//down
	for(int i=0;i<=lattice->nux+1;i++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[i][0][k]=2.*u0-lattice->u[i][1][k];
		}
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[i][0][k]=v0;
		}
	}
	for(int i=0;i<=lattice->nwx+1;i++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[i][0][k]=2.*w0-lattice->w[i][1][k];
		}
	}
	//right
	for(int j=1;j<=lattice->nuy;j++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[lattice->nux+1][j][k]=u1;
		}
	}
	for(int j=1;j<=lattice->nvy;j++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[lattice->nvx+1][j][k]=2.*v1-lattice->v[lattice->nvx][j][k];
		}
	}
	for(int j=1;j<=lattice->nwy;j++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[lattice->nwx+1][j][k]=2.*w1-lattice->w[lattice->nwx][j][k];
		}
	}
	//up
	for(int i=0;i<=lattice->nux+1;i++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[i][lattice->nuy+1][k]=2.*u2-lattice->u[i][lattice->nuy][k];
		}
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[i][lattice->nvy+1][k]=v2;
		}
	}
	for(int i=0;i<=lattice->nwx+1;i++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[i][lattice->nwy+1][k]=2.*w2-lattice->w[i][lattice->nwy][k];
		}
	}
	//left
	for(int j=1;j<=lattice->nuy;j++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[0][j][k]=u3;
		}
	}
	for(int j=1;j<=lattice->nvy;j++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[0][j][k]=2.*v3-lattice->v[1][j][k];
		}
	}
	for(int j=1;j<=lattice->nwy;j++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[0][j][k]=2.*w3-lattice->w[1][j][k];
		}
	}
	//front
	for(int i=1;i<=lattice->nux;i++){
		for(int j=1;j<=lattice->nuy;j++){
			lattice->u[i][j][lattice->nuz+1]=2.*u4-lattice->u[i][j][lattice->nuz];
		}
	}
	for(int i=1;i<=lattice->nvx;i++){
		for(int j=1;j<=lattice->nvy;j++){
			lattice->v[i][j][lattice->nvz+1]=2.*v4-lattice->v[i][j][lattice->nvz];
		}
	}
	for(int i=1;i<=lattice->nwx;i++){
		for(int j=1;j<=lattice->nwy;j++){
			lattice->w[i][j][lattice->nwz+1]=w4;
		}
	}
	//back
	for(int i=1;i<=lattice->nux;i++){
		for(int j=1;j<=lattice->nuy;j++){
			lattice->u[i][j][0]=2.*u5-lattice->u[i][j][1];
		}
	}
	for(int i=1;i<=lattice->nvx;i++){
		for(int j=1;j<=lattice->nvy;j++){
			lattice->v[i][j][0]=2.*v5-lattice->v[i][j][1];
		}
	}
	for(int i=1;i<=lattice->nwx;i++){
		for(int j=1;j<=lattice->nwy;j++){
			lattice->w[i][j][0]=w5;
		}
	}
}
////////////////////////////////////////////////////
Hsmac3dChannel::Hsmac3dChannel(StaggeredLattice3d* lattice,
                               double _u0,double _v0,double _w0,
                               double _u2,double _v2,double _w2,
                               double _u4,double _v4,double _w4,
                               double _u5,double _v5,double _w5,
                               double _delp)
:Hsmac3dBC(lattice),u0(_u0),v0(_v0),w0(_w0),u2(_u2),v2(_v2),w2(_w2),u4(_u4),v4(_v4),w4(_w4),u5(_u5),v5(_v5),w5(_w5),delp(_delp){
}
Hsmac3dChannel::~Hsmac3dChannel(){
}

void Hsmac3dChannel::update(){
	//down
	for(int i=0;i<=lattice->nux+1;i++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[i][0][k]=2.*u0-lattice->u[i][1][k];
		}
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[i][0][k]=v0;
		}
	}
	for(int i=0;i<=lattice->nwx+1;i++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[i][0][k]=2.*w0-lattice->w[i][1][k];
		}
	}
	//right
	for(int j=1;j<=lattice->ny;j++){
		for(int k=0;k<=lattice->nz+1;k++){
			lattice->p[lattice->nx+1][j][k]=lattice->p[1][j][k]-delp;
		}
	}
	for(int j=1;j<=lattice->nuy;j++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[lattice->nux+1][j][k]=lattice->u[1][j][k];
		}
	}
	for(int j=1;j<=lattice->nvy;j++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[lattice->nvx+1][j][k]=lattice->v[1][j][k];
		}
	}
	for(int j=1;j<=lattice->nwy;j++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[lattice->nwx+1][j][k]=lattice->w[1][j][k];
		}
	}
	//up
	for(int i=0;i<=lattice->nux+1;i++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[i][lattice->nuy+1][k]=2.*u2-lattice->u[i][lattice->nuy][k];
		}
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[i][lattice->nvy+1][k]=v2;
		}
	}
	for(int i=0;i<=lattice->nwx+1;i++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[i][lattice->nwy+1][k]=2.*w2-lattice->w[i][lattice->nwy][k];
		}
	}
	//left
	for(int j=1;j<=lattice->nuy;j++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[0][j][k]=lattice->u[lattice->nux][j][k];
		}
	}
	for(int j=1;j<=lattice->nvy;j++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[0][j][k]=lattice->v[lattice->nvx][j][k];
		}
	}
	for(int j=1;j<=lattice->nwy;j++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[0][j][k]=lattice->w[lattice->nwx][j][k];
		}
	}
	//front
	for(int i=1;i<=lattice->nux;i++){
		for(int j=1;j<=lattice->nuy;j++){
			lattice->u[i][j][lattice->nuz+1]=2.*u4-lattice->u[i][j][lattice->nuz];
		}
	}
	for(int i=1;i<=lattice->nvx;i++){
		for(int j=1;j<=lattice->nvy;j++){
			lattice->v[i][j][lattice->nvz+1]=2.*v4-lattice->v[i][j][lattice->nvz];
		}
	}
	for(int i=1;i<=lattice->nwx;i++){
		for(int j=1;j<=lattice->nwy;j++){
			lattice->w[i][j][lattice->nwz+1]=w4;
		}
	}
	//back
	for(int i=1;i<=lattice->nux;i++){
		for(int j=1;j<=lattice->nuy;j++){
			lattice->u[i][j][0]=2.*u5-lattice->u[i][j][1];
		}
	}
	for(int i=1;i<=lattice->nvx;i++){
		for(int j=1;j<=lattice->nvy;j++){
			lattice->v[i][j][0]=2.*v5-lattice->v[i][j][1];
		}
	}
	for(int i=1;i<=lattice->nwx;i++){
		for(int j=1;j<=lattice->nwy;j++){
			lattice->w[i][j][0]=w5;
		}
	}
}
////////////////////////////////////////////////////
Hsmac3dLeesEdwards::Hsmac3dLeesEdwards(StaggeredLattice3d* lattice,double s,double d,int* ite)
:Hsmac3dBC(lattice),shear_rate(s),dt(d),iteration(ite){
}

Hsmac3dLeesEdwards::~Hsmac3dLeesEdwards(){
}

void Hsmac3dLeesEdwards::update(){
	//about shear deformation
	// Dn=1
	//  Dx    +---+---+---+---+
	// <----->|   |   |   |   |
	//        +---+---+---+---+
	//   Dnx<>|   |   |   |   |
	// +---+---+---+---+---+-
	// |   |   |   |   |   |
	// +---+---+---+---+---+
	// |   |   |   |   |   |
	// +---+---+---+---+---+
	double a=shear_rate*dt*(*iteration)*lattice->ly;
	double Dx=a-floor(a/lattice->lx)*lattice->lx;//displacement between original and imaginary cells
	int Dn=int(Dx/lattice->dx);//difference of lattice number
	double Dnx=Dx-Dn;
	a=Dnx/lattice->dx;
	//down
	for(int i=0;i<=lattice->nux+1;i++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[i][0][k]=lattice->u[i][lattice->nuy][k]-shear_rate*lattice->ly;
		}
	}
	for(int i=0;i<lattice->nvx+1;i++){
		for(int k=0;k<lattice->nvz+1;k++){
			lattice->v[i][0][k]=lattice->v[i][lattice->nvy][k];
		}
	}
	for(int i=0;i<lattice->nwx+1;i++){
		for(int k=0;k<lattice->nwz+1;k++){
			lattice->w[i][0][k]=lattice->w[i][lattice->nwy][k];
		}
	}
	//right
	for(int j=1;j<=lattice->ny;j++){
		for(int k=0;k<=lattice->nz;k++){
			lattice->p[lattice->nx+1][j][k]=lattice->p[1][j][k];
		}
	}
	for(int j=1;j<=lattice->nuy;j++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[lattice->nux+1][j][k]=lattice->u[1][j][k];
		}
	}
	for(int j=1;j<=lattice->nvy;j++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[lattice->nvx+1][j][k]=lattice->v[1][j][k];
		}
	}
	for(int j=1;j<=lattice->nwy;j++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[lattice->nwx+1][j][k]=lattice->w[1][j][k];
		}
	}
	//up
	for(int i=0;i<=lattice->nx+1;i++){
		for(int k=0;k<=lattice->nz+1;k++){
			int i0=i-Dn-1;int i1=i-Dn;
			if(i0<0)i0+=lattice->nx;if(i1<0)i1+=lattice->nx;
			lattice->p[i][lattice->ny+1][k]=(1.-a)*lattice->p[i0][1][k]+a*lattice->p[i1][1][k];
		}
	}
	for(int i=0;i<=lattice->nux+1;i++){
		for(int k=0;k<=lattice->nuz+1;k++){
			int i0=i-Dn-1;int i1=i-Dn;
			if(i0<0)i0+=lattice->nux;if(i1<0)i1+=lattice->nux;
			lattice->u[i][lattice->nuy+1][k]=(1.-a)*lattice->u[i0][1][k]+a*lattice->u[i1][1][k]+shear_rate*lattice->ly;
		}
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		for(int k=0;k<=lattice->nvz+1;k++){
			int i0=i-Dn-1;int i1=i-Dn;
			if(i0<0)i0+=lattice->nvx;if(i1<0)i1+=lattice->nvx;
			lattice->v[i][lattice->nvy+1][k]=(1.-a)*lattice->v[i0][1][k]+a*lattice->v[i1][1][k];
		}
	}
	for(int i=0;i<=lattice->nwx+1;i++){
		for(int k=0;k<=lattice->nwz+1;k++){
			int i0=i-Dn-1;int i1=i-Dn;
			if(i0<0)i0+=lattice->nwx;if(i1<0)i1+=lattice->nwx;
			lattice->w[i][lattice->nwy+1][k]=(1.-a)*lattice->w[i0][1][k]+a*lattice->w[i1][1][k];
		}
	}
	//left
	//for(int j=1;j<=lattice->ny;j++){
	//	for(int k=0;k<=lattice->nz;k++){
	//		lattice->p[0][j][k]=lattice->p[lattice->nx][j][k];
	//	}
	//}
	for(int j=1;j<=lattice->nuy;j++){
		for(int k=0;k<=lattice->nuz+1;k++){
			lattice->u[0][j][k]=lattice->u[lattice->nux][j][k];
		}
	}
	for(int j=1;j<=lattice->nvy;j++){
		for(int k=0;k<=lattice->nvz+1;k++){
			lattice->v[0][j][k]=lattice->v[lattice->nvx][j][k];
		}
	}
	for(int j=1;j<=lattice->nwy;j++){
		for(int k=0;k<=lattice->nwz+1;k++){
			lattice->w[0][j][k]=lattice->w[lattice->nwx][j][k];
		}
	}
	//back
	//for(int i=1;i<=lattice->nx;i++){
	//	for(int j=1;j<=lattice->ny;j++){
	//		lattice->p[i][j][0]=lattice->p[i][j][lattice->nz+1];
	//	}
	//}
	for(int i=1;i<=lattice->nux;i++){
		for(int j=1;j<=lattice->nuy;j++){
			lattice->u[i][j][0]=lattice->u[i][j][lattice->nuz];
		}
	}
	for(int i=1;i<=lattice->nvx;i++){
		for(int j=1;j<=lattice->nvy;j++){
			lattice->v[i][j][0]=lattice->v[i][j][lattice->nvz];
		}
	}
	for(int i=1;i<=lattice->nwx;i++){
		for(int j=1;j<=lattice->nwy;j++){
			lattice->w[i][j][0]=lattice->w[i][j][lattice->nwz];
		}
	}
	//front
	for(int i=1;i<=lattice->nx;i++){
		for(int j=1;j<=lattice->ny;j++){
			lattice->p[i][j][lattice->nz+1]=lattice->p[i][j][1];
		}
	}
	for(int i=1;i<=lattice->nux;i++){
		for(int j=1;j<=lattice->nuy;j++){
			lattice->u[i][j][lattice->nuz+1]=lattice->u[i][j][1];
		}
	}
	for(int i=1;i<=lattice->nvx;i++){
		for(int j=1;j<=lattice->nvy;j++){
			lattice->v[i][j][lattice->nvz+1]=lattice->v[i][j][1];
		}
	}
	for(int i=1;i<=lattice->nwx;i++){
		for(int j=1;j<=lattice->nwy;j++){
			lattice->w[i][j][lattice->nwz+1]=lattice->w[i][j][1];
		}
	}
}
