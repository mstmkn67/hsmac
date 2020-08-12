#include "BoundaryCondition.h"

BoundaryCondition::BoundaryCondition(StaggeredLattice2d* l)
:lattice(l){
}

BoundaryCondition::~BoundaryCondition(){
}

/////////////////////////////////////////////

Wall::Wall(StaggeredLattice2d* lattice,
           double _u0,double _v0,double _u1,double _v1,
           double _u2,double _v2,double _u3,double _v3)
:BoundaryCondition(lattice),u0(_u0),v0(_v0),u1(_u1),v1(_v1),
                            u2(_u2),v2(_v2),u3(_u3),v3(_v3){
}

Wall::~Wall(){
}

void Wall::update(){
	//down
	for(int i=0;i<=lattice->nux+1;i++){
		lattice->u[i][0]=2.*u0-lattice->u[i][1];
	}
	for(int i=0;i<lattice->nvx+1;i++){
		lattice->v[i][0]=v0;
	}
	//right
	for(int j=1;j<=lattice->nuy;j++){
		lattice->u[lattice->nux+1][j]=u1;
	}
	for(int j=1;j<=lattice->nvy;j++){
		lattice->v[lattice->nvx+1][j]=2.*v1-lattice->v[lattice->nvx][j];
	}
	//up
	for(int i=0;i<=lattice->nux+1;i++){
		lattice->u[i][lattice->nuy+1]=2.*u2-lattice->u[i][lattice->nuy];
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		lattice->v[i][lattice->nvy+1]=v2;
	}
	//left
	for(int j=1;j<=lattice->nuy;j++){
		lattice->u[0][j]=u3;
	}
	for(int j=1;j<=lattice->nvy;j++){
		lattice->v[0][j]=2.*v3-lattice->v[1][j];
	}
}
////////////////////////////////////////////////////
Channel::Channel(StaggeredLattice2d* lattice,
                 double _u0,double _v0,double _u2,double _v2,double _delp)
:BoundaryCondition(lattice),u0(_u0),v0(_v0),u2(_u2),v2(_v2),delp(_delp){
}
Channel::~Channel(){
}

void Channel::update(){
	//down
	for(int i=0;i<=lattice->nux+1;i++){
		lattice->u[i][0]=2.*u0-lattice->u[i][1];
	}
	for(int i=0;i<lattice->nvx+1;i++){
		lattice->v[i][0]=v0;
	}
	//right
	for(int j=1;j<=lattice->ny;j++){
		lattice->p[lattice->nx+1][j]=lattice->p[1][j]-delp;
	}
	for(int j=1;j<=lattice->nuy;j++){
		lattice->u[lattice->nux+1][j]=lattice->u[1][j];
	}
	for(int j=1;j<=lattice->nvy;j++){
		lattice->v[lattice->nvx+1][j]=lattice->v[1][j];
	}
	//up
	for(int i=0;i<=lattice->nux+1;i++){
		lattice->u[i][lattice->nuy+1]=2.*u2-lattice->u[i][lattice->nuy];
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		lattice->v[i][lattice->nvy+1]=v2;
	}
	//left
	for(int j=1;j<=lattice->nuy;j++){
		lattice->u[0][j]=lattice->u[lattice->nux][j];
	}
	for(int j=1;j<=lattice->nvy;j++){
		lattice->v[0][j]=lattice->v[lattice->nvx][j];
	}
}
////////////////////////////////////////////////////
Oscillation::Oscillation(StaggeredLattice2d* lattice,
                         double _u2,double _o2,double _dt,int* _ite)
:BoundaryCondition(lattice),u2(_u2),omega2(_o2),dt(_dt),iteration(_ite){
}
Oscillation::~Oscillation(){
}

void Oscillation::update(){
	//down
	for(int i=0;i<=lattice->nux+1;i++){
		lattice->u[i][0]=-lattice->u[i][1];
	}
	for(int i=0;i<lattice->nvx+1;i++){
		lattice->v[i][0]=0.0;
	}
	//right
	for(int j=1;j<=lattice->ny;j++){
		lattice->p[lattice->nx+1][j]=lattice->p[1][j];
	}
	for(int j=1;j<=lattice->nuy;j++){
		lattice->u[lattice->nux+1][j]=lattice->u[1][j];
	}
	for(int j=1;j<=lattice->nvy;j++){
		lattice->v[lattice->nvx+1][j]=lattice->v[1][j];
	}
	//up
	for(int i=0;i<=lattice->nux+1;i++){
		double a=u2*sin(omega2*(*iteration)*dt);
		lattice->u[i][lattice->nuy+1]=2.*a-lattice->u[i][lattice->nuy];
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		lattice->v[i][lattice->nvy+1]=0.0;
	}
	//left
	for(int j=1;j<=lattice->nuy;j++){
		lattice->u[0][j]=lattice->u[lattice->nux][j];
	}
	for(int j=1;j<=lattice->nvy;j++){
		lattice->v[0][j]=lattice->v[lattice->nvx][j];
	}
}
////////////////////////////////////////////////////
LeesEdwards::LeesEdwards(StaggeredLattice2d* lattice,double s,double d,int* ite)
:BoundaryCondition(lattice),shear_rate(s),dt(d),iteration(ite){
}

LeesEdwards::~LeesEdwards(){
}

void LeesEdwards::update(){
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
		lattice->u[i][0]=lattice->u[i][lattice->nuy]-shear_rate*lattice->ly;
	}
	for(int i=0;i<lattice->nvx+1;i++){
		lattice->v[i][0]=lattice->v[i][lattice->nvy];
	}
	//right
	for(int j=1;j<=lattice->ny;j++){
		lattice->p[lattice->nx+1][j]=lattice->p[1][j];
	}
	for(int j=1;j<=lattice->nuy;j++){
		lattice->u[lattice->nux+1][j]=lattice->u[1][j];
	}
	for(int j=1;j<=lattice->nvy;j++){
		lattice->v[lattice->nvx+1][j]=lattice->v[1][j];
	}
	//up
	for(int i=0;i<=lattice->nx+1;i++){
		int i0=i-Dn-1;int i1=i-Dn;
		if(i0<0)i0+=lattice->nx;if(i1<0)i1+=lattice->nx;
		lattice->p[i][lattice->ny+1]=(1.-a)*lattice->p[i0][1]+a*lattice->p[i1][1];
	}
	for(int i=0;i<=lattice->nux+1;i++){
		int i0=i-Dn-1;int i1=i-Dn;
		if(i0<0)i0+=lattice->nux;if(i1<0)i1+=lattice->nux;
		lattice->u[i][lattice->nuy+1]=(1.-a)*lattice->u[i0][1]+a*lattice->u[i1][1]+shear_rate*lattice->ly;
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		int i0=i-Dn-1;int i1=i-Dn;
		if(i0<0)i0+=lattice->nvx;if(i1<0)i1+=lattice->nvx;
		lattice->v[i][lattice->nvy+1]=(1.-a)*lattice->v[i0][1]+a*lattice->v[i1][1];
	}
	//left
	for(int j=1;j<=lattice->nuy;j++){
		lattice->u[0][j]=lattice->u[lattice->nux][j];
	}
	for(int j=1;j<=lattice->nvy;j++){
		lattice->v[0][j]=lattice->v[lattice->nvx][j];
	}
}
