#include "Hsmac3d.h"

Hsmac3d::Hsmac3d(StaggeredLattice3d* l,Hsmac3dBC* b,
                 double v,double d,double t,
                 int pi,double c,double r)
:lattice(l),bc(b),viscosity(v),density(d),dt(t),pressure_iteration(pi),
convergence(c),relaxation_coef(r){
}

Hsmac3d::~Hsmac3d(){
}

void Hsmac3d::update(){
	advance();
	calc_velocity();
	calc_pressure();
}

void Hsmac3d::initial(){
	for(int i=0;i<=lattice->nux+1;i++){
		for(int j=0;j<=lattice->nuy+1;j++){
			for(int k=0;k<=lattice->nuz+1;k++){
				lattice->u[i][j][k]=0.0;
				lattice->fx[i][j][k]=0.0;
			}
		}
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		for(int j=0;j<=lattice->nvy+1;j++){
			for(int k=0;k<=lattice->nvz+1;k++){
				lattice->v[i][j][k]=0.0;
				lattice->fy[i][j][k]=0.0;
			}
		}
	}
	for(int i=0;i<=lattice->nwx+1;i++){
		for(int j=0;j<=lattice->nwy+1;j++){
			for(int k=0;k<=lattice->nwz+1;k++){
				lattice->w[i][j][k]=0.0;
				lattice->fz[i][j][k]=0.0;
			}
		}
	}
	for(int i=0;i<=lattice->nx+1;i++){
		for(int j=0;j<=lattice->ny+1;j++){
			for(int k=0;k<=lattice->nz+1;k++){
				lattice->p[i][j][k]=0.0;
			}
		}
	}
}

void Hsmac3d::advance(){
	for(int i=0;i<=lattice->nux+1;i++){
		for(int j=0;j<=lattice->nuy+1;j++){
			for(int k=0;k<=lattice->nuz+1;k++){
				lattice->u0[i][j][k]=lattice->u[i][j][k];
			}
		}
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		for(int j=0;j<=lattice->nvy+1;j++){
			for(int k=0;k<=lattice->nvz+1;k++){
				lattice->v0[i][j][k]=lattice->v[i][j][k];
			}
		}
	}
	for(int i=0;i<=lattice->nwx+1;i++){
		for(int j=0;j<=lattice->nwy+1;j++){
			for(int k=0;k<=lattice->nwz+1;k++){
				lattice->w0[i][j][k]=lattice->w[i][j][k];
			}
		}
	}
}

void Hsmac3d::calc_velocity(){
	//update of u[i][j][k]
	for(int i=1;i<=lattice->nux;i++){
		for(int j=1;j<=lattice->nuy;j++){
			for(int k=1;k<=lattice->nuz;k++){
				double vv=0.25*(lattice->v0[i][j][k]  +lattice->v0[i+1][j][k]
				               +lattice->v0[i][j-1][k]+lattice->v0[i+1][j-1][k]);
				double ww=0.25*(lattice->w0[i][j][k]  +lattice->w0[i+1][j][k]
			                 +lattice->w0[i][j][k-1]+lattice->w0[i+1][j][k-1]);
				//convection term
				double cnu,cnv,cnw;
				if(lattice->u0[i][j][k]>0.0){
					cnu=lattice->u0[i][j][k]*(lattice->u0[i][j][k]-lattice->u0[i-1][j][k])/lattice->dx;
				}else{
					cnu=lattice->u0[i][j][k]*(lattice->u0[i+1][j][k]-lattice->u0[i][j][k])/lattice->dx;
				}
				if(vv>0.0){
					cnv=vv*(lattice->u0[i][j][k]-lattice->u0[i][j-1][k])/lattice->dy;
				}else{
					cnv=vv*(lattice->u0[i][j+1][k]-lattice->u0[i][j][k])/lattice->dy;
				}
				if(ww>0.0){
					cnw=ww*(lattice->u0[i][j][k]-lattice->u0[i][j][k-1])/lattice->dz;
				}else{
					cnw=ww*(lattice->u0[i][j][k+1]-lattice->u0[i][j][k])/lattice->dz;
				}
				//difusion term
				double dif=viscosity*(
				          (lattice->u0[i-1][j][k]-2.*lattice->u0[i][j][k]+lattice->u0[i+1][j][k])/lattice->dx/lattice->dx
			           +(lattice->u0[i][j-1][k]-2.*lattice->u0[i][j][k]+lattice->u0[i][j+1][k])/lattice->dy/lattice->dy
			           +(lattice->u0[i][j][k-1]-2.*lattice->u0[i][j][k]+lattice->u0[i][j][k+1])/lattice->dz/lattice->dz);
				//transient velocity 
				lattice->u[i][j][k]=lattice->u0[i][j][k]
			                 +dt*(-cnu-cnv-cnw+(dif-(lattice->p[i+1][j][k]-lattice->p[i][j][k])/lattice->dx)/density+lattice->fx[i][j][k]/density);
			}
		}
	}
	//update of v[i][j][k]
	for(int i=1;i<=lattice->nvx;i++){
		for(int j=1;j<=lattice->nvy;j++){
			for(int k=1;k<=lattice->nvz;k++){
				double uu=0.125*(lattice->u0[i-1][j+1][k]+lattice->u0[i][j+1][k]
				                +lattice->u0[i-1][j][k]  +lattice->u0[i][j][k]);
				double ww=0.125*(lattice->w0[i][j+1][k]  +lattice->w0[i][j][k]
				                +lattice->w0[i][j+1][k-1]+lattice->w0[i][j][k-1]);
				//convection term
				double cnu,cnv,cnw;
				if(uu>0.0){
					cnu=uu*(lattice->v0[i][j][k]-lattice->v0[i-1][j][k])/lattice->dx;
				}else{
					cnu=uu*(lattice->v0[i+1][j][k]-lattice->v0[i][j][k])/lattice->dx;
				}
				if(lattice->v0[i][j][k]>0.0){
					cnu=lattice->v0[i][j][k]*(lattice->v0[i][j][k]-lattice->v0[i][j-1][k])/lattice->dy;
				}else{
					cnu=lattice->v0[i][j][k]*(lattice->v0[i][j+1][k]-lattice->v0[i][j][k])/lattice->dy;
				}
				if(ww>0.0){
					cnw=ww*(lattice->v0[i][j][k]-lattice->v0[i][j][k-1])/lattice->dz;
				}else{
					cnw=ww*(lattice->v0[i][j][k+1]-lattice->v0[i][j][k])/lattice->dz;
				}
				//difusion term
				double dif=viscosity*(
					          (lattice->v0[i-1][j][k]-2.*lattice->v0[i][j][k]+lattice->v0[i+1][j][k])/lattice->dx/lattice->dx
				           +(lattice->v0[i][j-1][k]-2.*lattice->v0[i][j][k]+lattice->v0[i][j+1][k])/lattice->dy/lattice->dy
				           +(lattice->v0[i][j][k-1]-2.*lattice->v0[i][j][k]+lattice->v0[i][j][k+1])/lattice->dz/lattice->dz);
				//transient velocity 
				lattice->v[i][j][k]=lattice->v0[i][j][k]
				                 +dt*(-cnu-cnv-cnw+(dif-(lattice->p[i][j+1][k]-lattice->p[i][j][k])/lattice->dy)/density+lattice->fy[i][j][k]/density);
			}
		}
	}
	//update of w[i][j][k]
	for(int i=1;i<=lattice->nwx;i++){
		for(int j=1;j<=lattice->nwy;j++){
			for(int k=1;k<=lattice->nwz;k++){
				double uu=0.125*(lattice->u0[i-1][j][k+1]+lattice->u0[i][j][k+1]
				                +lattice->u0[i-1][j][k]  +lattice->u0[i][j][k]);
				double vv=0.125*(lattice->v0[i][j-1][k]  +lattice->v0[i][j][k]
				                +lattice->v0[i][j-1][k+1]+lattice->v0[i][j][k+1]);
				//convection term
				double cnu,cnv,cnw;
				if(uu>0.0){
					cnu=uu*(lattice->w0[i][j][k]-lattice->w0[i-1][j][k])/lattice->dx;
				}else{
					cnu=uu*(lattice->w0[i+1][j][k]-lattice->w0[i][j][k])/lattice->dx;
				}
				if(vv>0.0){
					cnu=vv*(lattice->w0[i][j][k]-lattice->w0[i][j-1][k])/lattice->dy;
				}else{
					cnu=vv*(lattice->w0[i][j+1][k]-lattice->w0[i][j][k])/lattice->dy;
				}
				if(lattice->w0[i][j][k]>0.0){
					cnw=lattice->w0[i][j][k]*(lattice->w0[i][j][k]-lattice->w0[i][j][k-1])/lattice->dz;
				}else{
					cnw=lattice->w0[i][j][k]*(lattice->w0[i][j][k+1]-lattice->w0[i][j][k])/lattice->dz;
				}
				//difusion term
				double dif=viscosity*(
					          (lattice->w0[i-1][j][k]-2.*lattice->w0[i][j][k]+lattice->w0[i+1][j][k])/lattice->dx/lattice->dx
				           +(lattice->w0[i][j-1][k]-2.*lattice->w0[i][j][k]+lattice->w0[i][j+1][k])/lattice->dy/lattice->dy
				           +(lattice->w0[i][j][k-1]-2.*lattice->w0[i][j][k]+lattice->w0[i][j][k+1])/lattice->dz/lattice->dz);
				//transient velocity 
					lattice->w[i][j][k]=lattice->w0[i][j][k]
				                 +dt*(-cnu-cnv-cnw+(dif-(lattice->p[i][j][k+1]-lattice->p[i][j][k])/lattice->dz)/density+lattice->fz[i][j][k]/density);
			}
		}
	}
	bc->update();
}

void Hsmac3d::calc_pressure(){
	for(int iteration=1;iteration<=pressure_iteration;iteration++){
		int imax=0,jmax=0;
		double dmax=1.e-10;
		double del=2.0*dt/density*(1./lattice->dx/lattice->dx+1./lattice->dy/lattice->dy+1./lattice->dz/lattice->dz);
		for(int i=1;i<=lattice->nx;i++){
			for(int j=1;j<=lattice->ny;j++){
				for(int k=1;k<=lattice->nz;k++){
					double div=(lattice->u[i][j][k]-lattice->u[i-1][j][k])/lattice->dx
				            +(lattice->v[i][j][k]-lattice->v[i][j-1][k])/lattice->dy
				            +(lattice->w[i][j][k]-lattice->w[i][j][k-1])/lattice->dz;
					if(abs(div)>abs(dmax)){
						imax=i;jmax=j;dmax=div;
					}
					double delp=-relaxation_coef*div/del;
					lattice->p[i][j][k]+=delp;
					lattice->u[i][j][k]+=1.0/density*dt/lattice->dx*delp;
					lattice->u[i-1][j][k]-=1.0/density*dt/lattice->dx*delp;
					lattice->v[i][j][k]+=1.0/density*dt/lattice->dy*delp;
					lattice->v[i][j-1][k]-=1.0/density*dt/lattice->dy*delp;
					lattice->w[i][j][k]+=1.0/density*dt/lattice->dz*delp;
					lattice->w[i][j][k-1]-=1.0/density*dt/lattice->dz*delp;
				}
			}
		}
		bc->update();
		if(iteration%10==0)cout << "iteration: " << iteration << " max_divergence: " << dmax << endl;
		if(abs(dmax)<convergence){
			return;//converge
		}
	}
	cout << "pressure not converge" << endl;
	exit(1);
}

