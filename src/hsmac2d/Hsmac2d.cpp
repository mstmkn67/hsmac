#include "Hsmac2d.h"

Hsmac2d::Hsmac2d(UDFManager* _in,UDFManager* _out)
:in(_in),out(_out){
}

Hsmac2d::~Hsmac2d(){
	delete bc;
	delete lattice;
}

void Hsmac2d::update(){
	input();
	initial();
	int total=in->i("input.time_evolution.simulation_steps");
	int report=in->i("input.time_evolution.report_steps");
	for(iteration=1;iteration<=total;iteration++){
		advance();
		calc_velocity();
		calc_pressure();
		if(iteration%report==0){
			cout << "report step at " << iteration << endl;
			output();
		}
	}
}

void Hsmac2d::input(){
	dt=in->d("input.time_evolution.dt");
	viscosity=in->d("input.viscosity");
	density=in->d("input.density");
	pressure_iteration=in->i("input.algorithm.pressure_iteration");
	convergence=in->d("input.algorithm.convergence");
	relaxation_coef=in->d("input.algorithm.relaxation_coefficient");
	//
	int nx=in->i("input.size.lattice.x");
	int ny=in->i("input.size.lattice.y");
	double lx=in->d("input.size.length.x");
	double ly=in->d("input.size.length.y");
	string type=in->s("input.flow_type.type");
	if(type=="wall"){
		lattice=new StaggeredLattice2d(nx,ny,lx,ly);
		double u0=in->d("input.flow_type.wall.down.u");
		double v0=in->d("input.flow_type.wall.down.v");
		double u1=in->d("input.flow_type.wall.right.u");
		double v1=in->d("input.flow_type.wall.right.v");
		double u2=in->d("input.flow_type.wall.up.u");
		double v2=in->d("input.flow_type.wall.up.v");
		double u3=in->d("input.flow_type.wall.left.u");
		double v3=in->d("input.flow_type.wall.left.v");
		bc=new Wall(lattice,u0,v0,u1,v1,u2,v2,u3,v3);
	}else if(type=="channel"){
		lattice=new StaggeredLattice2d(nx,ny,lx,ly,nx,ny);
		double u0=in->d("input.flow_type.channel.down.u");
		double v0=in->d("input.flow_type.channel.down.v");
		double u2=in->d("input.flow_type.channel.up.u");
		double v2=in->d("input.flow_type.channel.up.v");
		double p=in->d("input.flow_type.channel.pressure_difference");
		bc=new Channel(lattice,u0,v0,u2,v2,p);
	}else if(type=="LeesEdwards"){
		lattice=new StaggeredLattice2d(nx,ny,lx,ly,nx,ny,nx,ny);
		double gamma=in->d("input.flow_type.LeesEdwards.shear_rate");
		bc=new LeesEdwards(lattice,gamma,dt,&iteration);
	}else if(type=="oscillation"){
		lattice=new StaggeredLattice2d(nx,ny,lx,ly,nx,ny);
		double u2=in->d("input.flow_type.oscillation.u");
		double omega2=in->d("input.flow_type.oscillation.omega");
		bc=new Oscillation(lattice,u2,omega2,dt,&iteration);
	}else{
		cout << "input.flow_type.type is not good" << endl;
		exit(1);
	}
}

void Hsmac2d::initial(){
	for(int i=0;i<=lattice->nux+1;i++){
		for(int j=0;j<=lattice->nuy+1;j++){
			lattice->u[i][j]=0.0;
		}
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		for(int j=0;j<=lattice->nvy+1;j++){
			lattice->v[i][j]=0.0;
		}
	}
	for(int i=0;i<=lattice->nx+1;i++){
		for(int j=0;j<=lattice->ny+1;j++){
			lattice->p[i][j]=0.0;
		}
	}
}

void Hsmac2d::advance(){
	for(int i=0;i<=lattice->nux+1;i++){
		for(int j=0;j<=lattice->nuy+1;j++){
			lattice->u0[i][j]=lattice->u[i][j];
		}
	}
	for(int i=0;i<=lattice->nvx+1;i++){
		for(int j=0;j<=lattice->nvy+1;j++){
			lattice->v0[i][j]=lattice->v[i][j];
		}
	}
}

void Hsmac2d::calc_velocity(){
	//update of u[i][j]
	for(int i=1;i<=lattice->nux;i++){
		for(int j=1;j<=lattice->nuy;j++){
			double vv=0.25*(lattice->v0[i][j]  +lattice->v0[i+1][j]
			               +lattice->v0[i][j-1]+lattice->v0[i+1][j-1]);
			//convection term
			double cnu,cnv;
			if(lattice->u0[i][j]>0.0){
				cnu=lattice->u0[i][j]*(lattice->u0[i][j]-lattice->u0[i-1][j])/lattice->dx;
			}else{
				cnu=lattice->u0[i][j]*(lattice->u0[i+1][j]-lattice->u0[i][j])/lattice->dx;
			}
			if(vv>0.0){
				cnv=vv*(lattice->u0[i][j]-lattice->u0[i][j-1])/lattice->dy;
			}else{
				cnv=vv*(lattice->u0[i][j+1]-lattice->u0[i][j])/lattice->dy;
			}
			//difusion term
			double dif=viscosity*(
				          (lattice->u0[i-1][j]-2.*lattice->u0[i][j]+lattice->u0[i+1][j])/lattice->dx/lattice->dx
			           +(lattice->u0[i][j-1]-2.*lattice->u0[i][j]+lattice->u0[i][j+1])/lattice->dy/lattice->dy);
			//transient velocity 
			lattice->u[i][j]=lattice->u0[i][j]
			                 +dt*(-cnu-cnv+(dif-(lattice->p[i+1][j]-lattice->p[i][j])/lattice->dx)/density );
		}
	}
	//update of v[i][j]
	for(int i=1;i<=lattice->nvx;i++){
		for(int j=1;j<=lattice->nvy;j++){
			double uu=0.25*(lattice->v0[i-1][j+1]+lattice->v0[i][j+1]
			               +lattice->v0[i-1][j]  +lattice->v0[i][j]);
			//convection term
			double cnu,cnv;
			if(uu>0.0){
				cnu=uu*(lattice->v0[i][j]-lattice->v0[i-1][j])/lattice->dx;
			}else{
				cnu=uu*(lattice->v0[i+1][j]-lattice->v0[i][j])/lattice->dx;
			}
			if(lattice->v0[i][j]>0.0){
				cnu=lattice->v0[i][j]*(lattice->v0[i][j]-lattice->v0[i][j-1])/lattice->dy;
			}else{
				cnu=lattice->v0[i][j]*(lattice->v0[i][j+1]-lattice->v0[i][j])/lattice->dy;
			}
			//difusion term
			double dif=viscosity*(
				          (lattice->v0[i-1][j]-2.*lattice->v0[i][j]+lattice->v0[i+1][j])/lattice->dx/lattice->dx
			           +(lattice->v0[i][j-1]-2.*lattice->v0[i][j]+lattice->v0[i][j+1])/lattice->dy/lattice->dy);
			//transient velocity 
			lattice->v[i][j]=lattice->v0[i][j]
			                 +dt*(-cnu-cnv+(dif-(lattice->p[i][j+1]-lattice->p[i][j])/lattice->dy)/density );
		}
	}
	bc->update();
}

void Hsmac2d::calc_pressure(){
	for(int t=1;t<=pressure_iteration;t++){
		int imax=0,jmax=0;
		double dmax=1.e-10;
		double del=2.0*dt/density*(1./lattice->dx/lattice->dx+1./lattice->dy/lattice->dy);
		for(int i=1;i<=lattice->nx;i++){
			for(int j=1;j<=lattice->ny;j++){
				double div=(lattice->u[i][j]-lattice->u[i-1][j])/lattice->dx
				          +(lattice->v[i][j]-lattice->v[i][j-1])/lattice->dy;
				if(abs(div)>abs(dmax)){
					imax=i;jmax=j;dmax=div;
				}
				double delp=-relaxation_coef*div/del;
				lattice->p[i][j]+=delp;
				lattice->u[i][j]+=1.0/density*dt/lattice->dx*delp;
				lattice->u[i-1][j]-=1.0/density*dt/lattice->dx*delp;
				lattice->v[i][j]+=1.0/density*dt/lattice->dy*delp;
				lattice->v[i][j-1]-=1.0/density*dt/lattice->dy*delp;
			}
		}
		bc->update();
		cout << "iteration: " << t << " max_divergence: " << dmax << endl;
		if(abs(dmax)<convergence){
			return;//converge
		}
	}
	cout << "pressure not converge" << endl;
	exit(1);
}

void Hsmac2d::output(){
	out->newRecord();
	for(int i=1;i<=lattice->nx;i++){
		for(int j=1;j<=lattice->ny;j++){
			double u=0.5*(lattice->u[i][j]+lattice->u[i-1][j]);
			out->put("output.lattice[][].u",u,INDEX(i,j));
			double v=0.5*(lattice->v[i][j]+lattice->v[i][j-1]);
			out->put("output.lattice[][].v",v,INDEX(i,j));
			out->put("output.lattice[][].p",lattice->p[i][j],INDEX(i,j));
		}
	}
	out->write();
}
