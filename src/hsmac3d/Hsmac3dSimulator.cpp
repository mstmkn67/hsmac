#include "Hsmac3dSimulator.h"

Hsmac3dSimulator::Hsmac3dSimulator(UDFManager* _in,UDFManager* _out)
:in(_in),out(_out){
}

Hsmac3dSimulator::~Hsmac3dSimulator(){
	delete flow;
	delete fbc;
	delete lattice;
}

void Hsmac3dSimulator::update(){
	input_flow();
	flow->initial();
	int total=in->i("input.time_evolution.simulation_steps");
	int report=in->i("input.time_evolution.report_steps");
	for(iteration=1;iteration<=total;iteration++){
		flow->update();
		if(iteration%report==0){
			cout << "report step at " << iteration << endl;
			output();
		}
	}
}

void Hsmac3dSimulator::input_flow(){
	double dt=in->d("input.time_evolution.dt");
	double viscosity=in->d("input.viscosity");
	double density=in->d("input.density");
	int pressure_iteration=in->i("input.flow.algorithm.pressure_iteration");
	double convergence=in->d("input.flow.algorithm.convergence");
	double relaxation_coef=in->d("input.flow.algorithm.relaxation_coefficient");
	//
	int nx=in->i("input.size.lattice.x");
	int ny=in->i("input.size.lattice.y");
	int nz=in->i("input.size.lattice.z");
	double lx=in->d("input.size.length.x");
	double ly=in->d("input.size.length.y");
	double lz=in->d("input.size.length.z");
	string type=in->s("input.flow.type");
	if(type=="wall"){
		lattice=new StaggeredLattice3d(nx,ny,nz,lx,ly,lz);
		double u0=in->d("input.flow.wall.down.u");
		double v0=in->d("input.flow.wall.down.v");
		double w0=in->d("input.flow.wall.down.w");
		double u1=in->d("input.flow.wall.right.u");
		double v1=in->d("input.flow.wall.right.v");
		double w1=in->d("input.flow.wall.right.w");
		double u2=in->d("input.flow.wall.up.u");
		double v2=in->d("input.flow.wall.up.v");
		double w2=in->d("input.flow.wall.up.w");
		double u3=in->d("input.flow.wall.left.u");
		double v3=in->d("input.flow.wall.left.v");
		double w3=in->d("input.flow.wall.left.w");
		double u4=in->d("input.flow.wall.front.u");
		double v4=in->d("input.flow.wall.front.v");
		double w4=in->d("input.flow.wall.front.w");
		double u5=in->d("input.flow.wall.back.u");
		double v5=in->d("input.flow.wall.back.v");
		double w5=in->d("input.flow.wall.back.w");
		fbc=new Hsmac3dWall(lattice,u0,v0,w0,u1,v1,w1,u2,v2,w2,u3,v3,w3,u4,v4,w4,u5,v5,w5);
	}else if(type=="channel"){
		lattice=new StaggeredLattice3d(nx,ny,nz,lx,ly,lz,nx,ny,nz);
		double u0=in->d("input.flow.channel.down.u");
		double v0=in->d("input.flow.channel.down.v");
		double w0=in->d("input.flow.channel.down.w");
		double u2=in->d("input.flow.channel.up.u");
		double v2=in->d("input.flow.channel.up.v");
		double w2=in->d("input.flow.channel.up.w");
		double u4=in->d("input.flow.channel.front.u");
		double v4=in->d("input.flow.channel.front.v");
		double w4=in->d("input.flow.channel.front.w");
		double u5=in->d("input.flow.channel.back.u");
		double v5=in->d("input.flow.channel.back.v");
		double w5=in->d("input.flow.channel.back.w");
		double p=in->d("input.flow.channel.pressure_difference");
		fbc=new Hsmac3dChannel(lattice,u0,v0,w0,u2,v2,w2,u4,v4,w4,u5,v5,w5,p);
	}else if(type=="LeesEdwards"){
		lattice=new StaggeredLattice3d(nx,ny,nz,lx,ly,lz,nx,ny,nz,nx,ny,nz,nx,ny,nz);
		double gamma=in->d("input.flow.LeesEdwards.shear_rate");
		fbc=new Hsmac3dLeesEdwards(lattice,gamma,dt,&iteration);
	}else{
		cout << "input.flow_type.type is not good" << endl;
		exit(1);
	}
	flow=new Hsmac3d(lattice,fbc,viscosity,density,dt,
	                 pressure_iteration,convergence,relaxation_coef);
}

void Hsmac3dSimulator::output(){
	out->newRecord();
	for(int i=1;i<=lattice->nx;i++){
		for(int j=1;j<=lattice->ny;j++){
			for(int k=1;k<=lattice->nz;k++){
				double u=0.5*(lattice->u[i][j][k]+lattice->u[i-1][j][k]);
				out->put("output.lattice[][][].u",u,INDEX(i,j,k));
				double v=0.5*(lattice->v[i][j][k]+lattice->v[i][j-1][k]);
				out->put("output.lattice[][][].v",v,INDEX(i,j,k));
				double w=0.5*(lattice->w[i][j][k]+lattice->w[i][j][k-1]);
				out->put("output.lattice[][][].w",w,INDEX(i,j,k));
				out->put("output.lattice[][][].p",lattice->p[i][j][k],INDEX(i,j,k));
			}
		}
	}
	out->write();
}
