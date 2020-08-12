#ifndef _HSMAC3D_SIMULATOR_H_
#define _HSMAC3D_SIMULATOR_H_

#include "udfmanager.h"
#include "flow/Hsmac3d.h"
#include <cmath>
using namespace std;

class Hsmac3dSimulator{
public:
	Hsmac3dSimulator(UDFManager* in,UDFManager* out);
	virtual ~Hsmac3dSimulator();
	virtual void update();
protected:
	virtual void input_flow();
	virtual void output();
private:
	UDFManager *in,*out;
	int iteration;
	//
	Hsmac3d *flow;
	StaggeredLattice3d *lattice;
	Hsmac3dBC* fbc;
};

#endif // _HSMAC3D_SIMULATOR_H_
