#include "udf/gourmain.h"
#include "Hsmac3dSimulator.h"

void udfHeaderCheck()
{
	string version("1.0"),engine("hsmac3d");
	cout << "**************************************************************" << endl;
	cout <<  "              " <<  engine << "  " << version << "           " << endl;
	cout << "                                        Masato MAKINO         " << endl;
	cout << "**************************************************************" << endl;
	cout <<  endl;
}

void error_massage(){
	cout << "usage: hsmac3d -I input_udf [-O output_udf] " << endl;
}


int gourmain(UDFManager* in,UDFManager* out){
	udfHeaderCheck();
	Hsmac3dSimulator* sim=new Hsmac3dSimulator(in,out);
	sim->update();
	delete sim;
	return 0;
}
