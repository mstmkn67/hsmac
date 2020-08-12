#include "udf/gourmain.h"
#include "Hsmac2d.h"

void udfHeaderCheck()
{
	string version("1.0"),engine("hsmac2d");
	cout << "**************************************************************" << endl;
	cout <<  "              " <<  engine << "  " << version << "           " << endl;
	cout << "                                        Masato MAKINO         " << endl;
	cout << "**************************************************************" << endl;
	cout <<  endl;
}

void error_massage(){
	cout << "usage: hsmac2d -I input_udf [-O output_udf] " << endl;
}


int gourmain(UDFManager* in,UDFManager* out){
	udfHeaderCheck();
	Hsmac2d* sim=new Hsmac2d(in,out);
	sim->update();
	delete sim;
	return 0;
}
