#include <iostream>
#include <fstream>
#include "MnSi.h"
using namespace std;

int main(int argc, char **argv)
{
	int mc_step = 200000;
	int trash= 100000;
	int Lx = 36, Ly = 36, Lz = 1;

	double period = 6;
	double H, J, K, A1, A2;

    H = 2.0;
	J = 1;
	K = J * sqrt(2.0) * tan(2 * PI / period);
	A1 = 0.0;
	A2 = 0.0;
	
	MnSi model(J, K, A1, A2, Lx, Ly, Lz);

	char file[50];
	ofstream out;
	
	
	//for(double H = 0.0; H < 3.51; H += 0.5)
	//{
	for(double tem = 0.9; tem > 0.09; tem -= 0.1)
	{
		
		model.run_MC(tem, H, mc_step, trash);

		sprintf(file, "conf_H%.2f_T%.2f_dec.txt", H, tem);
		out.open(file);
		model.store_to_file_xyz(out);
		out.close();

		sprintf(file, "thetaphi_H%.2f_T%.2f_dec.txt", H, tem);
		out.open(file);
		model.store_to_file_angle(out);
		out.close();

		out.open("data.txt", ios_base::app);
		model.store_data(out, tem);
		out.close();
	}

	return 0;
}
