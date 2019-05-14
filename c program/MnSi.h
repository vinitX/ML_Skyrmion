#ifndef _MnSi_h
#define _MnSi_h

#include <cmath>
#include <ctime>
#include <iostream>
#include "ran3.h"
using namespace std;

const double PI = 4 * atan(1.0);

struct Spin;
struct Spin
{
	double theta, phi;
	
	Spin *mx, *px, *my, *py, *mz, *pz;
};

class MnSi
{
public:
	MnSi(double J, double K, double A1, double A2, int Lx, int Ly, int Lz);
	~MnSi();

	void run_MC(double t, double H, int mc_step, int trash);

	void allocate_randomly();
	void allocate_from_file_angle(istream &in);
	void allocate_from_file_xyz(istream &in);

	void store_to_file_angle(ostream &out);
	void store_to_file_xyz(ostream &out);

	void store_data(ostream &out, double temperature);

private:
	double cal_energy(Spin *pos, double H);
	void cal_energy_mag_chiral(Spin *pos, double &energy, 
							double &mx, double &my, double &mz, double &chirality, double H);

private:
	int Lx, Ly, Lz, size;
	Spin *m;
	double J, K, A1, A2;
	double energy, magnetization, cv, msus, chirality, chirality2;
};

#define M(x,y,z)	m[Ly * Lz * (x) + Lz * (y) + (z)]


MnSi::MnSi(double J, double K, double A1, double A2, int Lx, int Ly, int Lz)
{
	this->J = J;
	this->K = K;
	this->A1 = A1;
	this->A2 = A2;

	this->Lx = Lx;
	this->Ly = Ly;	
	this->Lz = Lz;	
	size = Lx*Ly*Lz;

	m = new Spin[size];

	for(int x = 0; x != Lx; x++)
	{
		for(int y = 0; y != Ly; y++)
		{
			for(int z = 0; z != Lz; z++)
			{
				if(x != 0)
					M(x,y,z).mx = &M(x-1,y,z);
				else
					M(x,y,z).mx = &M(Lx-1,y,z);

				if(x != Lx - 1)
					M(x,y,z).px = &M(x+1,y,z);
				else
					M(x,y,z).px = &M(0,y,z);

				if(y != 0)
					M(x,y,z).my = &M(x,y-1,z);
				else
					M(x,y,z).my = &M(x,Ly-1,z);

				if(y != Ly - 1)
					M(x,y,z).py = &M(x,y+1,z);
				else
					M(x,y,z).py = &M(x,0,z);

				if(z != 0)
					M(x,y,z).mz = &M(x,y,z-1);
				else
					M(x,y,z).mz = &M(x,y,Lz-1);

				if(z != Lz - 1)
					M(x,y,z).pz = &M(x,y,z+1);
				else
					M(x,y,z).pz = &M(x,y,0);
			}
		}
	}

	allocate_randomly();
	/*
	ifstream fin("SS.txt");
	for(int i = 0; i != size; i++)
	{
		fin >> m[i].phi >> m[i].theta;
	}
	
	fin.close();
    */


	
}

MnSi::~MnSi()
{
	delete[] m;
}

void MnSi::run_MC(double temperature, double H, int mc_step, int trash)
{
	long idum = (long)0 - (unsigned)time(NULL);

	double d_energy, d_sine;
	double d_phi = 2*PI, d_theta = PI;
	int x,y,z;
	Spin new_spin;
	int accepted_num_per_mc;

	double ener, mx, my, mz, chiral;
	double summ=0, summ2=0;
	double sumE=0, sumE2=0;
	double sumchirality = 0, sumchirality2 = 0;

	double temp_ener, temp_mx, temp_my, temp_mz, temp_chiral;

	for(int step = 0; step != mc_step; step++)
	{
		accepted_num_per_mc = 0;
		for(x = 0; x != Lx; x++)   
		{   
			for(y = 0; y != Ly; y++)
			{
				for(z = 0; z != Lz; z++)
				{
					new_spin = M(x,y,z);

					new_spin.phi += d_phi * (ran3(&idum)-0.5);
					new_spin.theta += d_theta * (ran3(&idum)-0.5);

					if(new_spin.phi < 0)			new_spin.phi += 2 * PI;
					else if(new_spin.phi > 2 * PI)	new_spin.phi -= 2 * PI;

					if(new_spin.theta < 0)			new_spin.theta += PI;
					else if(new_spin.theta > PI)	new_spin.theta -= PI;
					d_energy = cal_energy(&new_spin, H) 
						- cal_energy(&M(x,y,z), H);
					d_sine = log(fabs(sin(new_spin.theta))) 
						- log(fabs(sin(M(x,y,z).theta)));

					if(d_energy < 0 
						|| exp(- d_energy/temperature + d_sine) > ran3(&idum))
					{
						M(x,y,z) = new_spin;
						accepted_num_per_mc++;
					}
				}
			}
		}

		if(fabs((double)accepted_num_per_mc) / (size) < .5)
		{
			d_phi *= .9;
			d_theta *= .9;
		}

		if(step > trash - 1)
		{
			ener = mx = my = mz = chiral = 0;
			for(x = 0; x != Lx; x++)
			{
				for(y = 0; y != Ly; y++)   
				{
					for(z = 0; z != Lz; z++)   
					{
						cal_energy_mag_chiral(&M(x,y,z), temp_ener, 
							temp_mx, temp_my, temp_mz, temp_chiral, H);

						ener += temp_ener;

						mx += temp_mx;
						my += temp_my;
						mz += temp_mz;

						chiral += temp_chiral;
					}
				}
			}

			sumchirality += chiral;
			sumchirality2 += chiral*chiral;

			summ += sqrt(mx*mx + my*my + mz*mz);
			summ2 += mx*mx + my*my + mz*mz;

			sumE += ener;
			sumE2 += ener*ener;
		}
	}

	summ /= (mc_step-trash);
	summ2 /= (mc_step-trash);

	sumE /= (mc_step-trash); 
	sumE2 /= (mc_step-trash);

	sumchirality /= (mc_step-trash);
	sumchirality2 /= (mc_step-trash);

	energy = sumE/size;
	magnetization = summ/size;

	cv = (sumE2 - sumE*sumE)/temperature/temperature/size;

	msus= (summ2 - summ*summ)/temperature/size;

	chirality = sumchirality;
	chirality2 = sumchirality2;

}

void MnSi::allocate_randomly()
{
	long idum = (long)0 - (unsigned)time(NULL);

	for(int i = 0; i != size; i++)
	{
		m[i].phi = 2 * PI * ran3(&idum);
		m[i].theta = PI * ran3(&idum);
	}
}

void MnSi::allocate_from_file_angle(istream &in)
{
	for(int i = 0; i != size; i++)
	{
		in >> m[i].phi >> m[i].theta;
	}
}

void MnSi::allocate_from_file_xyz(istream &in) /// 수정해야함.
{
	double x,y,z;
	for(int i = 0; i != size; i++)
	{
			in >> x >> y >> z;
			m[i].phi = fabs(atan(y/x));
			if(x < 0 && y < 0)
				m[i].phi += PI;
			else if(x < 0 && y > 0)
				m[i].phi = PI - m[i].phi;
			else if(x > 0 && y < 0)
				m[i].phi = 2 * PI - m[i].phi;

			m[i].theta = fabs(atan(sqrt(x*x+y*y)/z));
			if(z < 0)
				m[i].theta = PI - m[i].theta;
	}
}

void MnSi::store_to_file_angle(ostream &out)
{
	for(int i = 0; i != size; i++)
	{
		out << m[i].phi << '\t' << m[i].theta << endl;
	}
}

void MnSi::store_to_file_xyz(ostream &out)
{
	for(int i = 0; i != size; i++)
	{
		out << sin(m[i].theta)*cos(m[i].phi) << '\t' 
			<< sin(m[i].theta)*sin(m[i].phi) << '\t' 
			<< cos(m[i].theta) << endl;
	}
}

void MnSi::store_data(ostream &out, double temperature)
{
	out << temperature << '\t'
		<< energy << '\t'
		<< cv << '\t'
		<< magnetization << '\t'
		<< msus << '\t'
		<< chirality << '\t'
		<< chirality2 << endl;
}

double MnSi::cal_energy(Spin *pos, double H)
{
	double m_x, m_y, m_z;
	double m_px_x, m_px_y, m_px_z;
	double m_mx_x, m_mx_y, m_mx_z;
	double m_py_x, m_py_y, m_py_z;
	double m_my_x, m_my_y, m_my_z;
	double m_pz_x, m_pz_y, m_pz_z;
	double m_mz_x, m_mz_y, m_mz_z;

	m_x = sin(pos->theta) * cos(pos->phi);
	m_y = sin(pos->theta) * sin(pos->phi);
	m_z = cos(pos->theta);

	m_px_x = sin(pos->px->theta) * cos(pos->px->phi);
	m_px_y = sin(pos->px->theta) * sin(pos->px->phi);
	m_px_z = cos(pos->px->theta);

	m_mx_x = sin(pos->mx->theta) * cos(pos->mx->phi);
	m_mx_y = sin(pos->mx->theta) * sin(pos->mx->phi);
	m_mx_z = cos(pos->mx->theta);

	m_py_x = sin(pos->py->theta) * cos(pos->py->phi);
	m_py_y = sin(pos->py->theta) * sin(pos->py->phi);
	m_py_z = cos(pos->py->theta);

	m_my_x = sin(pos->my->theta) * cos(pos->my->phi);
	m_my_y = sin(pos->my->theta) * sin(pos->my->phi);
	m_my_z = cos(pos->my->theta);

	m_pz_x = sin(pos->pz->theta) * cos(pos->pz->phi);
	m_pz_y = sin(pos->pz->theta) * sin(pos->pz->phi);
	m_pz_z = cos(pos->pz->theta);

	m_mz_x = sin(pos->mz->theta) * cos(pos->mz->phi);
	m_mz_y = sin(pos->mz->theta) * sin(pos->mz->phi);
	m_mz_z = cos(pos->mz->theta);

	double temp = m_x * m_px_x + m_y * m_px_y + m_z * m_px_z;
	temp += m_x * m_mx_x + m_y * m_mx_y + m_z * m_mx_z;
	temp += m_x * m_py_x + m_y * m_py_y + m_z * m_py_z;
	temp += m_x * m_my_x + m_y * m_my_y + m_z * m_my_z;
	if(Lz != 1)
	{
		temp += m_x * m_pz_x + m_y * m_pz_y + m_z * m_pz_z;
		temp += m_x * m_mz_x + m_y * m_mz_y + m_z * m_mz_z;
	}

	double retval = (-J) * temp;

	temp = m_z * m_px_x - m_px_z * m_x;
	temp += m_mx_z * m_x - m_z * m_mx_x;
	temp += m_y * m_py_z - m_py_y * m_z;
	temp += m_my_y * m_z - m_y * m_my_z;  /*    anti-skyrmion    */
	if(Lz != 1)
	{
		temp += m_x * m_pz_y - m_pz_x * m_y;
		temp += - m_x * m_mz_y + m_mz_x * m_y;
	}

	retval += (-K) * temp;

	temp = m_x * m_px_x + m_x * m_mx_x;
	temp += m_y * m_py_y + m_y * m_my_y;
	if(Lz != 1)
	{
		temp += m_z * m_pz_z + m_z * m_mz_z;
	}

	retval += (-A2) * temp;

	temp = m_x*m_x*m_x*m_x + m_y*m_y*m_y*m_y + m_z*m_z*m_z*m_z;

	retval += (A1) * temp;

    retval += (-H) * m_z;

	return retval;
}

void MnSi::cal_energy_mag_chiral(Spin *pos, double &energy, 
							 double &m_x, double &m_y, double &m_z, double &chirality, double H)
{
	double m_px_x, m_px_y, m_px_z;
	double m_mx_x, m_mx_y, m_mx_z;
	double m_py_x, m_py_y, m_py_z;
	double m_my_x, m_my_y, m_my_z;
	double m_pz_x, m_pz_y, m_pz_z;
	double m_mz_x, m_mz_y, m_mz_z;

	m_x = sin(pos->theta) * cos(pos->phi);
	m_y = sin(pos->theta) * sin(pos->phi);
	m_z = cos(pos->theta);

	m_px_x = sin(pos->px->theta) * cos(pos->px->phi);
	m_px_y = sin(pos->px->theta) * sin(pos->px->phi);
	m_px_z = cos(pos->px->theta);

	m_mx_x = sin(pos->mx->theta) * cos(pos->mx->phi);
	m_mx_y = sin(pos->mx->theta) * sin(pos->mx->phi);
	m_mx_z = cos(pos->mx->theta);

	m_py_x = sin(pos->py->theta) * cos(pos->py->phi);
	m_py_y = sin(pos->py->theta) * sin(pos->py->phi);
	m_py_z = cos(pos->py->theta);

	m_my_x = sin(pos->my->theta) * cos(pos->my->phi);
	m_my_y = sin(pos->my->theta) * sin(pos->my->phi);
	m_my_z = cos(pos->my->theta);

	m_pz_x = sin(pos->pz->theta) * cos(pos->pz->phi);
	m_pz_y = sin(pos->pz->theta) * sin(pos->pz->phi);
	m_pz_z = cos(pos->pz->theta);

	m_mz_x = sin(pos->mz->theta) * cos(pos->mz->phi);
	m_mz_y = sin(pos->mz->theta) * sin(pos->mz->phi);
	m_mz_z = cos(pos->mz->theta);

	double temp = m_x * m_px_x + m_y * m_px_y + m_z * m_px_z;
	temp += m_x * m_py_x + m_y * m_py_y + m_z * m_py_z;
	if(Lz != 1)
		temp += m_x * m_pz_x + m_y * m_pz_y + m_z * m_pz_z;

	energy = (-J) * temp;

	temp = m_z * m_px_x - m_px_z * m_x;
	temp += m_y * m_py_z - m_py_y * m_z; /*    anti-skyrmion    */
	if(Lz != 1)
		temp += m_x * m_pz_y - m_pz_x * m_y;

	energy += (-K) * temp;

	temp = m_x * m_px_x + m_y * m_py_y;
	if(Lz != 1)
		temp += m_z * m_pz_z;

	energy += (-A2) * temp;

	temp = m_x*m_x*m_x*m_x + m_y*m_y*m_y*m_y + m_z*m_z*m_z*m_z;

	energy += (A1) * temp;

    energy += (-H) * m_z;


	chirality = (m_px_y * m_py_z - m_px_z * m_py_y) * m_x;
	chirality += (m_px_z * m_py_x - m_px_x * m_py_z) * m_y;
	chirality += (m_px_x * m_py_y - m_px_y * m_py_x) * m_z;
	chirality += (m_mx_y * m_my_z - m_mx_z * m_my_y) * m_x;
	chirality += (m_mx_z * m_my_x - m_mx_x * m_my_z) * m_y;
	chirality += (m_mx_x * m_my_y - m_mx_y * m_my_x) * m_z;
}

#endif
