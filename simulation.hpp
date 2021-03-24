#ifndef SIMULATION_HPP
#define SIMULATION_HPP
#include "molecule.hpp"

class simulation{
	vector<molecule*> molptrs;
	box bx;
	ofstream xyz;
	flt eps_rlj, eps, sig, w, theta0, dtheta, mg, KbT, U, zmin, r_cut, sigma;
	flt phi_min, phi_max, w4, C0, C1, C2, C3, alpha, sig6; 
	vector<flt> params;
	Vector3d dr, e_strain;
	Quaterniond dq;
	ArrayXXd Uij;
	ArrayXXd BondTime;
	vector<uint> indsUij;
	vector<flt> dUs;
	int ConsecRej;
	uint N;
	bool BondFlag;
	flt BondRate;

	public:
		simulation(vector<molecule*>,flt (*UwallFunc)(vector<flt>,flt), flt, flt, flt, flt, flt, flt, vector<flt>, flt, Vector3d, Quaterniond, Vector3d,string* , Vector3d, uint, flt);
		simulation(vector<molecule*>,flt (*UwallFunc)(vector<flt>,flt), flt, flt, flt, flt, flt, flt, vector<flt>, flt, Vector3d, Quaterniond, Vector3d,string* , Vector3d, flt, int, uint, flt);
		void reset_Rej();
		int get_Rej();
		void equate_molptrs(vector<molecule*>);
		void resetU(uint);
		void set_N(uint);
		void apply_strain();
		void printRsQs2screen();
		void write_xyzfile(ofstream&, flt);
		void write_theta_horz(ofstream&, flt, uint);
		flt get_Ui(uint);
		flt get_Ui_new(uint, flt);
		flt get_Uij(flt, flt, flt, uint, uint);
		flt (*Uwall)(vector<flt>, flt);
		bool metMC(flt dU, flt z);
		void timestep(uint);
		void update_Uij(uint i);
		void timestep_new(uint);
		void update_BondTime(uint);
};



#endif 