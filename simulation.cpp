#include "simulation.hpp"


simulation::simulation(vector<molecule*> molptrs0, flt (*UwallFunc)(vector<flt>,flt), flt eps_rlj0, flt eps0, flt sig0, flt w0, flt theta00, flt dtheta0, vector<flt> params0, flt KbT0, Vector3d dr0, Quaterniond dq0, Vector3d e_strain0, string* bound_arr, Vector3d Ls0, uint Nseg, flt BondRate0){
	// Constructor for simulation

	this->Uwall = UwallFunc;
	equate_molptrs(molptrs0);
	Vector3d Ls;
	for(uint i=0;i<3;i++){
		dr[i] = dr0[i];
		e_strain[i] = e_strain0[i];
		Ls[i] = Ls0[i];
	}
	//molptrs = molptrs;
	bx = box(bound_arr, Ls);
	eps_rlj = eps_rlj0;
	eps = eps0;
	sig = sig0;
	w = w0;
	sigma = sig-w;
	theta0 = theta00;
	dtheta = dtheta0;
	params = params0;
	KbT = KbT0;
	dr = dr0;
	dq = dq0;
	e_strain = e_strain;
	zmin = 0;
	N=1;
    molecule &m = *molptrs[N];
    U = Uwall(params,m.get_zcom());
	//U = mg*m.get_zcom(); // Gravitational potential energy
	ConsecRej = 0;
    r_cut = sig + w;
    dUs.clear();
    indsUij.clear();

    Uij = ArrayXXd(Nseg,Nseg);

    if(eps0>0){
    	BondRate = BondRate0 / eps0;	
    }
    if(BondRate>0){
   		BondTime = ArrayXXd(Nseg,Nseg);
    	BondFlag = true;
    }else{ 
    	BondFlag = false;
    }

	phi_min = theta0-dtheta;
	phi_max = theta0+dtheta;
    sig6 = pow(sig,6.);
	w4 = pow(w,4.);
    C0 = eps/w4;
    C1 = 2.*(3.*sig*sig - w*w);
    C2 = 4.*(sig*sig-w*w)*sig;
    C3 = (sig*sig - w*w)*(sig*sig + w*w);
    alpha = M_PI/dtheta;
}
simulation::simulation(vector<molecule*> molptrs0, flt (*UwallFunc)(vector<flt>, flt), flt eps_rlj0, flt eps0, flt sig0, flt w0, flt theta00, flt dtheta0, vector<flt> params0, flt KbT0, Vector3d dr0, Quaterniond dq0, Vector3d e_strain0, string* bound_arr, Vector3d Ls0, flt zmin0, int N0, uint Nseg, flt BondRate0){
	// Constructor for simulation

	this->Uwall = UwallFunc;
	equate_molptrs(molptrs0);
	Vector3d Ls;
	for(uint i=0;i<3;i++){
		dr[i] = dr0[i];
		e_strain[i] = e_strain0[i];
		Ls[i] = Ls0[i];
	}
	//molptrs = molptrs;
	bx = box(bound_arr, Ls);
	eps_rlj = eps_rlj0;
	eps = eps0;
	sig = sig0;
	w = w0;
	sigma = sig-w;
	theta0 = theta00;
	dtheta = dtheta0;
	params = params0;
	//mg = mg0;
	//cout << "MG: " << mg << endl;
	KbT = KbT0;
	dr = dr0;
	dq = dq0;
	e_strain = e_strain;
	zmin = zmin0;
	N=N0;
    molecule &m = *molptrs[N];
    U = Uwall(params,m.get_zcom());
	//U = mg*m.get_zcom(); // Gravitational potential energy
	ConsecRej = 0;

    r_cut = sig + w;

    Uij = ArrayXXd(Nseg,Nseg);
    if(eps0>0){
    	BondRate = BondRate0 / eps0;	
    }
    if(BondRate>0){
   		BondTime = ArrayXXd(Nseg,Nseg);
    	BondFlag = true;
    }else{ 
    	BondFlag = false;
    }

	phi_min = theta0-dtheta;
	phi_max = theta0+dtheta;
    sig6 = pow(sig,6.);
	w4 = pow(w,4.);
    C0 = eps/w4;
    C1 = 2.*(3.*sig*sig - w*w);
    C2 = 4.*(sig*sig-w*w)*sig;
    C3 = (sig*sig - w*w)*(sig*sig + w*w);
    alpha = M_PI/dtheta;
}

void simulation::reset_Rej(){
	ConsecRej = 0;
	return;
}

int simulation::get_Rej(){
	return ConsecRej;
}

void simulation::equate_molptrs(vector<molecule*> new_molptrs){
	uint Nseg = new_molptrs.size();
    molptrs.reserve(Nseg);
    for(uint i=0;i<Nseg;i++){
    	//cout << "In Equate molptrs: i: " << i << " Address: " << new_molptrs[i] << endl;
    	molptrs[i] = new_molptrs[i];
    	//molecule &m = *molptrs[i];
    	//m.print_rq2screen();
    }
}

void simulation::set_N(uint N0){
	N = N0;
}

void simulation::resetU(uint i){
	// Reset U to the mg of the index of interest
    molecule &m = *molptrs[i];
    U = Uwall(params,m.get_zcom());
	//U = mg*m.get_zcom();
}

void simulation::apply_strain(){
	// Apply strain to box
	bx.apply_strain(e_strain);
	// Apply strain to molecules
	for(uint i=0;i<N;i++){
        molecule &m = *molptrs[i];
		m.strain(e_strain);
	}

}

void simulation::printRsQs2screen(){
	Vector3d Ls = bx.get_Ls();
	cout << "\nBox dimensions: " << Ls[0] << ", " << Ls[1] << ", " << Ls[2] << endl;
	cout << "N: " << N <<endl;
	cout<< "i: x\ty\tz\tq.w\tq.x\tq.y\tq.z\n";

	for(uint i=0; i<N; i++){
		cout<< i <<":\t";
        molecule &m = *molptrs[i];
		m.print_rq2screen();
	}
	cout << "i:\t ux\tuy\tuz\n";
	for(uint i=0; i<N; i++){
		cout<< i <<":\t";
        molecule &m = *molptrs[i];
		m.print_u2screen();
	}
	cout << endl;
	return;
}

void simulation::write_xyzfile(ofstream& xyz, flt time){
	// Get current box dimensions
	Vector3d Ls = bx.get_Ls();
	
	// Print out header for timestep
	xyz << N <<"\nLattice=\""<< Ls[0] << " 0. 0. 0. " << Ls[1] << " 0. 0. 0. " << Ls[2] << "\" Origin=\"" 
		<< -0.5*Ls[0] << " " << -0.5*Ls[1] << " " << -0.5*Ls[2] << "\" Time="<< time << endl;
	// Print out each molecules coords
	for(uint i=0; i<N; i++){
        molecule &m = *molptrs[i];
		m.write2xyz(xyz, bx);
	}
	return;
}

void simulation::write_theta_horz(ofstream& thetafile, flt time, uint Nseg){
	// Get current box dimensions
	Vector3d Ls = bx.get_Ls();
	// Initialial 4 columns are time Lx Ly Lz
	// Followed by the angle from the horizontal 
	thetafile << time << "\t" << Ls[0] << "\t" << Ls[1] << "\t" << Ls[2];
	for(uint i=0; i<Nseg; i++){
		if (i<N){	
	        molecule &m = *molptrs[i];
	        thetafile << "\t";
			m.write2thetas_horz(thetafile);
		}else{
			thetafile << "\tNAN";
		}
	}
	thetafile << "\n";
	return;
}

void simulation::print_Uijs(){
	for(uint i=0;i<N;i++){
		for(uint j=i+1;j<N;j++){
			if(Uij>1e-8){
				cout << "\t(" << i << ", " << j << "): " << Uij(i,j);
			}
		}
	}
}


flt simulation::get_Uij(flt dist, flt phi, flt coeff, uint i, uint j){
	flt U_rad = 0;
	flt U_rlj = 0;
	flt U_ang = 0;
	flt U_pair = 0;
	//cout <<"dist: " << dist << ", phi: " << phi << ", coef: " << coeff << endl;
	//cout << "sig: " << sig <<", r_cut: " << r_cut << ", sigma: " << sigma << "\n"; 
	if((dist<sig)&&(dist>sigma)){
		flt dist2 = dist*dist;
		flt dist4 = dist2*dist2;
		flt dist3 = pow(dist,3.);
		flt inv_dist6 = 1./(dist4*dist2);
		U_rad = - C0 *(dist4 - 4.*dist3 + C1 * dist2 - C2 * dist + C3);
		U_rlj = eps_rlj * (1. - sig6*inv_dist6) * (1. - sig6*inv_dist6);
		if(((fmod((phi-phi_min+0.5*M_PI), M_PI) - 0.5*M_PI)>0) && ((fmod((phi-phi_max+0.5*M_PI), M_PI) - 0.5*M_PI)<0)){
			U_ang = -0.5*(cos(alpha*(phi-theta0)) + 1.);
		}
	}else if((dist<r_cut)&&(dist>sig)){
		flt dist2 = dist*dist;
		flt dist4 = dist2*dist2;
		flt dist3 = pow(dist,3.);
		U_rad = - C0 *(dist4 - 4.*dist3 + C1 * dist2 - C2 * dist + C3);
		if(((fmod((phi-phi_min+0.5*M_PI), M_PI) - 0.5*M_PI)>0) && ((fmod((phi-phi_max+0.5*M_PI), M_PI) - 0.5*M_PI)<0)){
			U_ang = -0.5*(cos(alpha*(phi-theta0)) + 1.);
		}
	}else if(dist<sigma){
		flt inv_dist6 = 1./pow(dist,6.0);
		U_rlj = eps_rlj * (1. - sig6*inv_dist6) * (1. - sig6*inv_dist6);
	}
	/*
	cout << "dist: " << dist << ", phi: " << phi << " rad., " << phi*180./M_PI << " deg." << endl;
	cout << "phi_min: " << phi_min <<", phi_max: " << phi_max << endl;
	cout << "U_rlj: " << U_rlj << ", U_rad " << U_rad << ", U_ang: " << U_ang << endl;
	*/
	if(BondFlag){
		U_pair = U_rlj - U_rad * U_ang * coeff * (1. + BondRate * BondTime(i,j));
	}else{
		U_pair = U_rlj - U_rad * U_ang * coeff;
	}
	return U_pair;
}

flt simulation::get_Ui_new(uint ind, flt Uw){
	Vector3d dr_vec;
	flt dist, rijsq, phi, coeff, Uijval, dU;
    molecule &mi = *molptrs[ind];
	Vector3d ri = mi.get_r();

	// Initialize U with wall potential energy
	dU = Uwall(params,mi.get_zcom()) - Uw;
	//dU = mg*mi.get_zcom() - Uw;
	/*
	cout << "z: " << mi.get_zcom() << ", mg: " << mg;
	cout << ", U from wall: " << U << endl;
	*/
	for(uint j=0;j<N;j++){
		if(j==ind){
			continue;
		}
		//cout << "(" << ind <<", " << j <<")\n";
		// Get values from molecule j 
    	molecule &mj = *molptrs[j];
    	Vector3d rj = mj.get_r();
    	flt li = mi.get_l();
    	flt lj = mj.get_l();
    	// Only calculate true distance vector if two coms are within 1.5 mean segment length from each other
    	flt first_cutsq = 9./16.*(li+lj)*(li+lj);

    	// As a first check see if boxed distance between rcoms
    	Vector3d rij = rj - ri;
    	rij = bx.dist(rij);
    	rijsq = rij.dot(rij);

    	if(rijsq<first_cutsq){
    		// Get directional vector and calculate displacement vector
	    	Vector3d uj = mj.get_u();

	    	// Get displacement vector and distance between 
	    	//     nearest two points on two line segments
	    	dr_vec = mi.scl_disp(uj, rj, lj, bx, &dist, &coeff);
	    	// Make dr_vec a unit vector
	    	dr_vec = 1./dist*dr_vec;
	    	if((dist<r_cut)&&(coeff>0)){
	    		phi = mi.calc_phi_chiral(uj, dr_vec);
	    	}else{
	    		phi = 0.;
	    		coeff = 0.;
	    	}
	    	// Add incremental potential energy
	    	Uijval = get_Uij(dist, phi, coeff, ind, j);
	    	//Uij[ind,j] = Uijval;
	    	//Uij[j,ind] = Uijval;
	    	if(Uijval!=0){
	    		flt dUval= Uijval - Uij(ind, j);
	    		if(abs(dUval) > 1e-8){
		    		dUs.push_back(dUval);
		    		indsUij.push_back(j);	
	    		}
	    	}

	    	dU += Uijval;
    	}
	}
	return dU;
}


flt simulation::get_Ui(uint ind){
	Vector3d dr_vec;
	flt dist, rijsq, phi, coeff, Uijval;
    molecule &mi = *molptrs[ind];
	Vector3d ri = mi.get_r();

	// Initialize U with wall potential energy
    U = Uwall(params,mi.get_zcom());
	//U = mg*mi.get_zcom();
	/*
	cout << "z: " << mi.get_zcom() << ", mg: " << mg;
	cout << ", U from wall: " << U << endl;
	*/
	for(uint j=0;j<N;j++){
		if(j==ind){
			continue;
		}
		//cout << "(" << ind <<", " << j <<")\n";
		// Get values from molecule j 
    	molecule &mj = *molptrs[j];
    	Vector3d rj = mj.get_r();
    	flt li = mi.get_l();
    	flt lj = mj.get_l();
    	// Only calculate true distance vector if two coms are within 1.5 mean segment length from each other
    	flt first_cutsq = 9./16.*(li+lj)*(li+lj);

    	// As a first check see if boxed distance between rcoms
    	Vector3d rij = rj - ri;
    	rij = bx.dist(rij);
    	rijsq = rij.dot(rij);

    	if(rijsq<first_cutsq){
    		// Get directional vector and calculate displacement vector
	    	Vector3d uj = mj.get_u();

	    	// Get displacement vector and distance between 
	    	//     nearest two points on two line segments
	    	dr_vec = mi.scl_disp(uj, rj, lj, bx, &dist, &coeff);
	    	// Make dr_vec a unit vector
	    	dr_vec = 1./dist*dr_vec;
	    	if((dist<r_cut)&&(coeff>0)){
	    		phi = mi.calc_phi_chiral(uj, dr_vec);
	    	}else{
	    		phi = 0.;
	    		coeff = 0.;
	    	}
	    	// Add incremental potential energy
	    	Uijval = get_Uij(dist, phi, coeff, ind, j);
	    	//Uij[ind,j] = Uijval;
	    	//Uij[j,ind] = Uijval;
	    	U += Uijval;
    	}
	}
	return U;
}

bool simulation::metMC(flt dU, flt z){
	if(z<zmin){
		return false;
	}else if(dU<0){
		return true;
	}else{
		// Update with probability according to Boltzmann factor
		return randUniform()<exp(-dU/KbT);
	}
}

void simulation::timestep(uint i){
    Vector3d r0;
    Vector3d disp_vec;
    Quaterniond q0 = Quaterniond(0.,0.,0.,0.);
    Quaterniond delta_q = Quaterniond(0.,0.,0.,0.);
    flt U0, dU, U;
    molecule &mi = *molptrs[i];

    // Store the original r0 and q0
    r0 = mi.get_r();
    q0 = mi.get_q();
    
    // Collect initial potential energy
    U0 = get_Ui(i);
    //cout << "U0: " << U0 << endl;
    
    ////  Make translational MC move
    // Generate displacement vector and update vector com
    disp_vec = dr.array() * randUnitVec().array();
    // Make move
    mi.update_r(disp_vec);

    // Calculate change in energy due to updated displacement
    U = get_Ui(i);
    if(metMC(U-U0,mi.get_minz())){
    	//cout << "Move accepted\n";
    	U0 = U; // set U0 to U
    	ConsecRej = 0;
    }else{
    	// Move not accepted revert rcom back to old one
    	//cout << "Move failed\n";
    	mi.set_r(r0);
    	ConsecRej++;
    }


    //// Make rotational MC move
    // Generate displacement vector and update vector com
    delta_q.vec() = dq.vec().array() * randUnitVec().array();
    // Make move
    mi.update_q(delta_q);

    // Calculate change in energy due to updated displacement
    U = get_Ui(i);
    if(metMC(U-U0,mi.get_minz())){
    	//cout << "Rotational Move accepted\n";
    	U0 = U; // set U0 to U
    	ConsecRej = 0;
    }else{
    	// Move not accepted revert rcom back to old one
    	//cout << "Rotational Move failed\n";
    	mi.set_q(q0);
    	ConsecRej++;
    }
}
void simulation::init_Uij(){
	cout << "In init_Uij\t";
	for(uint i=0;i<N;i++){
		cout << i << "\t";
		molecule &mi = *molptrs[i];

		flt Uw = Uwall(params,mi.get_zcom());

		// Calculate change in energy due to updated displacement
		flt dU = get_Ui_new(i, Uw);
		update_Uij(i);
	}
	cout << "\n";
}

void simulation::update_Uij(uint i){
	for(uint j=0; j<dUs.size(); j++){
		Uij(i, indsUij[j]) += dUs[j];
		Uij(indsUij[j], i) += dUs[j];
	}
	indsUij.clear();
	dUs.clear();
}

flt simulation::get_Utot(){
	flt Utot = 0.;
	for(uint i=0; i<N; i++){
		molecule &mi = *molptrs[i];
		Utot += Uwall(params,mi.get_zcom());
		for(uint j=i+1; j<N; j++){
			Utot += Uij(i,j);
		}
	}
	return Utot;
}

void simulation::timestep_new(uint i){
    Vector3d r0;
    Vector3d disp_vec;
    Quaterniond q0 = Quaterniond(0.,0.,0.,0.);
    Quaterniond delta_q = Quaterniond(0.,0.,0.,0.);
    flt U0, dU, U;
    molecule &mi = *molptrs[i];


    // Store the original r0 and q0
    r0 = mi.get_r();
    q0 = mi.get_q();
    
    //cout << "U0: " << U0 << endl;
    flt Uw = Uwall(params,mi.get_zcom());
	//flt Uw = mg*mi.get_zcom(); // Gravitational potential energy before move
    
    ////  Make translational MC move
    // Generate displacement vector and update vector com
    disp_vec = dr.array() * randUnitVec().array();
    // Make move
    mi.update_r(disp_vec);

    // Calculate change in energy due to updated displacement
    dU = get_Ui_new(i, Uw);
    if(metMC(dU,mi.get_minz())){
    	//cout << "Move accepted\n";
    	U0 = U; // set U0 to U
    	ConsecRej = 0;
    	update_Uij(i);
    	Uw = Uwall(params,mi.get_zcom());
		//Uw = mg*mi.get_zcom(); // Gravitational potential energy before move
    }else{
    	// Move not accepted revert rcom back to old one
    	//cout << "Move failed\n";
    	mi.set_r(r0);
    	ConsecRej++;
	    dUs.clear();
	    indsUij.clear();
    }


    //// Make rotational MC move
    // Generate displacement vector and update vector com
    delta_q.vec() = dq.vec().array() * randUnitVec().array();
    // Make move
    mi.update_q(delta_q);

    // Calculate change in energy due to updated displacement
    dU = get_Ui_new(i, Uw);
    if(metMC(dU,mi.get_minz())){
    	//cout << "Rotational Move accepted\n";
    	U0 = U; // set U0 to U
    	ConsecRej = 0;
    	update_Uij(i);
    }else{
    	// Move not accepted revert rcom back to old one
    	//cout << "Rotational Move failed\n";
    	mi.set_q(q0);
    	ConsecRej++;
	    dUs.clear();
	    indsUij.clear();
    }
    if(BondFlag){
    	update_BondTime(i);
    }
}

void simulation::update_BondTime(uint i){
	for(uint j=i+1;j<N;j++){
		flt Uijval = Uij(i,j);
	 	if(Uijval<0){
			BondTime(i,j) += 1;
			BondTime(j,i) += 1;
		}else{
			BondTime(i,j) = 0;
			BondTime(j,i) = 0;
		}
    }
}