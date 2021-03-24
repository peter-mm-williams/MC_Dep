#include <time.h>
#include "simulation.hpp"
#include <getopt.h>
#include <iostream>
#include <fstream>

#define CLOCKS_PER_MS (CLOCKS_PER_SEC / 1000)


int main(int argc, char **argv){

	int c;
	/*
	Set up optional arguments
	*/
	// Declare variables for optional arguments
	flt KbT = 1.; // Temperature for MC algorithm
	flt l_seg = 100.0/3.; // Length of segments
	int NperCSC = 21; // Number of segments on individual line
	int Nseg = 600; // Number of total segments deposited
	flt z_drop = 20.; // Height from which segments are deposited
    string xyzname = "test.xyz"; // outpuf filename
    string inputfileName = "input.coords"; // input file for coordinates
    flt eps = 1000.; // Well depth
    flt e_factor = 1.; // Factor for strain rate
    flt eps_rlj = 1000.; // Pre-factor strength for repulsive lennard-jones
    flt theta0 = 10*M_PI/180.; // Preferred Relative Angle
    flt dtheta = 20*M_PI/180.; // Width of well in angular dimension
    flt mg = 50.; // Strength of wall deposition force
    flt w = 1.0; // Range of attraction in radial direction
    int Nprint = 20000; // Number of MC steps between printing coordinates
	int MaxRejects = 100.; // Number of consecutive rejections before deposition
	flt rottheta = 0.2; // Mean size of quaternion MC step in degrees
	flt dr_factor = 0.05; // Mean size of displacement in MC step in sigma (width of rods)
	uint seedval = 1; // Seed for random number generator
	int maxCnt = 100000; // Max number of steps taken during deposition
	uint Nrelax = 100; // Number of Steps to relax system
	uint Nload = 0; // Initializing the number of molecules loaded from .coords file
	uint Loadfile=0; // Binary value as to whether a .coords file will be provided
    flt z_cut = 0.; // Cutoff above which wall interaction is non-zero
    flt BondRate = 0; // Rate at which well depths grow eps = eps*(1+r/eps*t)

	opterr = 0; 
    c = getopt (argc, argv, "T:l:M:N:z:x:e:E:j:t:d:g:w:P:R:q:r:S:C:L:f:c:B:");
    while(c != -1){
        switch(c){
            case 'T':
                KbT = atof(optarg);
                break;
            case 'l':
                l_seg = atof(optarg);
                break;
            case 'M':
            	NperCSC = atol(optarg);
            	break;
            case 'N':
            	Nseg = atol(optarg);
            	break;
            case 'z':
            	z_drop = atof(optarg);
            	break;
            case 'x':
            	xyzname = optarg;
            	break;
            case 'e':
            	eps = atof(optarg);
            	break;
            case 'E':
            	e_factor = atof(optarg);
            	break;
            case 'j':
            	eps_rlj = atof(optarg);
            	break;
            case 't':
            	theta0 = M_PI/180.0*atof(optarg);
            	break;
            case 'd':
            	dtheta = M_PI/180.0*atof(optarg);
            	break;
            case 'g':
            	mg = atof(optarg);
            	break;
           	case 'w':
           		w = atof(optarg);
           		break;
           	case 'P':
           		Nprint = atol(optarg);
           		break;
           	case 'R':
           		MaxRejects = atol(optarg);
           		break;
           	case 'q':
           		rottheta = atof(optarg);
           		break;
           	case 'r':
           		dr_factor = atof(optarg);
           		break;
           	case 'S':
           		seedval = atol(optarg);
           		break;
           	case 'C':
           		maxCnt = atol(optarg);
           		break;
           	case 'L':
           		Loadfile = atol(optarg);
           		break;
           	case 'f':
           		inputfileName = optarg;
           		break;
            case 'c':
                z_cut = atof(optarg);
                break;
            case 'B':
            	BondRate = atof(optarg);
            	break;
            default:
                abort ();
        }
        c = getopt (argc, argv, "T:l:M:N:z:x:e:E:j:t:d:g:w:P:R:q:r:S:C:L:f:c:B:");
    }
    
    // Seed random number generator    
    if (seedval==1){
        seed(seedval);
    }else{
        seed(seedval);
        srand(seedval);
    }
    
    flt sigma = 1.; // Width of rod

    flt rx, ry;
    flt a,x,y,z;
    a = cos(M_PI/4.);
    x = 0.;
    y = sin(M_PI/4.);
    z = 0.;
    Quaterniond q_arr = Quaterniond(a,x,y,z);

    // Set dr and dq
    Vector3d dr = Vector3d(dr_factor,5.*dr_factor,dr_factor);

    Quaterniond dq = Quaterniond(cos(0.5*rottheta)*cos(0.5*rottheta), 0.5*sin(rottheta),0.5*sin(rottheta),-sin(0.5*rottheta)*sin(0.5*rottheta));

    // Set Wall bounds and Ls
    string bounds[NDIM];
    for(uint i=0;i<NDIM-1;i++){bounds[i]="p";};
    bounds[NDIM-1] = "f";
    Vector3d Ls = Vector3d(70., 40., 100.);

    // Define strain rate
    Vector3d e_strain = e_factor*Vector3d(2.69e-4, 5.79e-4,0);

    // Open input files
    
    ifstream infile;
    if(Loadfile==1){
    	infile.open(inputfileName.c_str(), ios::in);
    	if (infile.fail()){
	    	cout << "Failed to open: " << inputfileName << endl;
	    }else{
	    	infile >> Nload >> Ls[0] >> Ls[1] >> Ls[2];
	    }
    }
    
    // Open output files
    ofstream xyzfile;
    xyzfile.open(xyzname.c_str(), ios::out);

    string thetaname = xyzname.substr(0,xyzname.size()-4)+"_thetas.dat";
    ofstream thetafile;
    thetafile.open(thetaname.c_str(), ios::out);
    thetafile << "time";
    for(uint i=0;i<Nseg;i++){
    	thetafile << "\t" << "theta" << i;
    }
    thetafile<<"\n";

	/* -----------------------------------Set Up Mol Vector and Simulation------------------------------------------------------ */
    
    // Allocate space for molecule vectors
    vector<molecule*> molptrs;
    molptrs.reserve(Nseg);
    vector<molecule> mols;
    mols.reserve(Nseg);

    // Fill mols
    for(uint i=0;i<Nseg;i++){
    	Quaterniond q_new;
    	Vector3d rcom;
    	if(i<Nload){
    		
    		infile >> rcom[0] >> rcom[1] >> rcom[2] >> a >> x >> y >>z;
    		Quaterniond q_load = Quaterniond(a,x,y,z);
    		q_new = q_load;
	    	cout << "Loaded in molecule " << i << "\n";
    		
    	}else{
    		rcom = randVec(); // Initialize random unitvector
	    	rx = l_seg*i - Ls[0]*round(l_seg*i/Ls[0]);
	    	if((i%NperCSC)==0){
	    		ry = Ls[1]*randUniform(); //  Set to new random y value
	    		ry = ry - Ls[1]*round(ry/Ls[1]);
			}
			q_new = q_arr;
	    	rcom = Vector3d(rx, ry, z_drop);
    	}
    	molecule mol = molecule(l_seg, rcom, q_new, sigma);
    	//mol.print_rq2screen();
    	
    	mols.push_back(mol);
    }
    if(Loadfile==1){
    	infile.close();
    }

    // Fill molptrs
    for(uint i=0;i<Nseg;i++){
	    molecule *m = &mols[i];
	    molptrs.push_back(m);
    }
    vector<flt> params;
    params.push_back(z_cut);
    params.push_back(mg);
    // Initialize simulation
    simulation sim = simulation(molptrs, ConstFieldPieceU, eps_rlj, eps, sigma, w, theta0, dtheta, params, KbT, dr, dq, e_strain, bounds, Ls, Nseg, BondRate);
    sim.set_N(Nload);
    sim.init_Uij();
    /* -----------------------------------Run Simulation------------------------------------------------------ */
    // set timer
    clock_t t0 = clock();
    clock_t t1;
    int time = 0;
    for(uint i=Nload;i<Nseg;i++){
    	sim.set_N(i+1);

    	// Strain ics of molecule to be placed
    	molecule &mol = *molptrs[i];
    	//molecule mol = mols[i];
    	mol.n_strain(e_strain, i);
    	sim.reset_Rej();
    	int cnt = 0; // cnt number of steps for single deposition
    	// Place molecule one by one
    	while((sim.get_Rej()<MaxRejects)&&(cnt<maxCnt)){
    		// increment cnt and time
    		cnt++;
    		time++;
    		// take timestep
    		sim.timestep_new(i);
    		// Print out in accordance to specified frequency
    		if((time%Nprint)==0){
                cout << "Placing moelcule: " << i <<", U of placed mol: "  << sim.get_Ui(i) << ", " << sim.get_Utot() << "\n";
    			sim.write_xyzfile(xyzfile,time);
    			sim.write_theta_horz(thetafile, time, Nseg);
    		}
    	}
    	cout << "Placed moelcule: " << i <<", U of placed mol: "  << sim.get_Ui(i);
    	cout << ", theta: " << mol.calc_theta_horizontal()*180./M_PI << "\n\n";
    	// Apply strain
    	if(e_factor>0){
    		sim.apply_strain();
    	}
    	// Run relaxation timesteps
    	for(uint j=0; j<i*Nrelax;j++){
    		time++;
    		int ind = round(randUniform()*i);
    		sim.timestep_new(ind);
    		// Print out in accordance to specified frequency
    		if((time%Nprint)==0){
    			sim.write_xyzfile(xyzfile,time);
    			sim.write_theta_horz(thetafile, time, Nseg);
    		}
    	}
    }

    t1 = clock();
    cout << "Total time: " << (t1-t0)/CLOCKS_PER_SEC << ", time per step: " << (t1-t0)/time/CLOCKS_PER_SEC << endl;
    
    // One last print statement
	sim.write_xyzfile(xyzfile,time);
	sim.write_theta_horz(thetafile, time, Nseg);

    // Close output files
    xyzfile.close();
    thetafile.close();
    
    /*
    // Initialize MPI environment
    MPI_Init(NULL, NULL);

    // Get number of processes
    int Nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);

    // Get rank of process
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
	*/


}
