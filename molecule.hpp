#ifndef MOLECULE_HPP
#define MOLECULE_HPP
#include "box.hpp"
#include <iostream>
#include <fstream>

class molecule{
    uint N;
    Quaterniond q;
    Vector3d r;
    Vector3d u;
    vector<Vector3d> rs;
    vector<Vector3d> rs0;
    flt l,sigma;

    public:
        // objects
        //vector<uint> neighbors, new_neighbors;
        // Member functions
        molecule(flt, Vector3d, Quaterniond, flt);
        void rotRs2W();
        void rotRs2M();
        void update_q(Quaterniond);
        void set_q(Quaterniond);
        void set_r(Vector3d);
        void update_r(Vector3d);
        void strain(Vector3d);
        void n_strain(Vector3d, int);
        flt get_zcom();
        flt get_minz();
        void print_rq2screen();
        void print_u2screen();
        flt calc_theta_horizontal();
        void write2xyz(ofstream&, box);
        void write2thetas_horz(ofstream&);
        Vector3d get_u();
        Vector3d get_r();
        flt get_l();
        Quaterniond get_q();
        Vector3d scl_disp(Vector3d, Vector3d, flt, box, flt*, flt*);
        flt calc_phi_chiral(Vector3d uj, Vector3d dr_vec);
};

#endif // MOLECULE_H
