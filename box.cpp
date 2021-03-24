#include "box.hpp"

box::box(){
    string bounds[NDIM];
    for(uint i=0;i<NDIM;i++){bounds[i]="p";};
    set_bounds(bounds);
    Ls = Vector3d(0.,0.,0.);
}

box::box(string* bound_arr, Vector3d Ls0){
    set_bounds(bound_arr);
    for(uint i=0; i<NDIM; i++){
        Ls[i] = Ls0[i];
    }
}

void box::set_bounds(string* bound_arr){
    for(uint i=0; i<NDIM; i++){
        bounds.push_back(bound_arr[i]);
    }
}

void box::set_Ls(flt* Larr){
    for (uint i=0;i<NDIM;i++){
        this->Ls[i]=Larr[i];
    }
}

Vector3d box::get_Ls(){
    // Getter for box dimensions
    return Ls;
}

Vector3d box::dist(Vector3d rdiff){
    for (uint i=0; i<NDIM; i++){
        if (bounds[i]=="p"){
            rdiff[i]= rdiff[i] - Ls[i]*round(rdiff[i]/Ls[i]);
        }
    }
    return rdiff;
}

void box::apply_strain(Vector3d e_strain){
    Ls = Ls.array()*(Vector3d(1.,1.,1.)+e_strain).array();
}

void box::rescale(Vector3d Lnew){
    Ls = Lnew;
    /*
    Vector3d Lrat = Lnew.array()/Ls.array();
    Vector3d dL = Lrat.array() - 1;
    Vector3d dR;
    for(uint i=0; i<molptrs.size(); i++){
        molecule &m = *molptrs[i];
        dR = dist(m.r).array()*dL.array();
        // shift location of all atoms in each molecule
        for(uint j=0; j<m.N; j++){
            m.rs[j]+=dR;
        }
        // Shift com position vector of molecule
        m.r+=dR;
    }
    */
}

flt ConstFieldU(vector<flt> params, flt z){
    /*
    Apply constant field to a dimension
    k = params[0]
    U = k*z
    */
    return params[0]*z;
}

flt ConstFieldPieceU(vector<flt> params, flt z){
    /*
    cut = params[0] 
    k = params[1]

    U = {
    0          for z<cut
    k*(z-cut)  for z>cut
    }
    */
    if(z<=params[0]){
        return 0;
    }
    else{
        return params[1]*(z-params[0]);
    }
}