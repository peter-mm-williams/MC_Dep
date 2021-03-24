#include "molecule.hpp"

molecule::molecule(flt l_seg, Vector3d rcom, Quaterniond q_arr, flt width){
    N=2;
    r=rcom;
    for(uint i=0;i<3;i++){
        r[i] = rcom[i];
    }
    //cout << "r: " << r[0] << "," << r[1] << "," << r[2] << endl;
    rs.push_back(Vector3d (0,0,-0.5*l_seg)); // Coordinates of endpoints of segment
    rs.push_back(Vector3d (0,0,0.5*l_seg)); // 
    rs0.push_back(Vector3d (0,0,-0.5*l_seg));
    rs0.push_back(Vector3d (0,0,0.5*l_seg));
    q=q_arr;
    l=l_seg;
    sigma=width;
    Vector3d u_i(0.0,0.0,1.0);
    u=u_i;
    rotRs2W();
}

void molecule::rotRs2W(){
    // Rotate coords to world frame
    for(uint i=0; i<N; i++){
        Quaterniond q_r;
        q_r.w() = 0;
        q_r.vec() = rs0[i];
        Quaterniond rot_r = q* q_r * q.inverse();
        rs[i] =  rot_r.vec()+ r;
    }
    Quaterniond q_u;
    Vector3d u_i(0.0,0.0,1.0);
    q_u.w()=0;
    q_u.vec()=u_i;
    u = (q* q_u * q.inverse()).vec();
    return;
}

void molecule::rotRs2M(){
    // translate coordinates from world frame to molecule frame
    for(uint i=0; i<N; i++){
        Quaterniond q_r;
        q_r.w() = 0;
        q_r.vec() = rs[i] - r;
        Quaterniond rot_r = q.inverse()* q_r * q;
        rs[i] =  rot_r.vec();
    }
    return;
}

void molecule::update_r(Vector3d dr){
    // update com vector
    r = r + dr;
    // update rs 
    for(uint i=0;i<N;i++){
        rs[i] = rs[i] + dr;
    }
    return;
}

void molecule::update_q(Quaterniond dq){
    // Perform displacement to quaternion and renormalize
    q.w() += dq.w();
    q.x() += dq.x();
    q.y() += dq.y();
    q.z() += dq.z();
    q.normalize();
    rotRs2W(); // update rs and u based upon new quaternion
    return;
}

void molecule::set_q(Quaterniond q_new){
    // Replace quaternion with new quaternion and update rs and u
    q.w() = q_new.w();
    q.vec() = q_new.vec();
    q.normalize();
    rotRs2W();
    return;
}

void molecule::set_r(Vector3d r0){
    r = r0;
    return;
}

void molecule::strain(Vector3d e_strain){
    // Apply straint to molecule
    Vector3d e_factor = Vector3d(1.,1.,1.)+e_strain;
    r = r.array()*e_factor.array();
    flt theta0 = atan(u[1]/u[0]);
    flt dtheta = atan((e_strain[1]+1.)/(e_strain[0]+1.)*tan(theta0)) - theta0;
    Quaterniond q_rot = Quaterniond(cos(0.5*dtheta),0.,0.,sin(0.5*dtheta));
    q = q_rot*q;
    rotRs2W();
    return;
}

void molecule::n_strain(Vector3d e_strain, int n){
    // Apply straint to molecule
    Vector3d e_factor = Vector3d(0.,0.,0.);
    for(uint i=0;i<NDIM;i++){
        e_factor[i] = pow(1. + e_strain[i], n);
    }
    r = r.array()*e_factor.array();
    flt theta0 = atan(u[1]/u[0]);
    flt dtheta = atan((e_strain[1]+1.)/(e_strain[0]+1.)*tan(theta0)) - theta0;
    Quaterniond q_rot = Quaterniond(cos(0.5*dtheta),0.,0.,sin(0.5*dtheta));
    q = q_rot*q;
    rotRs2W();
    return;
}

flt molecule::get_zcom(){
    return r[2];
}

flt molecule::get_minz(){
    return min(rs[0][2],rs[1][2]);
}

void molecule::print_rq2screen(){
    cout<< r[0] << ", " << r[1] << ", " << r[2] << ", " << q.w() << ", " << q.x() << ", " << q.y() << ", " << q.z() << "\n"; 
    return;
}

void molecule::print_u2screen(){
    cout << u[0] << ", " << u[1] << ", " << u[2] << endl;
    return;
}

flt molecule::calc_theta_horizontal(){
    return atan(u[1]/u[0]);
}

void molecule::write2xyz(ofstream& xyz, box bx){
    flt theta = this->calc_theta_horizontal();
    // Set 3 color values for particles (scale blue with angle relative to horizontal)
    flt color_b = (theta+M_PI/2.)/M_PI;
    flt color_g = 0.5;
    flt color_r = 0.1;
    Vector3d r_box = bx.dist(r) + 0.5*bx.get_Ls();
    // Print out particle type, r_com, radius of spherocylinder, and length of spherocylinder
    xyz << "X\t" << r_box[0] << "\t" << r_box[1] << "\t" << r_box[2] << "\t0.5\t" << l << "\t";
    // Print out quaternion and color triplet:
    xyz << q.x() << "\t" << q.y() << "\t" << q.z() << "\t" << q.w() <<"\t" << color_r << "\t" << color_g << "\t" << color_b << endl;
    return;
}

void molecule::write2thetas_horz(ofstream& thetasfile){
    flt theta = this->calc_theta_horizontal();
    thetasfile << theta;
    return;
}

Vector3d molecule::get_u(){
    return u;
}


Vector3d molecule::get_r(){
    return r;
}

flt molecule::get_l(){
    return l;
}


Quaterniond molecule::get_q(){
    return q;
}

Vector3d molecule::scl_disp(Vector3d uj, Vector3d rj, flt lj, box bx, flt* dr, flt* coeff){
    Vector3d ui, ri, rij;
    Vector3d dr_v;
    //Vector3d &dr_v = dr_vec;
    flt uij, rui, ruj, A, lambda_i, lambda_j, li, delta_i, delta_j, rijsq;
    int i_case, j_case;
    i_case = 1;
    j_case = 1;
    // Extract u and r for the molecule
    ui = u;//get_u();
    ri = r;//get_r();
    li = l;

    // Calculate boxed displacement vector of coms
    rij = rj - ri;
    rij = bx.dist(rij);

    // Calculate dot product of two orientation vectors
    uij = ui.dot(uj);
    rui = rij.dot(ui);
    ruj = rij.dot(uj);
    A = 1./(1.-uij*uij); // Calculate pre-factor for lambdas

    // Calculate initial scaling factors for nearest points on two segments
    lambda_i = A*(rui - ruj*uij);
    lambda_j = -A*(ruj - rui*uij);

    // Assess special case where segments are nearly parallel
    if(abs(uij*uij-1.)<1e-8){
        flt si, sj, li_mid, li_end, lj_mid, lj_end;
        if(abs(rui)<1e-8){
            si = 0.;
        }else{
            si = rui/abs(rui);
        }
        if(abs(ruj)<1e-8){
            sj = 0.;
        }else{
            sj = -ruj/abs(ruj);
        }
        li_mid = rui + 0.5*sj*lj;
        li_end = 0.5*si*li;
        lj_end = 0.5*sj*lj;
        lj_mid = -ruj + 0.5*si*li;
        i_case = 4;
        j_case = 4;

        lambda_i = 0.5*(li_mid + li_end);
        lambda_j = 0.5*(lj_end + lj_mid);

        if(abs(lj_mid)>0.5*lj){
            lambda_j = 0;
            lambda_i = rui;
            j_case = 5;
            i_case = 6;
        }else if(abs(li_mid)>0.5*li){
            lambda_i = 0;
            lambda_j = -ruj;
            j_case = 6;
            i_case = 5;
        }
    }


    // Calculate distance from 
    delta_i = abs(lambda_i) - 0.5*li;
    delta_j = abs(lambda_j) - 0.5*lj;

    if((delta_i>0)||(delta_j>0)){
        if((delta_i>delta_j)){
            lambda_i = 0.5*li*lambda_i/abs(lambda_i); // Set i to endpoint
            i_case = 2;
            lambda_j = -ruj+lambda_i*uij;
            j_case = 3;

            if(abs(lambda_j)>0.5*lj){
                lambda_j = 0.5*lj*lambda_j/abs(lambda_j); // Endpoint
                j_case = 2;
            }
        }else{
            lambda_j = 0.5*lj*lambda_j/abs(lambda_j);
            j_case = 2;
            lambda_i = rui + lambda_j*uij;
            i_case = 3;

            if(abs(lambda_i)>0.5*li){
                lambda_i = 0.5*li*lambda_i/abs(lambda_i);
                i_case = 2;
            }
        }
    }
    /*
    cout << "i_case: " << i_case << ", j_case: " << j_case << "\n";
    cout << "lambda_i: " << lambda_i << ", lambda_j: " << lambda_j << ", rui: " << rui << ", ruj: " << ruj << ", uij: " << uij << endl;
    // Calculate distance squared and displacement vector between nearest two points on segment
    cout << " tanh((|lam1|-l1/2)^2): " << tanh((lambda_i*lambda_i - 0.25*li*li)*(lambda_i*lambda_i - 0.25*li*li));
    cout << ", tanh((|lam2|-l2/2)^2): " << tanh((lambda_j*lambda_j - 0.25*lj*lj)*(lambda_j*lambda_j - 0.25*lj*lj)) << endl;
    cout << ", Argument 1: "<< (lambda_i*lambda_i - 0.25*li*li)*(lambda_i*lambda_i - 0.25*li*li) << endl;
    cout << ", Argument 1 take 2: "<< abs(abs(lambda_i) - 0.5*li) << ", tanh(arg1_2) : " << tanh(abs(abs(lambda_i) - 0.5*li)) << endl;
    cout << ", Argument 2: "<< abs(abs(lambda_j) - 0.5*lj) << endl;
    cout << "new_coeff: " << tanh(2.*abs(abs(lambda_j) - 0.5*lj))*tanh(2.*abs(abs(lambda_i) - 0.5*li)) << endl;
    */
    rijsq = rij.dot(rij);
    //*dr = sqrt(rijsq + lambda_i*lambda_i + lambda_j*lambda_j + 2.0*lambda_j*ruj - 2.0*lambda_i*rui - 2.0*lambda_i*lambda_j*uij);
    dr_v = rij + lambda_j*uj - lambda_i*ui;
    dr_v = bx.dist(dr_v);
    *dr = sqrt(dr_v.dot(dr_v));
    *coeff = tanh(2.*abs(abs(lambda_j) - 0.5*lj))*tanh(2.*abs(abs(lambda_i) - 0.5*li));
    return dr_v;
}

flt molecule::calc_phi_chiral(Vector3d uj, Vector3d dr_vec){
    Vector3d ui, u_perp, u_plane;
    flt ujup, sqrtcoeff;
    ui = u;
    // Calc unit vector perpandicular to ui and dr_vec
    u_perp = dr_vec.cross(ui);
    u_perp = u_perp*1./sqrt(u_perp.dot(u_perp));

    // Calc dot product of 
    ujup = uj.dot(u_perp);
    /*
    cout << "ui dot uj: " << ui.dot(uj) << ", ui dot -uj: " << ui.dot(-uj) << endl;
    cout << "ui cross uj: \n" << ui.cross(uj) << "\n  ui cross -uj\n" << ui.cross(-uj) << endl; 
    cout << "dr cross ui: \n" << dr_vec.cross(ui) << "\n  dr cross uj\n" << dr_vec.cross(uj) << endl; 
    // Calc unit vector of plane defined by unit vectors of two vectors
    if ((ui.dot(uj)-1.)<1e-8){
        u_plane = ui.cross(-uj+ Vector3d(1.,1e-4,1e-4));
    }else{
        u_plane = ui.cross(uj);
    }
    u_plane = u_plane*1./sqrt(u_plane.dot(u_plane));

    sqrtcoeff = dr_vec.dot(u_plane)*dr_vec.dot(u_plane);
    *coeff = sqrtcoeff*sqrtcoeff;
    */
    return atan(ujup/sqrt(1-ujup*ujup));
}