#ifndef BOX_HPP
#define BOX_HPP
#include "global.hpp"

//using namespace std;

class box{
	// Define boundary vector of strings
    vector<string> bounds;
    Vector3d Ls;

    public:
    	box();
        box(string*, Vector3d);
        void set_bounds(string*);
        void set_Ls(flt*);
        Vector3d get_Ls();
        Vector3d dist(Vector3d);
        void rescale(Vector3d);
        void apply_strain(Vector3d);
};

#endif //BOX_HPP
