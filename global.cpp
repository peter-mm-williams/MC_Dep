#include "global.hpp"


uint NDIM=3;
engine randengine;
normdistribution mynormaldistribution(0, 1);
normgenerator gauss(randengine, mynormaldistribution);
lindistribution mylineardistribution;
lingenerator uniformrand(randengine, mylineardistribution);

unsigned int seed(unsigned int n) {
    randengine.seed(n);
    return n;
}

unsigned int seed() {
    unsigned int n = static_cast<unsigned int>(time(0));
    randengine.seed(n);
    return n;
}

flt randUniform(){
	return uniformrand();
}

Vector3d randVec(){
    return Vector3d(uniformrand(),uniformrand(), uniformrand());
}

Vector3d randNormVec(){
    return Vector3d(gauss(), gauss(), gauss());
}

Vector3d randUnitVec(){
    Vector3d r = Vector3d(gauss(), gauss(), gauss());
    flt r_norm = sqrt(r.dot(r));
    r = 1./r_norm*r;
    return r;
}
