#include <iostream>
#ifndef GLOBAL_HPP
#define GLOBAL_HPP
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <Eigen/Dense>
#include <math.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <mpi.h>

using namespace Eigen;

using namespace std;

typedef unsigned int uint;
typedef long double flt;

typedef boost::mt19937 engine;
typedef boost::normal_distribution<flt> normdistribution;
typedef boost::uniform_01<flt, flt> lindistribution;
typedef boost::variate_generator<engine &, normdistribution> normgenerator;
typedef boost::variate_generator<engine &, lindistribution> lingenerator;

//int Nprocs;
//int ME;

extern uint NDIM;

unsigned int seed(unsigned int n);
unsigned int seed();

//Prototypes for useful function gnee
flt randUniform();
Vector3d randVec();
Vector3d randNormVec();
Vector3d randUnitVec();
flt ConstFieldU(vector<flt>, flt);
flt ConstFieldPieceU(vector<flt>, flt);



#endif // GLOBAL_HPP