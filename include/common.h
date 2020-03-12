/**
 * last edited on 10th Mar, xz
 * this is the head files used to define necessary dependencies
 * rules:
 * e.g.
 * vector means type std::vector
 * Vector means type Eigen::Vector
 * NN means matrix or V(v)ector with n rows and n columns
 * d means all the items in that specific matrix are double types
*/

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <fstream>


using namespace Eigen;
using namespace std;

typedef vector<Vector3d, Eigen::aligned_allocator<Vector3d>> vector31d;
typedef Matrix<double, Dynamic, 4> MatrixN4d;
typedef Matrix<double, 4, 4> Matrix44d;
typedef Matrix<double, 4, 1> Matrix41d;
typedef Matrix<double, Dynamic, 1> MatrixN1d;
typedef Matrix<double, Dynamic, Dynamic> MatrixNNd;