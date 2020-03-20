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
#ifndef _COMMON_H_
#define _COMMON_H_
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

typedef Matrix<double, Dynamic, 7> MatrixN7d;
typedef Matrix<double, 7, 1> Matrix71d;

typedef Matrix<double, Dynamic, 13> MatrixN13d;
typedef Matrix<double, 13, 1> Matrix131d;
#endif