#include "common.h"


#ifndef _UTILS_H_
#define _UTILS_H_
/**
 * last edited at 9th March by xz
 * function    to judge if index out of range
 * 
 * **************************************************
 * param a     the one we want to convert to double
 * return      double type number
 * **************************************************
 * 
 */
double string2double(string a);

/**
 * last edited at 9th March by xz
 * function    to get coordinate data from .xyz files
 * 
 * ***************************************************
 * param path  the path of the data file
 * ***************************************************
 */
MatrixN1d SetData(string path);


#endif