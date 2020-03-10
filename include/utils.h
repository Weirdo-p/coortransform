#include <iostream>
#include <string>
#include <fstream>
#include "common.h"



/**
 * last edited at 9th March by xz
 * function    to judge if index out of range
 * 
 * **********************************************************************************************
 * param i     index of row
 * param j     index of column
 * param n     range of column
 * **********************************************************************************************
 * 
 */
void isOut(int &i, int &j, int num = 2);


/**
 * last edited at 9th March by xz
 * function    to judge if index out of range
 * 
 * **********************************************************************************************
 * param a     the one we want to convert to double
 * return      result
 * **********************************************************************************************
 * 
 */
double string2double(string a);

/**
 * function    to get how many lines in the data file which determines rows in matrix
 * 
 * ***********************************************************
 * param path  data file path
 * ***********************************************************
 *
 */
int GetLines(string path);


/**
 * last edited at 9th March by xz
 * function    to get coordinate data from .xyz files
 * 
 * ************************************************************************************************
 * param path  the path of the data file
 * param data  to store data from the file and also the returns from this function
 * param n     how many lines you want to get in the file. you will get all data when it is -1.
 * ************************************************************************************************
 * 
 */
vector31d GetData(string path);


