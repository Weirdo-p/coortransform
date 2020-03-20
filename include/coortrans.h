#include "common.h"

#ifndef _COORTRANS_
#define _COORTRANS_
/**
 * last edited on 13 March
 * base class provides functions and datas may use in 4/7/13 params transforming
 * *****************************************************
 * pointA    origin data
 * pointB    target data
 * P         weight matrix
 * savepath   path to save file
 * L         matrix l in Bx - l = 0;
 * ************************************************
 */
class coortrans
{
    protected:
        MatrixN1d pointA;
        MatrixN1d pointB;
        MatrixNNd P;
        string savepath;
        MatrixN1d L;
        MatrixN1d V;
        MatrixN1d residual;
        int n;
    public:
        virtual MatrixNNd SetP(int n) = 0;
        virtual MatrixN1d SetV() = 0;
        virtual MatrixN1d SetResidual(int n) = 0;
        virtual void WriteToFile(string savepath) = 0;

};
#endif