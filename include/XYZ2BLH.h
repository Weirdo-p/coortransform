#include "common.h"
#include "utils.h"

#ifndef _XYZ2BLH_H_
#define _XYZ2BLH_H_

class XYZ2BLH
{
    private:
        MatrixNNd point;
        double W;
        double N;
        double e;
        double B;
        double L;
        double H;
    
    public:
        XYZ2BLH(string datapath, string savepath);
        XYZ2BLH();

        double Sete();
        double SetW(double B);
        double RAD2DEG(double rad);
        double DEG2RAD(double deg);
        double SetB(double X, double Y, double Z);
        double SetL(double X, double Y);
        double SetH(double X, double Y, double Z, double B);
        double SetN(double W);
        MatrixNNd GetPoint();
        void WriteToFile(string savepath, MatrixNNd result);

        friend MatrixN1d SetData(string path);

};

#endif