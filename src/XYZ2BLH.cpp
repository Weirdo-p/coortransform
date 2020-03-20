#include "../include/XYZ2BLH.h"
#include <iomanip>
#include <math.h>
using namespace std;

#define alpha 1 / 298.257223563
#define a 6378137.0

XYZ2BLH::XYZ2BLH()
{
    
}
XYZ2BLH::XYZ2BLH(string datapath, string savepath)
{
    MatrixN1d point = SetData(datapath);
    this->point = point;
    this->point.resize(3, point.rows()/3);
    this->e = this->Sete();
}

double XYZ2BLH::Sete()
{
    double e = 2 * alpha - alpha*alpha;
    return e;
}

double XYZ2BLH::SetW(double B)
{
    B = DEG2RAD(B);
    double sinB2 = sin(B) * sin(B);
    double e2 = e * e;
    double W = sqrt(1 - e2 * sinB2);
    return W;
}

double XYZ2BLH::SetN(double W)
{
    double N = a / W;
    return N;
}

double XYZ2BLH::RAD2DEG(double rad)
{
    double sec = rad * 206265;
    double deg = sec / 3600.0;
    return deg;
}

double XYZ2BLH::DEG2RAD(double deg)
{
    double sec = deg * 3600;
    double rad = sec / 206265.0;
    return rad;
}

double XYZ2BLH::SetL(double X, double Y)
{
    double L;
    if(X != 0)
    {
        L = atan(Y / X);
        L = this->RAD2DEG(L);
        if(X > 0 && Y < 0)
            L += 360;
        else if(X < 0 && Y > 0)
            L += 180;
        else if(X < 0 && Y < 0)
            L = abs(L) + 180;
    }
    else
        if(Y > 0)
            L = 90;
        else if(Y < 0)
            L = 180;
        else
            L = 0;

    return L;
}

double XYZ2BLH::SetB(double X, double Y, double Z)
{
    double B;
    double B0 = 1;
    while(true)
    {
        double W = this->SetW(B0);
        double N = this->SetN(W);
        double H = this->SetH(X, Y, Z, B0);
        double e2 = this->e * this->e;
        B0 = this->DEG2RAD(B0);
        double upper = Z + N * e2 *sin(B0);
        double under = sqrt(X * X + Y * Y);
        B = atan(upper / under);
        B0 = this->RAD2DEG(B0);
        B = RAD2DEG(B);

        double error = abs(B - B0);
        if(error > 1e-7)
            B0 = B;
        else
            break;
    }
    return B;
}

double XYZ2BLH::SetH(double X, double Y, double Z, double B)
{
    double W = this->SetW(B);

    double N = this->SetN(W);
    double e2 = this->e * this->e;
    double upper = Z;
    B = DEG2RAD(B);
    double under = sin(B);
    double H = upper / under - N * (1 - e2);
    return H;
}

MatrixNNd XYZ2BLH::GetPoint()
{
    return this->point;
}

void XYZ2BLH::WriteToFile(string savepath, MatrixNNd result)
{
    fstream out(savepath, ios::out);
    if(!out)
    {
        cerr << "无法写入文件！" << endl;
        exit(1);
    }

    for(int i = 0; i < result.cols(); ++i)
    {
        out << setprecision(15) << i + 1 << "," << result(0, i) << "," << result(1, i)<< "," << result(2, i) << "\n";
    }
    out.close();
}