#include "../include/XYZ2NEU.h"
#include <iomanip>
#include <math.h>
using namespace std;

XYZ2NEU::XYZ2NEU(string datapath, string savepath, MatrixNNd Cpoint)
{
    this->helper = XYZ2BLH(datapath, savepath);
    this->point = helper.GetPoint();
    this->savepath = savepath;
    this->Cxyz = Cpoint;
    double B = helper.SetB(Cxyz(0, 0), Cxyz(1, 0), Cxyz(2, 0));
    double L = helper.SetL(Cxyz(0, 0), Cxyz(1, 0));
    double H = helper.SetH(Cxyz(0, 0), Cxyz(1, 0), Cxyz(2, 0), B);
    Cblh.resize(3, 1);
    this->Cblh(0, 0) = B;
    this->Cblh(1, 0) = L;
    this->Cblh(2, 0) = H;
    this->rotation = this->SetRotation();
}

XYZ2NEU::XYZ2NEU(string datapath, string savepath)
{
    this->helper = XYZ2BLH(datapath, savepath);
    this->savepath = savepath;
    this->point = helper.GetPoint();
    this->Cxyz.resize(3, 1);
    for(int i = 0; i < Cxyz.rows(); ++i)
        this->Cxyz(i, 0) = this->point(i, 0);
    double B = helper.SetB(Cxyz(0, 0), Cxyz(1, 0), Cxyz(2, 0));
    double L = helper.SetL(Cxyz(0, 0), Cxyz(1, 0));
    double H = helper.SetH(Cxyz(0, 0), Cxyz(1, 0), Cxyz(2, 0), B);
    Cblh.resize(3, 1);
    this->Cblh(0, 0) = B;
    this->Cblh(1, 0) = L;
    this->Cblh(2, 0) = H;
    this->rotation = this->SetRotation();
    this->NEU = this->Compute();
    WriteToFile(savepath);
}

MatrixNNd XYZ2NEU::SetRotation()
{
    MatrixNNd rotation(3, 3);
    double B = this->Cblh(0, 0);
    double L = this->Cblh(1, 0);
    B = helper.DEG2RAD(B);
    L = helper.DEG2RAD(L);
    for(int i = 0; i < rotation.rows(); ++i)
    {
        for(int j = 0; j < rotation.cols(); ++j)
        {
            if((i + 1) % 3 == 1)
            {
                switch(j)
                {
                    case 0: rotation(i, j) = -sin(B) * cos(L); break;
                    case 1: rotation(i, j) = -sin(B) * sin(L); break;
                    case 2: rotation(i, j) = cos(B); break;
                    default: break;
                }
            }
            else if((i + 1) % 3 == 2)
            {
                switch(j)
                {
                    case 0: rotation(i, j) = -sin(L); break;
                    case 1: rotation(i ,j) = cos(L); break;
                    case 2: rotation(i, j) = 0; break;
                    default: break;
                }
            }
            else
            {
                switch(j)
                {
                    case 0: rotation(i, j) = cos(B) * cos(L); break;
                    case 1: rotation(i, j) = cos(B) * sin(L); break;
                    case 2: rotation(i, j) = sin(B); break;
                    default: break;
                }
            }
        }
    }
    return rotation;
}

MatrixNNd XYZ2NEU::Compute()
{
    MatrixNNd rotation = this->rotation;
    int n = this->point.cols();
    MatrixNNd xyz(3, n);
    for(int j = 0; j < xyz.cols(); ++j)
        for(int i = 0; i < xyz.rows(); ++i)
            xyz(i ,j) = this->Cxyz(i, 0);

    MatrixNNd NEU = rotation * (this->point - xyz);
    return NEU;
}

void XYZ2NEU::WriteToFile(string savepath)
{
    fstream out(savepath, ios::out);
    if(!out)
    {
        cerr << "无法写入文件！" << endl;
        exit(1);
    }

    MatrixNNd NEU = this->NEU;
    for(int i = 0; i < NEU.cols(); ++i)
    {
        out << setprecision(15) << i + 1 << "," << NEU(0, i) << "," << NEU(1, i)<< "," << NEU(2, i) << "\n";
    }
    out.close();

}

