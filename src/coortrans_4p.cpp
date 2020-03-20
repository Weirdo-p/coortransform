#include "../include/coortrans_4p.h"
#include <iomanip>
#include <math.h>

coortrans_4p::coortrans_4p(string pathA, string pathB, const int n, string savepath)
{
    cout << "开始初始化" << endl;
    this->n = n;
    this->savepath = savepath;
    cout << "读取坐标文件中..." << endl;
    this->pointA = SetData(pathA);
    this->pointB = SetData(pathB);
    cout << "计算参数估计值..." << endl;
    this->EstiParam = this->SetEstiParam();
    cout << "计算B矩阵..." << endl;
    this->B = this->SetB(n);
    cout << "初始化P矩阵..." <<endl;
    this->P = this->SetP(n);
    cout << "初始化L阵..." <<endl;
    this->L = this->SetL(n);
    cout << "计算结果..." <<endl;
    this->x0 = this->Setx0();
    this->V = this->SetV();
    this->residual = this->SetResidual(n);
    cout << "四参数转换成功！ 相关数据已经写入    " << savepath << endl;
    this->WriteToFile(this->savepath);
}

MatrixN4d coortrans_4p::SetCoefficient(MatrixN1d pointA)
{
    int rows = pointA.rows();
    MatrixN4d coefficient(rows, 4);
    for(int i = 0; i < pointA.rows(); ++i)
    {
        for(int j = 0; j < 4; j++)
        {
            if((i+1) % 2 == 1)
            {
                switch(j)
                {
                    case 0: coefficient(i, j) = 1; break;
                    case 1: coefficient(i, j) = 0; break;
                    case 2: coefficient(i, j) = pointA(i + 1, 0); break;
                    case 3: coefficient(i, j) = pointA(i, 0); break;
                    default: break;
                }
            }
            else
            {
                switch(j)
                {
                    case 0: coefficient(i, j) = 0; break;
                    case 1: coefficient(i, j) = 1; break;
                    case 2: coefficient(i, j) = -pointA(i - 1, 0); break;
                    case 3: coefficient(i, j) = pointA(i, 0); break;
                    default: break;
                }
            }
            
        }
    }
    return coefficient;
}
Vector4d coortrans_4p::SetEstiParam()
{
    Matrix<double, 4, 4> Coefficient;
    Vector4d errorBA;
    Vector4d EstiParam;
    int count = 0;
    int j = 0;
    MatrixN1d pointA(4, 1), pointB(4, 1);
    
    for(int i = 0; i < 6; ++i)
    {
        if(i % 3 == 2)
            continue;
        pointA(j, 0) = this->pointA(i, 0);
        pointB(j, 0) = this->pointB(i, 0);
        j++;
    }
    errorBA = pointB - pointA;
    Coefficient = this->SetCoefficient(pointA);
    EstiParam = Coefficient.inverse() * errorBA;
    return EstiParam;
}

MatrixN4d coortrans_4p::SetB(const int n)
{
    MatrixN1d pointA(n * 2, 1);
    MatrixN4d B;
    B.resize(2 * n, 4);
    int count = 0, j = 0;
    for(int i = 0; i < 12; ++i)
    {
        if(i % 3 == 2)
            continue;
        pointA(j, 0) = this->pointA(i, 0);
        j++;
    }

    B = this->SetCoefficient(pointA);
    return B;
}

MatrixNNd coortrans_4p::SetP(const int n)
{
    Matrix<double, Dynamic, Dynamic> P(2*n, 2*n);
    const int t = 2 * n;
    for(int i = 0; i < 2 * n; ++i)
    {
        for(int j = 0; j < 2 *n; ++j)
        {
            if(i == j)
                P(i, j) = 1;
            else
                P(i, j) = 0;
        }
    }
    return P;
}

MatrixN1d coortrans_4p::SetL(const int n)
{
    
    MatrixN1d pointA;
    MatrixN1d pointB;
    pointA.resize(2 * n, 1);
    pointB.resize(2 * n, 1);
    int j = 0;
    for(int i = 0; i < n * 3; i++)
    {
        if(i % 3 == 2)
            continue;

        pointA(j, 0) = this->pointA(i, 0);
        pointB(j, 0) = this->pointB(i, 0);
        j++;
    }
    MatrixN1d L;
    L.resize(2*n,1);
    
    L = pointB - this->B * this->EstiParam - pointA;
    return L;
}

Vector4d coortrans_4p::Setx0()
{
    MatrixNNd NBB = this->B.transpose() * this->P * this->B;
    MatrixNNd W = this->B.transpose() * this->P * this->L;
    Vector4d x0 = NBB.inverse() * W;

    return x0;
}

MatrixN1d coortrans_4p::SetV()
{
    MatrixN1d V = this->B * this->x0 - this->L;
    return V;
}

//TODO:
MatrixN1d coortrans_4p::SetResidual(int n)
{
    MatrixN1d pointA((this->pointA.rows() + 1) / 3 * 2);
    MatrixN1d pointB((this->pointA.rows() + 1) / 3 * 2);
    int j = 0, i = 0;
    for(; i < this->pointA.rows(); i++)
    {
        if(i % 3 == 2)
            continue;

        pointA(j, 0) = this->pointA(i, 0);
        pointB(j, 0) = this->pointB(i, 0);

        j++;
    }
    MatrixN4d B1(i - 2 * n, 4);
    B1 = this->SetCoefficient(pointA);
    MatrixN1d residual1 = pointA + B1 * (this->EstiParam + this->x0) - pointB;

    return residual1;
}

void coortrans_4p::WriteToFile(string savepath)
{
    fstream out(savepath, ios::out);
    if(!out)
    {
        cerr << "无法写入文件！" << endl;
        exit(1);
    }
    int lines = this->pointA.rows() / 2;
    out << "*********************************4参数转换结果*************************************" << endl << endl;
    out << "必要观测t\t\t\t\t" << 4 << endl;
    out << "总 观 测n\t\t\t\t" << n * 2 << endl;
    out << "多余观测r\t\t\t\t" << n *2 - 4 << endl;
    out << endl;
    MatrixNNd v = this->V.transpose() * this->P * this->V /(lines - 4);
    cout << v << endl;
    double vv = pow(v(0, 0), 0.5);
    out << "单位权中误差\t\t\t\t" << vv << endl;
    out << endl;
    out << "****************转换参数*****************" << endl;
    MatrixNNd Param = this->EstiParam + this->x0;
    out << "DX\t\t\t\t" << Param(0, 0) << endl;
    out << "DY\t\t\t\t" << Param(1, 0) << endl;
    out << "DZ\t\t\t\t" << Param(2, 0) << endl;
    out << "DK\t\t\t\t" << Param(3, 0) << endl;
    out << endl;


    out << "**************************非公共点残差*************************" << endl;
    for(int i=0, k = 0; i < this->residual.rows(); ++i, ++k)
    {
        out << k << "\t\t\t\t" << residual(i, 0) << "\t\t\t\t" << residual(i + 1, 0) << endl;
        i++;
    }
    out << endl;
    out << "********** 115165.237, 118165,826, 0.000 转换结果 ***********" << endl;
    MatrixNNd point(2, 1);
    point << 115165.237, 
    118165.826;
    MatrixNNd result = this->SetCoefficient(point) * Param + point;
    out << setprecision(15) << result(0, 0) << "\t" <<setiosflags(ios::right)<< setw(10) << result(1, 0);
    out.close();
}