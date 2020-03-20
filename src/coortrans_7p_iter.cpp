#include "../include/coortrans_7p_iter.h"
#include <iomanip>
#include <math.h>

coortrans_7p_iter::coortrans_7p_iter(string pathA, string pathB, const int n, string savepath)
{
    cout << "开始初始化" << endl;
    this->n = n;
    this->savepath = savepath;
    cout << "读取坐标文件中..." << endl;
    this->pointA = SetData(pathA);
    this->pointB = SetData(pathB);
    this->EstiParam = this->SetEstiParam();
    cout << "初始化中..." << endl;
    this->B = this->SetB(this->n);
    this->P = this->SetP(n);
    this->L = this->SetL (n, this->EstiParam);
    cout << "开始迭代" << endl;
    this->x0 = this->Setx0(this->L);
    this->V = this->SetV();
    this->residual = this->SetResidual(n);
    this->iter(n);
    this->SetResidual(n);
    this->WriteToFile(savepath);
    cout << "成功写入" << savepath << endl;
}

MatrixN7d coortrans_7p_iter::SetCoefficient(MatrixN1d pointA)
{
    int rows = pointA.rows();
    MatrixN7d coefficient(rows, 7);
    for(int i = 0; i < rows; ++i)
    {
        for(int j = 0; j < 7; ++j)
        {
            if((i+1)%3 == 1)
            {
                switch(j)
                {
                    case 0: coefficient(i, j) = 1; break;
                    case 1: coefficient(i, j) = 0; break;
                    case 2: coefficient(i, j) = 0; break;
                    case 3: coefficient(i, j) = 0; break;
                    case 4: coefficient(i, j) = -pointA(i + 2, 0); break;
                    case 5: coefficient(i, j) = pointA(i + 1, 0); break;
                    case 6: coefficient(i, j) = pointA(i, 0); break;
                    default: break;
                }
            }
            else if((i+1)%3 == 2)
            {
                switch(j)
                {
                    case 0: coefficient(i, j) = 0; break;
                    case 1: coefficient(i, j) = 1; break;
                    case 2: coefficient(i, j) = 0; break;
                    case 3: coefficient(i, j) = pointA(i + 1, 0); break;
                    case 4: coefficient(i, j) = 0; break;
                    case 5: coefficient(i, j) = -pointA(i - 1, 0); break;
                    case 6: coefficient(i, j) = pointA(i, 0); break;
                    default: break;
                }
            }
            else
            {
                switch(j)
                {
                    case 0: coefficient(i, j) = 0; break;
                    case 1: coefficient(i, j) = 0; break;
                    case 2: coefficient(i, j) = 1; break;
                    case 3: coefficient(i, j) = -pointA(i - 1, 0); break;
                    case 4: coefficient(i, j) = pointA(i - 2, 0); break;
                    case 5: coefficient(i, j) = 0; break;
                    case 6: coefficient(i, j) = pointA(i, 0); break;
                    default: break;
                }
            }
            
        }
    }
    return coefficient;
}

Matrix71d coortrans_7p_iter::SetEstiParam()
{
    MatrixN7d coefficient;
    MatrixN1d pointA = this->pointA;
    MatrixN1d pointB = this->pointB;
    coefficient.resize(7, 7);
    int j = 0;
    MatrixNNd tmp = this->SetCoefficient(pointA);
    for(int i = 0; i < coefficient.rows(); ++i)
        for(int k = 0; k < coefficient.cols(); ++k)
            coefficient(i, k) = tmp(i, k);
    Matrix71d A;
    Matrix71d B;

    for(int i=0;i < A.rows(); i++)
    {
        A(i, 0) = pointA(i, 0);
        B(i, 0) = pointB(i, 0);
    }

    Matrix71d EstiParam = coefficient.inverse() * (B - A);
    return EstiParam;
}

MatrixN7d coortrans_7p_iter::SetB(int n)
{
    MatrixN1d pointA;
    MatrixN7d B;
    pointA.resize(n * 3, 1);
    B.resize(n * 3, 7);
    int j = 0;
    for(int i = 0; i < n * 3; ++i)
    {
        pointA(i, 0) = this->pointA(i, 0);
    }
    B = this->SetCoefficient(pointA);
    return B;
}

MatrixNNd coortrans_7p_iter::SetP(int n)
{
    MatrixNNd P1(3*n, 3*n);

    for(int i = 0; i < 3 * n; ++i)
    {
        for(int j = 0; j < 3 *n; ++j)
        {
            if(i == j)
                P1(i, j) = 1;
            else
                P1(i, j) = 0;
        }
    }
    return P1;
}

MatrixN1d coortrans_7p_iter::SetL (int n, MatrixNNd EstiParam1)
{
    MatrixN1d pointA;
    MatrixN1d pointB;
    pointA.resize(n * 3, 1);
    pointB.resize(n * 3, 1);
    int j = 0;
    for(int i = 0; i < 3 * n; ++i)
    {
        pointA(i, 0) = this->pointA(i, 0);
        pointB(i, 0) = this->pointB(i, 0);
    }
    MatrixN1d L = pointB - pointA - this->B * EstiParam1;
    return L;
}

Matrix71d coortrans_7p_iter::Setx0(MatrixN1d L1)
{
    MatrixNNd NBB = this->B.transpose() * this->P * this->B;
    MatrixNNd W = this->B.transpose() * this->P * L1;
    Matrix71d x0 = NBB.inverse() * W;

    return x0;
}

MatrixN1d coortrans_7p_iter::SetV()
{
    MatrixN1d V = this->B * this->x0 - this->L;
    return V;
}

MatrixN1d coortrans_7p_iter::SetResidual (int n)
{
    MatrixNNd coe = this->SetCoefficient(this->pointA);
    return (this->pointA - this->pointB + coe*(this->x0 + this->EstiParam));
}

int coortrans_7p_iter::iter(int n)
{
    Matrix71d x0 = this->x0;
    double thres = 1e-18;
    MatrixNNd Param = this->EstiParam;
    MatrixNNd X = Param + x0;
    
    int count = 0;
    MatrixN1d L;
    MatrixNNd dx;

    while(true)
    {
        L = this->SetL(n, X);
        dx = this->Setx0(L);
        x0 = dx;
        X += dx;
        double max = dx(0,0);
        for(int i =0; i < dx.rows(); ++i)
            if(dx(i, 0) > max)
                max = dx(i, 0);
        count ++;
        if(max < thres)
            break;
    }
    this->x0 = x0;
    this->EstiParam = X;
    this->L = this->SetL(n, X);
    this->V = this->SetV();
    return count;
}

void coortrans_7p_iter::WriteToFile(string savepath)
{
    fstream out(savepath, ios::out);
    if(!out)
    {
        cerr << "无法写入文件！" << endl;
        exit(1);
    }
    int lines = (pointA.rows()+1) / 3;
    out << "*********************************7参数转换结果*************************************" << endl << endl;
    out << "必要观测t\t\t\t\t" << 7 << endl;
    out << "总 观 测n\t\t\t\t" << n * 3 << endl;
    out << "多余观测r\t\t\t\t" << n *3 - 7 << endl;
    out << endl;
    MatrixNNd v = this->V.transpose() * this->P * this->V /(lines*3 - 7);
    double vv = pow(v(0, 0), 0.5);
    out << "单位权中误差\t\t\t\t" << vv << endl;
    out << endl;
    out << "****************转换参数*****************" << endl;
    MatrixNNd Param = this->EstiParam + this->x0;
    out << "DX\t\t\t\t" << Param(0, 0) << endl;
    out << "DY\t\t\t\t" << Param(1, 0) << endl;
    out << "DZ\t\t\t\t" << Param(2, 0) << endl;
    out << "RX\t\t\t\t" << Param(3, 0) << endl;
    out << "RY\t\t\t\t" << Param(4, 0) << endl;
    out << "RZ\t\t\t\t" << Param(5 ,0) << endl;
    out << "DK\t\t\t\t" << Param(6, 0) << endl;
    out << endl;

    out << "******************各点位中误差********************" << endl;
    out << "点号\t\t\t\t" << "中误差" << endl;
    MatrixNNd NBB = this->B.transpose() * this->P * this->B;

    MatrixNNd midErrors = this->B * NBB.inverse() * B.transpose();
    for(int i = 0; i < midErrors.rows(); ++i)
        for(int j = 0; j < midErrors.cols(); ++j)
            if((i+1) % 2 == 0)
                if(i==j)
                    out << int((i+1)/2) << "\t\t\t\t" << sqrt(midErrors(i-1, j-1)+midErrors(i, j)) << endl;
    out << endl;

    out << "**************************非公共点残差*************************" << endl;
    for(int i=0; i < this->residual.rows(); ++i)
    {
        if(i < 3 * n)
            continue;
        out << i - 3 * n << "\t\t\t\t" << residual(i, 0) << endl;
    }
    out.close();
}
