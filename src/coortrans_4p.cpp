#include "../include/coortrans_4p.h"
#include <iomanip>

coortrans_4p::coortrans_4p(string pathA, string pathB, const int n, string savepath)
{
    cout << "开始初始化" << endl;
    this->n = n;
    this->savepath = savepath;
    cout << "读取坐标文件中..." << endl;
    this->dataA = GetData(pathA);
    this->dataB = GetData(pathB);
    cout << "计算参数估计值..." << endl;
    this->EstiParam = this->GetEstiParam();
    cout << "计算B矩阵..." << endl;
    this->B = this->GetB(n);
    cout << "初始化P矩阵..." <<endl;
    this->P = this->GetP(n);
    cout << "初始化L阵..." <<endl;
    this->L = this->GetL(n);
    cout << "计算结果..." <<endl;
    this->x0 = this->Getx0();
    this->VB = this->GetVB();
    this->residual = this->GetResidual(n);
    cout << "四参数转换成功！ 相关数据已经写入    " << savepath << endl;
    this->WriteToFile(this->savepath);
}

Vector4d coortrans_4p::GetEstiParam()
{
    Matrix<double, 4, 4> Coefficient;
    Vector4d pointA;
    Vector4d pointB;
    Vector4d errorBA;
    Vector4d EstiParam;
    int count = 0;
    int j = 0;
    for(auto data : this->dataA)
    {
        if(count == 2)
            break;

        for(int k = 0; k <=1; k++, j++)
        {
            pointA(j, 0) = data(k, 0);
        }
        count++;
    }
    count = 0;
    j = 0;
    for(auto data : this->dataB)
    {
        if(count == 2)
            break;

        for(int k = 0; k <=1; k++, j++)
        {
            pointB(j, 0) = data(k, 0);
        }
        count++;
    }
    

    errorBA = pointB - pointA;
    Coefficient << 1, 0, pointA[1], pointA[0],
                   0, 1, -pointA[0], pointA[1],
                   1, 0, pointA[3], pointA[2],
                   0, 1, -pointA[2], pointA[3];
    
    EstiParam = Coefficient.inverse() * errorBA;
    return EstiParam;
}

MatrixN4d coortrans_4p::GetB(const int n)
{
    MatrixN1d pointA;
    MatrixN1d pointB;
    MatrixN4d B;
    B.resize(2 * n, 4);
    int count = 0, j = 0;
    pointA.resize(2 * n, 1);
    pointB.resize(2 * n, 1);
    
    for(auto data : dataA)
    {
        if(count == n)
            break;

        for(int k = 0; k <= 1; k++, j++)
        {
            pointA(j, 0) = data(k, 0);
        }
        count++;
    }

    for(int i = 0; i < pointA.rows(); ++i)
    {
        for(int j = 0; j < 4; j++)
        {
            if((i+1) % 2 == 1)
            {
                switch(j)
                {
                    case 0: B(i, j) = 1; break;
                    case 1: B(i, j) = 0; break;
                    case 2: B(i, j) = pointA((int(i+1)/2) * 2 + 1, 0); break;
                    case 3: B(i, j) = pointA(int((i+1)/2) * 2, 0); break;
                    default: break;
                }
            }
            else
            {
                switch(j)
                {
                    case 0: B(i, j) = 0; break;
                    case 1: B(i, j) = 1; break;
                    case 2: B(i, j) = -pointA(i - 1, 0); break;
                    case 3: B(i, j) = pointA(i, 0); break;
                    default: break;
                }
            }
            
        }
    }

    return B;
}

MatrixNNd coortrans_4p::GetP(const int n)
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

MatrixN1d coortrans_4p::GetL(const int n)
{
    
    MatrixN1d pointA;
    MatrixN1d pointB;
    pointA.resize(2 * n, 1);
    pointB.resize(2 * n, 1);
    int count = 0, j = 0;
    for(auto data : dataA)
    {
        if(count == n)
            break;

        for(int k = 0; k <= 1; k++, j++)
        {
            pointA(j, 0) = data(k, 0);
        }
        count++;
    }
    j = 0, count = 0;
    for(auto data : dataB)
    {
        if(count == n)
            break;

        for(int k = 0; k <= 1; k++, j++)
        {
            pointB(j, 0) = data(k, 0);
        }
        count++;
    }
    MatrixN1d L;
    L.resize(2*n,1);
    
    L = pointB - this->B * this->EstiParam - pointA;
    return L;
}

Vector4d coortrans_4p::Getx0()
{
    MatrixNNd NBB = this->B.transpose() * this->P * this->B;
    MatrixNNd W = this->B.transpose() * this->P * this->L;
    Vector4d x0 = NBB.inverse() * W;

    return x0;
}

MatrixN1d coortrans_4p::GetVB()
{
    MatrixN1d VB = this->B * this->x0 - this->L;
    return VB;
}

//TODO:
MatrixN1d coortrans_4p::GetResidual(int n)
{
    MatrixN1d pointA;
    MatrixN1d pointB;
    pointA.resize(2 * (this->dataA.size() - n), 1);
    pointB.resize(2 * (this->dataB.size() - n), 1);
    int count = 0, j = 0;
    for(auto data : dataA)
    {
        if(count < n)
        {
            count++;
            continue;
        }
        for(int k = 0; k <= 1; k++, j++)
        {
            pointA(j, 0) = data(k, 0);
        }
        count++;
    }
    j = 0, count = 0;
    for(auto data : dataB)
    {
        if(count <n)
        {
            count++;
            continue;
        }

        for(int k = 0; k <= 1; k++, j++)
        {
            pointB(j, 0) = data(k, 0);
        }
        count++;
    }

    MatrixN4d B1(((dataA.size() - n)*2 ), 4);
    for(int i = 0; i < pointA.rows(); ++i)
    {
        for(int j = 0; j < 4; j++)
        {
            if((i+1) % 2 == 1)
            {
                switch(j)
                {
                    case 0: B1(i, j) = 1; break;
                    case 1: B1(i, j) = 0; break;
                    case 2: B1(i, j) = pointA((int(i+1)/2) * 2 + 1, 0); break;
                    case 3: B1(i, j) = pointA(int((i+1)/2) * 2, 0); break;
                    default: break;
                }
            }
            else
            {
                switch(j)
                {
                    case 0: B1(i, j) = 0; break;
                    case 1: B1(i, j) = 1; break;
                    case 2: B1(i, j) = -pointA(i - 1, 0); break;
                    case 3: B1(i, j) = pointA(i, 0); break;
                    default: break;
                }
            }
            
        }
    }
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
    out << "*********************************4参数转换结果*************************************" << endl << endl;
    out << "必要观测t\t\t\t\t" << 4 << endl;
    out << "总 观 测n\t\t\t\t" << this->n * 2 << endl;
    out << "多余观测r\t\t\t\t" << this->n *2 - 4 << endl;
    out << endl;

    out << "单位权中误差\t\t\t\t" << this->VB.transpose() * this->P * this->VB << endl;
    out << endl;
    out << "****************转换参数*****************" << endl;
    MatrixNNd Param = this->EstiParam + this->x0;
    out << "DX\t\t\t\t" << Param(0, 0) << endl;
    out << "DY\t\t\t\t" << Param(1, 0) << endl;
    out << "DZ\t\t\t\t" << Param(2, 0) << endl;
    out << "DK\t\t\t\t" << Param(3, 0) << endl;
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
        out << i << "\t\t\t\t" << residual(i, 0) << endl;
    }
    out.close();
}