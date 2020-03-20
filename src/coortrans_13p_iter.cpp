#include "../include/coortrans_13p_iter.h"
#include <iomanip>
#include <math.h>

coortrans_13p_iter::coortrans_13p_iter(string pathA, string pathB, const int n, string savepath)
{
    cout << "读取文件中" << endl;
    this->pointA = SetData(pathA);
    this->n = n;
    this->pointB = SetData(pathB);
    this->savepath = savepath;
    cout << "进行初始化" << endl;
    this->EstiMiu = 1;
    this->EstiRotation = this->SetEstiRotation();
    this->EstiTranslation = this->SetEstiTranslation();
    this->P = this->SetP(n);
    this->B = this->SetB(n);
    this->C = this->SetC(n);
    this->Wx = this->SetWx();
    this->l = this->Setl(n);
    this->x0 = this->Setx0();
    cout << "开始迭代" << endl;
    this->iter(n);

    this->V = this->SetV();
    this->residual = this->SetResidual(n);
    this->WriteToFile(savepath);
    cout << "13参数转化完毕，成功写入文件" << endl;
}

MatrixNNd coortrans_13p_iter::SetEstiRotation()
{
    MatrixNNd rotation(3, 3);
    for(int i = 0; i < rotation.rows(); ++i)
        for(int j = 0; j < rotation.cols(); ++j)
            if(i == j)
                rotation(i, j) = 1;
            else
                rotation(i, j) = 0;
    return rotation;
}

MatrixN1d coortrans_13p_iter::SetEstiTranslation()
{
    MatrixNNd trans(3, 1);
    for(int i = 0; i < trans.rows(); ++i)
        trans(i, 0) = 0;
    return trans;
}

MatrixNNd coortrans_13p_iter::SetP(int n)
{
    MatrixNNd P(n * 3 + 6, n * 3 + 6);
    for(int i = 0; i < P.rows(); ++i)
        for(int j = 0; j < P.cols(); ++j)
            if(i == j)
                P(i, j) = 1;
            else
                P(i, j) = 0;

    return P;
}

MatrixNNd coortrans_13p_iter::SetB(const int n)
{
    MatrixNNd B(n * 3, 13);
    MatrixNNd point = this->pointA;
    double miu = this->EstiMiu;
    for(int i = 0; i < B.rows(); ++i)
    {
        for(int j = 0; j < B.cols(); ++j)
        {
            if((i + 1) % 3 == 1)
            {
                switch(j)
                {
                    case 0: B(i, j) = 1; break;
                    case 3: B(i, j) = EstiRotation(0, 0) * point(i, 0) + EstiRotation(0, 1) * point(i + 1, 0) + EstiRotation(0, 2) * point(i + 2); break; 
                    case 4: B(i, j) = miu * point(i, 0); break;
                    case 5: B(i, j) = miu * point(i + 1, 0); break;
                    case 6: B(i, j) = miu * point(i + 2, 0); break;
                    default: B(i, j) = 0; break;
                }
            }

            else if((i + 1) % 3 == 2)
            {
                switch(j)
                {
                    case 1: B(i, j) = 1; break;
                    case 3: B(i, j) = EstiRotation(1, 0) * point(i - 1, 0) + EstiRotation(1, 1) * point(i, 0) + EstiRotation(1, 2) * point(i + 1); break;
                    case 7: B(i, j) = miu * point(i - 1); break;
                    case 8: B(i, j) = miu * point(i, 0); break;
                    case 9: B(i, j) = miu * point(i + 1, 0); break;
                    default: B(i, j) = 0; break;
                }
            }
            else
            {
                switch(j)
                {
                    case 2: B(i, j) = 1; break;
                    case 3: B(i, j) = EstiRotation(2, 0) * point(i - 2, 0) + EstiRotation(2, 1) * point(i - 1, 0) + EstiRotation(2, 2) * point(i); break;
                    case 10: B(i, j) = miu * point(i - 2, 0); break;
                    case 11: B(i, j) = miu * point(i - 1, 0); break;
                    case 12: B(i, j) = miu * point(i, 0); break;
                    default: B(i, j) = 0; break;
                }
            }
            
        }
    }
    return B;
}

MatrixNNd coortrans_13p_iter::SetC(const int n)
{
    MatrixNNd C(6, 13);
    MatrixNNd rotation = this->EstiRotation;
    for(int i = 0; i < C.rows(); ++i)
    {
        for(int j = 0; j < C.cols(); ++j)
        {
            if(i == 0)
            {
                switch(j)
                {
                    case 4: C(i, j) = 2 * rotation(0, 0); break;
                    case 5: C(i, j) = 2 * rotation(0, 1); break;
                    case 6: C(i, j) = 2* rotation(0, 2); break;
                    default: C(i, j) = 0; break;
                }
            }
            else if(i == 1)
            {
                switch(j)
                {
                    case 7: C(i, j) = 2 * rotation(1, 0); break;
                    case 8: C(i, j) = 2 * rotation(1, 1); break;
                    case 9: C(i, j) = 2 * rotation(1, 2); break;
                    default: C(i, j) = 0; break;
                }
            }

            else if(i == 2)
            {
                switch(j)
                {
                    case 10: C(i, j) = 2 * rotation(2, 0); break;
                    case 11: C(i, j) = 2 * rotation(2, 1); break;
                    case 12: C(i, j) = 2 * rotation(2, 2); break;
                    default: C(i, j) = 0; break;
                }
            }
            else if(i == 3)
            {
                switch(j)
                {
                    case 4: C(i, j) = rotation(0, 1); break;
                    case 5: C(i, j) = rotation(0, 0); break;
                    case 7: C(i, j) = rotation(1, 1); break;
                    case 8: C(i, j) = rotation(1, 0); break;
                    case 10: C(i, j) = rotation(2, 1); break;
                    case 11: C(i, j) = rotation(2, 0); break;
                    default: C(i, j) = 0; break;
                }
            }
            else if(i == 4)
            {
                switch(j)
                {
                    case 4: C(i, j) = rotation(0, 2); break;
                    case 6: C(i, j) = rotation(0, 0); break;
                    case 7: C(i, j) = rotation(1, 2); break;
                    case 9: C(i, j) = rotation(1, 0); break;
                    case 10: C(i, j) = rotation(2, 2); break;
                    case 12: C(i, j) = rotation(2, 0); break;
                    default: C(i, j) = 0; break;
                }
            }
            else
            {
                switch(j)
                {
                    case 5: C(i, j) = rotation(0, 2); break;
                    case 6: C(i, j) = rotation(0 ,1); break;
                    case 8: C(i, j) = rotation(1, 2); break;
                    case 9: C(i, j) = rotation(1, 1); break;
                    case 11: C(i, j) = rotation(2 ,2); break;
                    case 12: C(i, j) = rotation(2, 1); break;
                    default: C(i, j) = 0; break;
                }
            }
        }
    }
    return C;
}

MatrixN1d coortrans_13p_iter::SetWx()
{
    MatrixNNd Wx(6, 1);
    MatrixNNd rotation = EstiRotation;
    double a = pow(rotation(0, 0), 2) + pow(rotation(0, 1), 2) + pow(rotation(0, 2), 2) - 1;
    double b = pow(rotation(1, 0), 2) + pow(rotation(1, 1), 2) + pow(rotation(1, 2), 2) - 1;
    double c = pow(rotation(2, 0), 2) + pow(rotation(2, 1), 2) + pow(rotation(2, 2), 2) - 1;
    double a12 = rotation(0, 0) * rotation(0, 1) + rotation(1, 0) * rotation(1, 1) + rotation(2, 0) * rotation(2, 1);
    double a13 = rotation(0, 0) * rotation(0, 2) + rotation(1, 0) * rotation(1, 2) + rotation(2, 0) * rotation(2, 2);
    double a23 = rotation(0, 1) * rotation(0, 2) + rotation(1, 1) * rotation(1, 2) + rotation(2, 1) * rotation(2, 2);
    Wx(0, 0) = a;
    Wx(1, 0) = b;
    Wx(2, 0) = c;
    Wx(3, 0) = a12;
    Wx(4, 0) = a13;
    Wx(5, 0) = a23;
    return -Wx;
}
MatrixN1d coortrans_13p_iter::Setl(const int n)
{
    MatrixNNd rotation = this->EstiRotation;
    double miu = this->EstiMiu;
    MatrixNNd pointA(n * 3, 1), pointB(n * 3, 1);
    MatrixNNd trans(3 , n);
    for(int i = 0; i < pointA.rows(); ++i)
    {
        pointA(i, 0) = this->pointA(i, 0);
        pointB(i, 0) = this->pointB(i, 0);
    }
    for(int i = 0; i < trans.cols(); ++i)
        for(int j = 0; j < trans.rows(); ++j)
            trans(j, i) = this->EstiTranslation(j, 0);
    pointA.resize(3, n);
    pointB.resize(3, n);
    MatrixNNd l = -(miu * rotation * pointA - pointB + trans);
    l.resize(3 * n, 1);
    return l;
}


MatrixN1d coortrans_13p_iter::Setx0()
{
    int n = this->n;
    MatrixNNd B_iter(3 * n + 6, 13);
    MatrixNNd l_iter(3 * n + 6, 1);
    B_iter.block(0, 0, this->B.rows(), this->B.cols()) = B;
    B_iter.block(this->B.rows(), 0, this->C.rows(), C.cols()) = C;
    l_iter.block(0, 0, this->l.rows(), l.cols()) = l;
    l_iter.block(l.rows(), 0, this->Wx.rows(), this->Wx.cols()) = Wx;
    MatrixNNd NBB = B_iter.transpose() * P * B_iter;
    MatrixNNd x0 = NBB.inverse() * B_iter.transpose() * P * l_iter;
    return x0;
}

int coortrans_13p_iter::iter(int n)
{
    MatrixNNd &rotation = this->EstiRotation;
    MatrixN1d &trans = this->EstiTranslation;
    double &miu = this->EstiMiu;
    MatrixN1d &x0 = this->x0;
    double thres = 1e-10;
    for(int i = 0; i < trans.rows(); ++i)
        trans(i, 0) += x0(i, 0);
    int k = 4;
    for(int i = 0; i < rotation.rows(); ++i)
        for(int j = 0; j < rotation.cols(); ++j)
            rotation(i, j) += x0(k++, 0);
    miu += x0(3, 0);
    int count = 0;
    while(true)
    {
        this->B = this->SetB(n);
        this->C = this->SetC(n);
        this->Wx = this->SetWx();
        this->l = this->Setl(n);
        x0 = this->Setx0();
        for(int i = 0; i < trans.rows(); ++i)
            trans(i, 0) += x0(i, 0);
        k = 4;

        for(int i = 0; i < rotation.rows(); ++i)
            for(int j = 0; j < rotation.cols(); ++j)
                rotation(i, j) += x0(k++, 0);
        miu += x0(3, 0);
        double max = x0(0, 0);
        count++;
        for(int i = 0; i < x0.rows(); ++i)
            if(x0(i, 0) > max)
                max = x0(i, 0);
        if(max < thres)
            break;
    }
    this->B = this->SetB(n);
    this->C = this->SetC(n);
    return count;
}

MatrixN1d coortrans_13p_iter::SetV()
{
    MatrixNNd V = B * x0 - l;
    return V;
}

MatrixN1d coortrans_13p_iter::SetResidual(int n)
{
    MatrixNNd rotation = this->EstiRotation;
    MatrixNNd pointA = this->pointA;
    MatrixNNd pointB = this->pointB;
    pointA.resize(3, this->pointA.rows()/3);
    pointB.resize(3, this->pointB.rows()/3);
    double miu = this->EstiMiu;
    MatrixNNd trans(3 , this->pointA.rows()/3);
    for(int i = 0; i < trans.cols(); ++i)
        for(int j = 0; j < trans.rows(); ++j)
            trans(j, i) = this->EstiTranslation(j, 0);
    MatrixNNd l = -(miu * rotation * pointA - pointB + trans);
    l.resize(this->pointA.rows(), 1);
    return l;
    
}

void coortrans_13p_iter::WriteToFile(string savepath)
{
    fstream out(savepath, ios::out);
    if(!out)
    {
        cerr << "无法写入文件！" << endl;
        exit(1);
    }
    int lines = (pointA.rows()+1) / 3;
    out << "*********************************13参数转换结果，单位mm*************************************" << endl << endl;
    out << "必要观测t\t\t\t\t" << 7 << endl;
    out << "总 观 测n\t\t\t\t" << n * 3 << endl;
    out << "多余观测r\t\t\t\t" << n * 3 - 7 << endl;
    out << endl;
    MatrixNNd v = this->V.transpose() * this->V /(lines*3 - 7);
    double vv = pow(v(0, 0), 0.5);
    out << "单位权中误差\t\t\t\t" << vv << endl;
    out << endl;
    out << "****************转换参数*****************" << endl;
    out << "旋转矩阵：" << endl;
    MatrixNNd rotation = this->EstiRotation;
    for(int i = 0; i < rotation.rows(); ++i)
    {
        for(int j = 0; j < rotation.cols(); ++j)
        {
            out << setprecision(15) << rotation(i, j) << "\t\t";
        }
        out << endl;
    }
    out << endl;
    
    double miu = this->EstiMiu;
    MatrixNNd trans = this->EstiTranslation;
    out << "尺度因子：" << endl;
    out << miu << endl;
    out << endl;
    out << "平移矩阵" << endl;
    for(int i = 0; i < trans.rows(); ++i)
        out << trans(i, 0) << endl;
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