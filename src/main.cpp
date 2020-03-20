#include "../include/utils.h"
#include "../include/coortrans_4p.h"
#include "../include/coortrans_7p_iter.h"
#include "../include/coortrans_13p_iter.h"
#include "../include/XYZ2BLH.h"
#include "../include/XYZ2NEU.h"
#include <iomanip>


int main()
{
   //  四参数
   
   string pathA_4p = "/home/xz/coding/coortransform/data/XYZ_origin_1.xyz";
   string pathB_4p = "/home/xz/coding/coortransform/data/XYZ_target_1.xyz";
   string savepath_4p = "/home/xz/coding/coortransform/data/4paramresult.txt";
   const int n_4p = 4;
   coortrans_4p trans_4p = coortrans_4p(pathA_4p, pathB_4p, n_4p, savepath_4p);
   

  //  七参数

  string pathA_7p = "/home/xz/coding/coortransform/data/XYZ_origin_2.xyz";
  string pathB_7p = "/home/xz/coding/coortransform/data/XYZ_target_2.xyz";
  string savepath_7p = "/home/xz/coding/coortransform/data/7paramresult.txt";
  const int n_7p = 4;
  coortrans_7p_iter trans_7p = coortrans_7p_iter(pathA_7p, pathB_7p, n_7p, savepath_7p);

  string pathA_13p = "/home/xz/coding/coortransform/data/XYZ_origin_3.xyz";
  string pathB_13p = "/home/xz/coding/coortransform/data/XYZ_target_3.xyz";
  string savepath_13p = "/home/xz/coding/coortransform/data/13paramresult.txt";
  const int n_13p = 11;
  coortrans_13p_iter trans_13p = coortrans_13p_iter(pathA_13p, pathB_13p, n_13p, savepath_13p);

  string data_xyz2blh = "/home/xz/coding/coortransform/data/XYZ2BLHNEU.xyz";
  string savepath_xyz2blh = "/home/xz/coding/coortransform/data/BLH.txt";
  
  XYZ2BLH xyz2blh(data_xyz2blh, savepath_xyz2blh);
  MatrixNNd point = xyz2blh.GetPoint();
  MatrixNNd BLH(point.rows(), point.cols());
  int i = 0;
  for(int j = 0; j < point.cols(); ++j)
  {
    double X = point(0, j);
    double Y = point(1, j);
    double Z = point(2, j);
    double B = xyz2blh.SetB(X, Y, Z);
    double L = xyz2blh.SetL(X, Y);
    double H = xyz2blh.SetH(X, Y, Z, B);
    BLH(0, i) = B;
    BLH(1, i) = L;
    BLH(2, i) = H;
    i++;
  }
  xyz2blh.WriteToFile(savepath_xyz2blh, BLH);


 string datapath = "/home/xz/coding/coortransform/data/XYZ2BLHNEU.xyz";
 string savepath = "/home/xz/coding/coortransform/data/NEU.txt";
 //MatrixNNd point(3, 1);
 XYZ2NEU xyz2neu(datapath, savepath);
}