#include "../include/utils.h"
#include "../include/coortrans_4p.h"
#include <iomanip>


int main()
{
   string pathA = "/home/xz/coding/coortransform/data/XYZ_origin_1.xyz";
   string pathB = "/home/xz/coding/coortransform/data/XYZ_target_1.xyz";
   string savepath = "/home/xz/coding/coortransform/data/4paramresult.txt";
   const int n_4p = 4;
   coortrans_4p trans_4p = coortrans_4p(pathA, pathB, n_4p, savepath);
}