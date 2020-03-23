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
   cout << "*********************************************************" << endl;
   cout << "*****欢迎使用坐标转换系统v1.0*****" << endl;
   string pathA;
   string pathB;
   string savepath;
   int num;
   int choice;

   while(true)
   {
      cout << "输入功能前相应序号以执行相应功能" << endl;
      cout << "1.四参数转换" << endl;
      cout << "2.七参数转换" << endl;
      cout << "3.十三参数转换" << endl;
      cout << "4.地心地固转大地坐标" << endl;
      cout << "5.地心地固转站心坐标" << endl;
      cout << "6.退出" << endl;
      cout << "************************************************************" << endl;

      cin >> choice;
      if(choice == 1)
      {
        string pathA_4p = "/home/xz/coding/coortransform/data/XYZ_origin_1.xyz";
        string pathB_4p = "/home/xz/coding/coortransform/data/XYZ_target_1.xyz";
        string savepath_4p = "/home/xz/coding/coortransform/data/4paramresult.txt";
        int n_4p = 4;

        cout << "请输入待转换坐标文件路径(输入m则使用默认路径)" << endl;
        cin >> pathA;
        cout << "请输入目标坐标系文件路径(输入m则使用默认路径)" << endl;
        cin >> pathB;
        cout << "请输入公共点个数" << endl;
        cin >> num;
        if(pathA != "m")
          pathA_4p = pathA;
        if(pathB != "m")
          pathB_4p = pathB;
        if(num < 4)
          cout << "公共点数量不足以计算！" << endl;
        else
          n_4p = num;
        cout << "请输入保存路径(输入m则使用默认保存路径)" << endl;
        cin >> savepath;
        if(savepath != "m")
          savepath_4p = savepath;
        coortrans_4p trans_4p = coortrans_4p(pathA_4p, pathB_4p, n_4p, savepath_4p);
      }   

      //  七参数
      if(choice == 2)
      {
        string pathA_7p = "/home/xz/coding/coortransform/data/XYZ_origin_2.xyz";
        string pathB_7p = "/home/xz/coding/coortransform/data/XYZ_target_2.xyz";
        string savepath_7p = "/home/xz/coding/coortransform/data/7paramresult.txt";
        
        int n_7p = 4;
        
        cout << "请输入待转换坐标文件路径(输入m则使用默认路径)" << endl;
        cin >> pathA;
        cout << "请输入目标坐标系文件路径(输入m则使用默认路径)" << endl;
        cin >> pathB;
        cout << "请输入公共点个数" << endl;
        cin >> num;
        cout << "请输入保存路径(输入m则使用默认保存路径)" << endl;
        cin >> savepath;
        if(savepath != "m")
          savepath_7p = savepath;

        if(pathA != "m")
          pathA_7p = pathA;
        if(pathB != "m")
          pathB_7p = pathB;
        if(num < 3)
          cout << "公共点数量不足以计算！" << endl;
        else
          n_7p = num;

        coortrans_7p_iter trans_7p = coortrans_7p_iter(pathA_7p, pathB_7p, n_7p, savepath_7p);

      }
      if(choice == 3)
      {
        string pathA_13p = "/home/xz/coding/coortransform/data/XYZ_origin_3.xyz";
        string pathB_13p = "/home/xz/coding/coortransform/data/XYZ_target_3.xyz";
        string savepath_13p = "/home/xz/coding/coortransform/data/13paramresult.txt";
        int n_13p = 11;
        cout << "请输入待转换坐标文件路径(输入m则使用默认路径)" << endl;
        cin >> pathA;
        cout << "请输入目标坐标系文件路径(输入m则使用默认路径)" << endl;
        cin >> pathB;
        cout << "请输入公共点个数" << endl;
        cin >> num;
        cout << "请输入保存路径(输入m则使用默认保存路径)" << endl;
        cin >> savepath;
        if(savepath != "m")
          savepath_13p = savepath;


        if(pathA != "m")
          pathA_13p = pathA;
        if(pathB != "m")
          pathB_13p = pathB;
        if(num < 3)
          cout << "公共点数量不足以计算！" << endl;
        else
          n_13p = num;
        coortrans_13p_iter trans_13p = coortrans_13p_iter(pathA_13p, pathB_13p, n_13p, savepath_13p);

      }
      if(choice == 4)
      {
        string data_xyz2blh = "/home/xz/coding/coortransform/data/XYZ2BLHNEU.xyz";
        string savepath_xyz2blh = "/home/xz/coding/coortransform/data/BLH.txt";
        cout << "请输入待转换坐标文件路径(输入m则使用默认路径)" << endl;
        cin >> pathA;
        cout << "请输入保存路径(输入m则使用默认保存路径)" << endl;
        cin >> savepath;
        if(savepath != "m")
          savepath_xyz2blh = savepath;
        if(pathA != "m")
          data_xyz2blh = pathA;

        cout << "开始转换BLH" << endl;
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
        cout << "转换成功，转换结果写入   " << savepath_xyz2blh << endl;

      }
      if(choice == 5)
      {
        string datapath_xyz2neu = "/home/xz/coding/coortransform/data/XYZ2BLHNEU.xyz";
        string savepath_xyz2neu = "/home/xz/coding/coortransform/data/NEU.txt";
        cout << "请输入待转换坐标文件路径(输入m则使用默认路径)" << endl;
        cin >> pathA;
        cout << "请输入保存路径(输入m则使用默认保存路径)" << endl;
        cin >> savepath;
        if(savepath != "m")
          savepath_xyz2neu = savepath;
        if(pathA != "m")
          datapath_xyz2neu = pathA;


        XYZ2NEU xyz2neu(datapath_xyz2neu, savepath_xyz2neu);
        cout << "转换成功，转换结果写入    " << savepath_xyz2neu << endl;

      }
      if(choice == 6)
      {
        break;
      }
      else
      {
        cout << "请输入有效数字" << endl;
      }

   }
 return 0;
}