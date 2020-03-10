#include "../include/utils.h"
#include <iomanip>


int main()
{
    
    string path = "/home/xz/coding/coortransform/data/XYZ_origin_1.xyz";
    vector31d data;
    data = GetData(path);
    for (auto i : data)
    {
        cout << setprecision(15) << i.transpose() << endl ;

    }
    cout << endl;
    cout << setprecision(15) << data[0] * data[1].transpose() << endl;
    return 0;
}