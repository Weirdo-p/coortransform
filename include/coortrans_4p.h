#include "common.h"
#include "utils.h"
#include <iomanip>
/**
 * uses: transform points from coordinate A to coordinate B by using Indirect Adjustment
 * 
 * param:
 * dataA  points in coordinate A
 * dataB  points in coordinate B
 * EstiParam  estimate params X0
 * x0     corrections of EstiParam
 * VB     corrections of points B
 * P      weight
 * n      number of public points
 * savepath   where to save the results
*/
class coortrans_4p
{
    private:
        vector31d dataA;
        vector31d dataB;
        MatrixN4d B;
        Vector4d EstiParam;
        MatrixNNd P;
        string savepath;
        MatrixN1d L;
        Vector4d x0;
        MatrixN1d VB;
        MatrixN1d residual;
        int n;

    public:
        coortrans_4p(string pathA, string pathB, const int n, string savepath);
        MatrixN4d GetB(const int n);
        MatrixNNd GetP(const int n);
        Vector4d GetEstiParam();
        Vector4d Getx0();
        MatrixN1d GetVB();
        MatrixN1d GetL(int n);
        MatrixN1d GetResidual(int n);
        void WriteToFile(string savepath);
        
        friend vector31d GetData(string path);
        
};