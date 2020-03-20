#include <iomanip>

#ifndef _COORTRANS_4P_ITER_H_
#define _COORTRANS_4P_ITER_H_
#include "common.h"
#include "utils.h"
#include "coortrans.h"

/**
 * uses: transform points from coordinate A to coordinate B by using Indirect Adjustment
 * basic formula   v = Bx - l, V'PV = min
 * [XB      [XA      [1 0 YA XA         [DX
 *       =        +                 *    DY
 *  YB]      YA]      0 1 -XA YA]        DR
 *                                       DK]
 * ***************************************************
 * param:
 * dataA  points in coordinate A
 * dataB  points in coordinate B
 * EstiParam  estimate params X0
 * x0     corrections of EstiParam
 * VB     corrections of points B
 * P      weight
 * n      number of public points
 * savepath   where to save the results
 * 
 * ***********************************
*/
class coortrans_4p : public coortrans
{
    private:
        MatrixN4d B;
        Vector4d EstiParam;
        Vector4d x0;

    public:
        /**
         * construct function     when you initialize the object, it will compute automatically
         * *******************************
         * pathA    the path of origin data
         * pathB    the path of target data
         * n        the number of public points
         * savepath     the path to save results
        */
        coortrans_4p(string pathA, string pathB, const int n, string savepath);

        /**
         * function    to set Matrix B in Bx - l = v
         * 
         * ********************************
         * n      number of public points
         * ********************************
         */
        MatrixN4d SetB(const int n);

        /**
         * function    to set weight matrix P, it is derived from base class
         * ******************************************************
         * n     number of public points
         * *********************************
        */
        MatrixNNd SetP(const int n) override;

        /**
         * function    to compute estiparams(DX, DY, DR, DK) according to orgin data
         * 
        */
        Vector4d SetEstiParam();

        /**
         * to compute the correction of the estiparam
        */
        Vector4d Setx0();

        /**
         * to compute the correction of target point to calculate V'PV
        */
        MatrixN1d SetV() override;

        /**
         * to calculate matrix l of Bx - l = v
         * *****************************************
         * n     the number of public points
        */
        MatrixN1d SetL(int n);

        /**
         * using the point transformed by using params we calculated minus the target point
        */
        MatrixN1d SetResidual(int n) override;

        /**
         * to help get matrix B
         */
        MatrixN4d SetCoefficient(MatrixN1d pointA);

        /**
         * write result into file
         * ********************************
         * savepath   path of result file
        */
        void WriteToFile(string savepath) override;
        

        /**
         * to get data from file
         * ***************************
         * path    path of data file
        */
        friend MatrixN1d SetData(string path);
        
};

#endif // !1_COORTRANS_4P_ITER_H
