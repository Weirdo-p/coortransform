#ifndef _COORTRANS_7P_ITER_H_
#define _COORTRANS_7P_ITER_H_

#include "utils.h"
#include <iomanip>
#include "coortrans.h"



/**
 * last edited on 13th March
 * fuction  to run the 7 params coordinates transform
 * 
 * basic formula   v = Bx - l, V'PV = min
 * [XB      [XA      [1 0 YA XA         [DX
 *       =        +                 *    DY
 *  YB]      YA]      0 1 -XA YA]        DR
 *                                       DK]
 * ***************************************************
 * params are same to param 4 transforming
 * here we use iteration to update our params, to be honest, it's unnecessary due to tiny angle
 */
class coortrans_7p_iter : public coortrans
{
    private:
        MatrixN7d B;
        Matrix71d EstiParam;
        Matrix71d x0;

    public:
        //functions are same to 4 params transforming
        coortrans_7p_iter(string pathA, string pathB, const int n, string savepath);
        MatrixN7d SetB(const int n);
        MatrixNNd SetP(const int n) override;
        Matrix71d SetEstiParam();
        Matrix71d Setx0(MatrixN1d L);
        MatrixN1d SetV() override;
        MatrixN1d SetL(int n, MatrixNNd Estiparam);
        MatrixN1d SetResidual(int n) override;
        MatrixN7d SetCoefficient(MatrixN1d pointA);

        /**
         * use iteration to update our params
         * here we only to update params not the observation
         * 
         * **************************************************
         * n    number of public points
        */
        int iter(int n);

        void WriteToFile(string savepath) override;
        friend MatrixN1d SetData(string path);
};






#endif