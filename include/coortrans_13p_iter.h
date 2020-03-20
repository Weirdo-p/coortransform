#include "utils.h"
#include <iomanip>
#include "coortrans.h"

#ifndef _COORTRANS_13P_ITER_H_
#define _COORTRANS_13P_ITER_H_

/**
 * Parametric Adjustment With Constraints
 * 
 * actually, i still use Parametric Adjust to simplify computation
 * basic formula Bx - l = V   Cx + Wx = 0  V'PV = min
 * 
 *                      rotation matrix
 *                             ↓                  [DX
 * [XB       [XA      [1 0 0 0 -ZA YA XA           DY
 *  YB   =    YA   +   0 1 0 ZA 0 -XA YA    *      DZ
 *  ZB]       ZA]      0 0 1 -YA XA 0 ZA]          RX                                        DY
 *                                                 RY
 *                                                 RZ
 *                                                 DK]              
*/
class coortrans_13p_iter : public coortrans
{
    private:
        MatrixNNd B;
        MatrixNNd EstiRotation;
        MatrixNNd ResultRotation;
        MatrixN1d EstiTranslation;
        MatrixN1d ResultTranslation;
        MatrixN1d x0;
        MatrixNNd C;
        MatrixN1d Wx;
        MatrixN1d l;
        MatrixN1d W;
        double EstiMiu;

    public:
    
        coortrans_13p_iter(string pathA, string pathB, const int n, string savepath);

        /**
         * to initialize rotation matrix as identity matrix
        */
        MatrixNNd SetEstiRotation();

        /**
         * to initialize translation matrix as zero matrix
        */
        MatrixN1d SetEstiTranslation();

        // same as param 4/7
        MatrixNNd SetB(const int n);

        // same as param 4/7
        MatrixNNd SetP(const int n) override;

        // same as param 4/7
        MatrixN1d Setx0();

        // same as param 4/7
        MatrixN1d SetV() override;

        // same as param 4/7
        MatrixN1d SetResidual(int n) override;

        /**
         * set matrix C in Cx + Wx = 0
         * *************************************
         * n    the public points 
        */
        MatrixNNd SetC(const int n);

        /**
         * set matrix Wx in Cx + Wx = 0
        */
        MatrixN1d SetWx();

        // same as param 4/7
        MatrixN1d Setl(const int n);

        // same as param 4/7
        int iter(int n);

        // same as param 4/7
        void WriteToFile(string savepath) override;
        friend MatrixN1d SetData(string path);
};


#endif