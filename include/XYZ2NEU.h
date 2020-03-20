#include "common.h"
#include "utils.h"
#include "XYZ2BLH.h"

#ifndef _XYZ2NEU_H_
#define _XYZ2NEU_H_

/**
 * to transfer from xyz to neu
 *                 rotation
 *                     ↓
 * [n      [-sinBcosL -sinBsinL cosB        [Xq - Xp
 *  e   =   -sinL      cosL     0       *    Yq - Yp
 *  u]      cosBcosL   cosBsinL  sinB]       Zq - Yp]
 * P is the origin point of neu, expressed in xyz
 * where B(latitude), L(longitude) means P in BLH
 * q is the point we want to transform
 * 
 * 
 * *******************************************
 * point    point we want to transform into neu
 * Cxyz     point P in xyz
 * Cblh     point P in blh
 * NEU      the results
 * helper   transform P(xzy) to P(blh)
*/
class XYZ2NEU
{
    private:
        string savepath;
        MatrixNNd point;
        MatrixNNd Cxyz;
        MatrixNNd Cblh;
        MatrixNNd rotation;
        MatrixNNd NEU;
        XYZ2BLH helper;
    public:
        XYZ2NEU(string datapath, string savepath, MatrixNNd Cpoint);
        XYZ2NEU(string datapath, string savepath);
        MatrixNNd SetRotation();
        MatrixNNd Compute();
        void WriteToFile(string savepath);
};

#endif