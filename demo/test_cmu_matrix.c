#include <stdio.h>
#include "cmu_static_matrix.h"


void testMatrix2D() 
{
    Matrix2D A = {{{1, 1}, {1, 1}}};
    Matrix2D B = {{{2, 2}, {2, 2}}};
    Vector2D v = {{1, 1}};
    Matrix2D result_mat = { 0.0 };
    Vector2D result_vec = { 0.0 };
    double reuslt_scalar = 0.0;

    printf("testMatrix2D()\n");

    printf("A = [[%f, %f], [%f, %f]]\n", 
            A.data[0][0], A.data[0][1], A.data[1][0], A.data[1][1]);
    printf("B = [[%f, %f], [%f, %f]]\n", 
            B.data[0][0], B.data[0][1], B.data[1][0], B.data[1][1]);
    printf("v = [%f, %f]\n", v.x, v.y);

    // 测试矩阵加法
    result_mat = AddMat2D(&A, &B);
    printf("AddMat2D: result = [[%f, %f], [%f, %f]]\n", 
            result_mat.data[0][0], result_mat.data[0][1], 
            result_mat.data[1][0], result_mat.data[1][1]);

    // 测试矩阵乘法
    result_mat = DotMat2D(&A, &B);
    printf("MatMul2D: result = [[%f, %f], [%f, %f]]\n", 
            result_mat.data[0][0], result_mat.data[0][1], 
            result_mat.data[1][0], result_mat.data[1][1]);

    // 测试矩阵与标量相乘
    result_mat = DotScalarMat2D(2.0, &A);
    printf("ScalarMulMat2D: result = [[%f, %f], [%f, %f]]\n", 
            result_mat.data[0][0], result_mat.data[0][1], 
            result_mat.data[1][0], result_mat.data[1][1]);
    // 测试矩阵乘以向量
    result_vec = DotMatVec2D(&A, &v);
    printf("MatMulVec2D: result_vec = [%f, %f]\n", result_vec.x, result_vec.y);

    // 测试向量数乘
    result_vec = DotScalarVec2D(3.0, &v);
    printf("DotScalarVec2D: result_vec = [%f, %f]\n", result_vec.x, result_vec.y);

    // 测试向量点乘
    reuslt_scalar = DotVec2D(&v, &v);
    printf("DotVec2D: result_scalar = %f\n", reuslt_scalar);

    // 测试向量转对角矩阵
    result_mat = DiagVecToMat2D(&v);
    printf("DiagVecToMat2D: result_mat = [[%f, %f], [%f, %f]]\n", 
            result_mat.data[0][0], result_mat.data[0][1], 
            result_mat.data[1][0], result_mat.data[1][1]);

    // 测试向量平方
    result_vec = SquareVec2D(&v);
    printf("SquareVec2D: result_vec = [%f, %f]\n", result_vec.x, result_vec.y);

    // 测试向量开方
    result_vec = SqrtVec2D(&v);
    printf("SqrtVec2D: result_vec = [%f, %f]\n", result_vec.x, result_vec.y);

    // 测试向量模长
    reuslt_scalar = NormVec2D(&v);
    printf("NormVec2D: result_scalar = %f\n", reuslt_scalar);
};


int main() {
    testMatrix2D();

    printf("Press Enter to continue...");
    while(getchar() != '\n');

    return 0;
}
