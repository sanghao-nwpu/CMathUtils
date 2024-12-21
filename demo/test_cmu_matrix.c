#include <stdio.h>
#include "cmu_static_matrix.h"
#include "cmu_dynamic_matrix.h"

void printMatrixXD(MatrixXD *matrix) 
{
    for (size_t i = 0; i < matrix->rows; i++) {
        for (size_t j = 0; j < matrix->cols; j++) {
            printf("%lf ", matrix->data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void testMatrixXD()
{
    // 创建三维矩阵和向量
    MatrixXD matA, matB, resultMat;
    VectorXD vecA, vecB, resultVec;
    MatrixXD matC, matD;
    CreateMatXD(4, 4, &matC);
    CreateMatXD(4, 4, &matD);
    CreateMatXD(3, 3, &matA);
    CreateMatXD(3, 3, &matB);
    CreateMatXD(3, 3, &resultMat);
    CreateVecXD(3, &vecA);
    CreateVecXD(3, &vecB);
    CreateVecXD(3, &resultVec);
    Status status;
    double traceResult;
    double dotResult;

    // 填充矩阵A
    matA.data[0][0] = 1; matA.data[0][1] = 0; matA.data[0][2] = 1;
    matA.data[1][0] = 0; matA.data[1][1] = 2; matA.data[1][2] = 0;
    matA.data[2][0] = 1; matA.data[2][1] = 0; matA.data[2][2] = 3;

    // 填充矩阵B
    matB.data[0][0] = 9; matB.data[0][1] = 8; matB.data[0][2] = 7;
    matB.data[1][0] = 6; matB.data[1][1] = 5; matB.data[1][2] = 4;
    matB.data[2][0] = 3; matB.data[2][1] = 2; matB.data[2][2] = 1;

    // 填充矩阵C
    matC.data[0][0] = 1; matC.data[0][1] = 2; matC.data[0][2] = 0; matC.data[0][3] = 0;
    matC.data[1][0] = 0; matC.data[1][1] = 2; matC.data[1][2] = 0; matC.data[1][3] = 0;
    matC.data[2][0] = 0; matC.data[2][1] = 0; matC.data[2][2] = 3; matC.data[2][3] = 0;
    matC.data[3][0] = 2; matC.data[3][1] = 0; matC.data[3][2] = 0; matC.data[3][3] = 4;

    // 填充向量A
    vecA.data[0] = 1; vecA.data[1] = 2; vecA.data[2] = 3;

    // 填充向量B
    vecB.data[0] = 4; vecB.data[1] = 5; vecB.data[2] = 6;

    // 矩阵求和
    AddMatXD(&matA, &matB, &resultMat);
    printf("Matrix A + Matrix B:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", resultMat.data[i][j]);
        }
        printf("\n");
    }

    // 矩阵求差
    SubMatXD(&matA, &matB, &resultMat);
    printf("Matrix A - Matrix B:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", resultMat.data[i][j]);
        }
        printf("\n");
    }

    // 矩阵乘法
    DotMatXD(&matA, &matB, &resultMat);
    printf("Matrix A * Matrix B:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", resultMat.data[i][j]);
        }
        printf("\n");
    }

    // 矩阵数乘
    DotScalarMatXD(2.0, &matA, &resultMat);
    printf("Matrix A * 2.0:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", resultMat.data[i][j]);
        }
        printf("\n");
    }

    // 矩阵转置
    TransposeMatXD(&matA, &resultMat);
    printf("Transpose of Matrix A:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", resultMat.data[i][j]);
        }
        printf("\n");
    }

    // 矩阵原地转置
    TransposeMatXD(&resultMat, &resultMat);
    printf("Transpose of Matrix resultMat:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", resultMat.data[i][j]);
        }
        printf("\n");
    }


    // 计算矩阵的迹
    TraceMatXD(&matA, &traceResult);
    printf("Trace of Matrix A: %lf\n", traceResult);

    // 测试逆矩阵计算
    status = InverseMatXD(&matA, &resultMat);
    if (status != CMU_STATUS_SUCCESS) {
        printf("Matrix Inversion failed with status: %d\n", status);
    }
    printf("Inverse Matrix:\n");
    printMatrixXD(&resultMat);

    status = InverseMatXD(&matC, &matD);
    if (status != CMU_STATUS_SUCCESS) {
        printf("Matrix Inversion failed with status: %d\n", status);
    }
    printf("Inverse Matrix:\n");
    printMatrixXD(&matD);

    // 矩阵与向量相乘
    DotMatVecXD(&matA, &vecA, &resultVec);
    printf("Matrix A * Vector A:\n");
    for (int i = 0; i < 3; i++) {
        printf("%lf ", resultVec.data[i]);
    }
    printf("\n");

    // 向量数乘
    DotScalarVecXD(3.0, &vecA, &resultVec);
    printf("Vector A * 3.0:\n");
    for (int i = 0; i < 3; i++) {
        printf("%lf ", resultVec.data[i]);
    }
    printf("\n");

    // 向量点乘
    DotVecXD(&vecA, &vecB, &dotResult);
    printf("Dot product of Vector A and Vector B: %lf\n", dotResult);

    // 向量转换为对角矩阵
    DiagVecToMatXD(&vecA, &resultMat);
    printf("Diagonal matrix from Vector A:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", resultMat.data[i][j]);
        }
        printf("\n");
    }

    // 向量模长计算
    NormVecXD(&vecA, &dotResult);
    printf("Norm of Vector A: %lf\n", dotResult);

    // 释放动态矩阵和向量
    FreeMatXD(&matA);
    FreeMatXD(&matB);
    FreeVecXD(&vecA);
    FreeVecXD(&vecB);
    FreeVecXD(&resultVec);
    FreeMatXD(&resultMat);

}

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
    // testMatrix2D();
    testMatrixXD();
    printf("Press Enter to continue...");
    while(getchar() != '\n');

    return 0;
}
