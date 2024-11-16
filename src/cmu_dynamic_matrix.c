
#include <stdlib.h>
#include <math.h>

#include "cmu_dynamic_matrix.h"

// 创建动态矩阵
MatrixXD CreateMatXD(int rows, int cols) {
    MatrixXD matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.data = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix.data[i] = (double*)malloc(cols * sizeof(double));
    }
    return matrix;
}

// 创建动态向量
VectorXD CreateVecXD(int size) {
    VectorXD vector;
    vector.size = size;
    vector.data = (double*)malloc(size * sizeof(double));
    return vector;
}

// 释放动态矩阵
void FreeMatXD(MatrixXD *matrix) {
    if (matrix && matrix->data) {
        for (int i = 0; i < matrix->rows; i++) {
            free(matrix->data[i]);
        }
        free(matrix->data);
        matrix->data = NULL; // 清空指针以避免悬挂指针
    }
}

// 释放向量
void FreeVecXD(VectorXD *vector) {
    if (vector && vector->data) {
        free(vector->data);
        vector->data = NULL;
    }
}

// 矩阵求和，结果存储在 result 中
Status AddMatXD(const MatrixXD *A, const MatrixXD *B, MatrixXD *result) {
    // 检查输入矩阵的维度是否匹配
    if (A->rows != B->rows || A->cols != B->cols || A->rows == 0 || B->rows == 0) {
        return CMU_ERROR_MATRIX_OPERATION_INVALID; 
    }
    if (A->rows != result->rows || A->cols != result->cols || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            result->data[i][j] = A->data[i][j] + B->data[i][j];
        }
    }

    return CMU_STATUS_SUCCESS;
}

// 矩阵求差
Status SubMatXD(const MatrixXD *A, const MatrixXD *B, MatrixXD *result) 
{
    if (A->rows != B->rows || A->cols != B->cols || A->rows == 0 || B->rows == 0) 
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID; 
    }
    if (A->rows != result->rows || A->cols != result->cols || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            result->data[i][j] = A->data[i][j] - B->data[i][j];
        }
    }
    return CMU_STATUS_SUCCESS;
}

// 矩阵乘法
Status DotMatXD(const MatrixXD *A, const MatrixXD *B, MatrixXD *result) {
    if (A->cols != B->rows || A->rows == 0 || B->rows == 0) 
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID; 
    }
    if (A->rows != result->rows || B->cols != result->cols || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->cols; j++) {
            result->data[i][j] = 0;
            for (int k = 0; k < A->cols; k++) {
                result->data[i][j] += A->data[i][k] * B->data[k][j];
            }
        }
    }
    return CMU_STATUS_SUCCESS;
}

// 矩阵数乘
Status DotScalarMatXD(const double scalar, const MatrixXD *matrix, MatrixXD *result) {
    if (matrix->rows != result->rows || matrix->cols != result->cols ||
        matrix->rows == 0 || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            result->data[i][j] = scalar * matrix->data[i][j];
        }
    }
    return CMU_STATUS_SUCCESS;
}

// 矩阵转置
Status TransposeMatXD(const MatrixXD *matrix, MatrixXD *result) {
    if (matrix->rows != result->cols || matrix->cols != result->rows ||
        matrix->rows == 0 || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            result->data[j][i] = matrix->data[i][j];
        }
    }
    return CMU_STATUS_SUCCESS;
}

// 计算矩阵的迹
Status TraceMatXD(const MatrixXD *matrix, double *result)
{
    double trace = 0.0;
    if (matrix->rows != matrix->cols) {
        // 如果矩阵不是方阵，迹无意义
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < matrix->rows; i++) {
        trace += matrix->data[i][i];
    }
    *result = trace;

    return CMU_STATUS_SUCCESS;
}

// 矩阵与向量相乘
Status DotMatVecXD(const MatrixXD *matrix, const VectorXD *vector, VectorXD *result) 
{
    if (matrix->cols != vector->size || vector->size != result->size || 
        matrix->rows == 0 || vector->size == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    
    for (int i = 0; i < matrix->rows; i++) {
        result->data[i] = 0.0;
        for (int j = 0; j < matrix->cols; j++) {
            result->data[i] += matrix->data[i][j] * vector->data[j];
        }
    }
    return CMU_STATUS_SUCCESS;
}


//向量数乘
Status DotScalarVecXD(double scalar, const VectorXD *vector, VectorXD *result)
{
    if (vector->size != result->size || vector->size == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < vector->size; i++) {
        result->data[i] = scalar * vector->data[i];
    }
    return CMU_STATUS_SUCCESS;
}

//向量点乘
Status DotVecXD(const VectorXD *vectorA, const VectorXD *vectorB, double *result)
{
    double dot = 0.0;
    if (vectorA->size != vectorB->size || vectorA->size == 0) 
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }

    for (int i = 0; i < vectorA->size; i++) {
        dot += vectorA->data[i] * vectorB->data[i];
    }
    *result = dot;

    return CMU_STATUS_SUCCESS;
}


Status DiagVecToMatXD(const VectorXD *vector, MatrixXD *result)
{
    if (vector->size != result->rows || vector->size != result->cols ||
        vector->size == 0 || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    
    for (int i = 0; i < vector->size; i++) {
        result->data[i][i] = vector->data[i];
    }
    return CMU_STATUS_SUCCESS;
}


Status NormVecXD(const VectorXD *vector, double *result)
{
    double norm = 0.0;
    if (vector->size == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < vector->size; i++) {
        norm += vector->data[i] * vector->data[i];
    }
    norm = sqrt(norm);
    *result = norm;

    return CMU_STATUS_SUCCESS;
}