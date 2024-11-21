
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
    for (int i = 0; i < matrix.rows; i++) 
    {
        for (int j = 0; j < matrix.cols; j++) 
        {
            matrix.data[i][j] = 0.0;
        }
    }
    return matrix;
}

// 创建动态向量
VectorXD CreateVecXD(int size) {
    VectorXD vector;
    vector.size = size;
    vector.data = (double*)malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) {
        vector.data[i] = 0.0;
    }
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

Status InitIdentityMatXD(MatrixXD *matrix) 
{
    if (matrix->rows != matrix->cols || matrix->rows == 0) {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < matrix->rows; i++) 
    {
        for (int j = 0; j < matrix->cols; j++) 
        {
            if (i == j) {
                matrix->data[i][j] = 1.0;
            } else {
                matrix->data[i][j] = 0.0;
            }
        }
    }
    return CMU_STATUS_SUCCESS;
}

Status InitZeroMatXD(MatrixXD *matrix)
{
    if (matrix->rows == 0 || matrix->cols == 0) {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            matrix->data[i][j] = 0.0;
        }
    }
    return CMU_STATUS_SUCCESS;
}

Status InitZeroVecXD(VectorXD *vector)
{
    if (vector->size == 0) {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < vector->size; i++) {
        vector->data[i] = 0.0;
    }
    return CMU_STATUS_SUCCESS;
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

    // 临时存储转置的结果
    double **tempData = (double **)malloc(result->rows * sizeof(double *));
    for (int i = 0; i < result->rows; i++) {
        tempData[i] = (double *)malloc(result->cols * sizeof(double));
    }

    // 进行转置
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            tempData[j][i] = matrix->data[i][j];
        }
    }

    // 将转置结果写回到 result 中
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            result->data[i][j] = tempData[i][j];
        }
    }

    // 释放临时存储
    for (int i = 0; i < result->rows; i++) {
        free(tempData[i]);
    }
    free(tempData);
    
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
    
    for (int i = 0; i < result->rows; i++)
    {
        for (int j = 0; j < result->cols; j++)
        {
            if (i == j)
            {
                result->data[i][j] = vector->data[i];
            }
            else
            {
                result->data[i][j] = 0.0;
            }
        }
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


Status LUDecomposition(const MatrixXD *matrix, MatrixXD* L, MatrixXD* U) {
    if (matrix->rows != matrix->cols || matrix->rows == 0 ||
        L->rows != L->cols || L->rows != matrix->rows ||
        U->rows != U->cols || U->rows != matrix->rows ||
        L->rows != matrix->rows || U->cols != matrix->rows)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }

    for (int i = 0; i < matrix->rows; i++) 
    {
        for (int j = i; j < matrix->cols; j++) 
        {
            U->data[i][j] = matrix->data[i][j];
            for (int k = 0; k < i; k++) {
                U->data[i][j] -= L->data[i][k] * U->data[k][j];
            }
        }
        
        for (int j = i + 1; j < matrix->cols; j++) 
        {
            L->data[j][i] = matrix->data[j][i];
            for (int k = 0; k < i; k++) 
            {
                L->data[j][i] -= L->data[j][k] * U->data[k][i];
            }
            L->data[j][i] /= U->data[i][i];
        }
    }
    return CMU_STATUS_SUCCESS;
}


Status InverseMatXD(const MatrixXD *matrix, MatrixXD *result) {
    int n = matrix->rows;
    Status status;
    MatrixXD L, U;
    VectorXD e, x, y;
    
    if (matrix->rows != matrix->cols || matrix->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    L = CreateMatXD(n, n);
    U = CreateMatXD(n, n);
    e = CreateVecXD(n);
    x = CreateVecXD(n);
    y = CreateVecXD(n);

    // 初始化 L 为单位矩阵
    status = InitIdentityMatXD(&L);
    if (CMU_STATUS_SUCCESS != status)
    {
        FreeMatXD(&L);
        FreeMatXD(&U);
        return status;
    }

    // LU分解
    status = LUDecomposition(matrix, &L, &U);
    if (CMU_STATUS_SUCCESS != status)
    {
        FreeMatXD(&L);
        FreeMatXD(&U);
        return status;
    }

    // 对每一列进行回代，构造逆矩阵
    for (int i = 0; i < n; i++) 
    {
        // 前向替代解 L * y = e
        status = InitZeroVecXD(&e);
        e.data[i] = 1.0;
        for (int j = 0; j < n; j++) 
        {
            y.data[j] = e.data[j];
            for (int k = 0; k < j; k++) {
                y.data[j] -= L.data[j][k] * y.data[k];
            }
        }
        
        // 后向替代解 U * x = y
        for (int i = n - 1; i >= 0; i--) {
            x.data[i] = y.data[i];
            for (int j = i + 1; j < n; j++) {
                x.data[i] -= U.data[i][j] * x.data[j];
            }
            x.data[i] /= U.data[i][i];
        }

        // 填充逆矩阵
        for (int j = 0; j < n; j++) {
            result->data[j][i] = x.data[j];
        }

        FreeVecXD(&e);
        FreeVecXD(&x);
        FreeVecXD(&y);
    }
    FreeMatXD(&L);
    FreeMatXD(&U);

    return CMU_STATUS_SUCCESS;
}