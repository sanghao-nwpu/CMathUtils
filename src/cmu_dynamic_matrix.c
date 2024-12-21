
#include <stdlib.h>
#include <math.h>

#include "cmu_dynamic_matrix.h"

// 创建动态矩阵
Status CreateMatXD(const int rows, const int cols, MatrixXD *matrixPtr)
{
    if (rows == 0 || cols == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    matrixPtr->rows = rows;
    matrixPtr->cols = cols;
    matrixPtr->data = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        matrixPtr->data[i] = (double *)malloc(cols * sizeof(double));
    }
    for (int i = 0; i < matrixPtr->rows; i++)
    {
        for (int j = 0; j < matrixPtr->cols; j++)
        {
            matrixPtr->data[i][j] = 0.0;
        }
    }
    return CMU_STATUS_SUCCESS;
}

// 创建动态向量
Status CreateVecXD(const int size, VectorXD *vectorPtr)
{
    if (size == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    vectorPtr->size = size;
    vectorPtr->data = (double *)malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
    {
        vectorPtr->data[i] = 0.0;
    }
    return CMU_STATUS_SUCCESS;
}

// 释放动态矩阵
void FreeMatXD(MatrixXD *matrix)
{
    if (matrix && matrix->data)
    {
        for (int i = 0; i < matrix->rows; i++)
        {
            free(matrix->data[i]);
        }
        free(matrix->data);
        matrix->data = NULL; // 清空指针以避免悬挂指针
    }
}

// 释放向量
void FreeVecXD(VectorXD *vector)
{
    if (vector && vector->data)
    {
        free(vector->data);
        vector->data = NULL;
    }
}

Status InitIdentityMatXD(MatrixXD *matrix)
{
    if (matrix->rows != matrix->cols || matrix->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            if (i == j)
            {
                matrix->data[i][j] = 1.0;
            }
            else
            {
                matrix->data[i][j] = 0.0;
            }
        }
    }
    return CMU_STATUS_SUCCESS;
}

Status InitZeroMatXD(MatrixXD *matrix)
{
    if (matrix->rows == 0 || matrix->cols == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            matrix->data[i][j] = 0.0;
        }
    }
    return CMU_STATUS_SUCCESS;
}

Status InitZeroVecXD(VectorXD *vector)
{
    if (vector->size == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < vector->size; i++)
    {
        vector->data[i] = 0.0;
    }
    return CMU_STATUS_SUCCESS;
}

// 矩阵求和，结果存储在 result 中
Status AddMatXD(const MatrixXD *A, const MatrixXD *B, MatrixXD *result)
{
    // 检查输入矩阵的维度是否匹配
    if (A->rows != B->rows || A->cols != B->cols || A->rows == 0 || B->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    if (A->rows != result->rows || A->cols != result->cols || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }

    for (int i = 0; i < A->rows; i++)
    {
        for (int j = 0; j < A->cols; j++)
        {
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

    for (int i = 0; i < A->rows; i++)
    {
        for (int j = 0; j < A->cols; j++)
        {
            result->data[i][j] = A->data[i][j] - B->data[i][j];
        }
    }
    return CMU_STATUS_SUCCESS;
}

// 矩阵乘法
Status DotMatXD(const MatrixXD *A, const MatrixXD *B, MatrixXD *result)
{
    if (A->cols != B->rows || A->rows == 0 || B->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    if (A->rows != result->rows || B->cols != result->cols || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }

    for (int i = 0; i < A->rows; i++)
    {
        for (int j = 0; j < B->cols; j++)
        {
            result->data[i][j] = 0;
            for (int k = 0; k < A->cols; k++)
            {
                result->data[i][j] += A->data[i][k] * B->data[k][j];
            }
        }
    }
    return CMU_STATUS_SUCCESS;
}

// 矩阵数乘
Status DotScalarMatXD(const double scalar, const MatrixXD *matrix, MatrixXD *result)
{
    if (matrix->rows != result->rows || matrix->cols != result->cols ||
        matrix->rows == 0 || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            result->data[i][j] = scalar * matrix->data[i][j];
        }
    }
    return CMU_STATUS_SUCCESS;
}

// 矩阵转置
Status TransposeMatXD(const MatrixXD *matrix, MatrixXD *result)
{
    if (matrix->rows != result->cols || matrix->cols != result->rows ||
        matrix->rows == 0 || result->rows == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }

    // 临时存储转置的结果
    double **tempData = (double **)malloc(result->rows * sizeof(double *));
    for (int i = 0; i < result->rows; i++)
    {
        tempData[i] = (double *)malloc(result->cols * sizeof(double));
    }

    // 进行转置
    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            tempData[j][i] = matrix->data[i][j];
        }
    }

    // 将转置结果写回到 result 中
    for (int i = 0; i < result->rows; i++)
    {
        for (int j = 0; j < result->cols; j++)
        {
            result->data[i][j] = tempData[i][j];
        }
    }

    // 释放临时存储
    for (int i = 0; i < result->rows; i++)
    {
        free(tempData[i]);
    }
    free(tempData);

    return CMU_STATUS_SUCCESS;
}

// 计算矩阵的迹
Status TraceMatXD(const MatrixXD *matrix, double *result)
{
    double trace = 0.0;
    if (matrix->rows != matrix->cols)
    {
        // 如果矩阵不是方阵，迹无意义
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < matrix->rows; i++)
    {
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

    for (int i = 0; i < matrix->rows; i++)
    {
        result->data[i] = 0.0;
        for (int j = 0; j < matrix->cols; j++)
        {
            result->data[i] += matrix->data[i][j] * vector->data[j];
        }
    }
    return CMU_STATUS_SUCCESS;
}

// 向量数乘
Status DotScalarVecXD(double scalar, const VectorXD *vector, VectorXD *result)
{
    if (vector->size != result->size || vector->size == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    for (int i = 0; i < vector->size; i++)
    {
        result->data[i] = scalar * vector->data[i];
    }
    return CMU_STATUS_SUCCESS;
}

// 向量点乘
Status DotVecXD(const VectorXD *vectorA, const VectorXD *vectorB, double *result)
{
    double dot = 0.0;
    if (vectorA->size != vectorB->size || vectorA->size == 0)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }

    for (int i = 0; i < vectorA->size; i++)
    {
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
    for (int i = 0; i < vector->size; i++)
    {
        norm += vector->data[i] * vector->data[i];
    }
    norm = sqrt(norm);
    *result = norm;

    return CMU_STATUS_SUCCESS;
}

/**
 * @brief 高斯消元法求逆
 * @note 缺点是数值上不够稳定
 */
// Status InverseMatXD(const MatrixXD *matrix, MatrixXD *result) {
//     if (matrix->rows != matrix->cols || matrix->rows == 0 ||
//         result->rows != result->cols || result->rows != matrix->rows)
//     {
//         return CMU_ERROR_MATRIX_OPERATION_INVALID;
//     }

//     int n = matrix->rows;
//     // 创建扩展矩阵
//     MatrixXD augmentedMatrix;
//     augmentedMatrix = CreateMatXD(n, 2 * n);

//     // 扩展矩阵的左部分为原矩阵，右部分为单位矩阵
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             augmentedMatrix.data[i][j] = matrix->data[i][j];
//         }
//         for (int j = n; j < 2 * n; j++) {
//             augmentedMatrix.data[i][j] = (i == j - n) ? 1.0 : 0.0;
//         }
//     }

//     // 高斯消元法
//     for (int i = 0; i < n; i++)
//     {
//         // 寻找主元
//         for (int k = i + 1; k < n; k++)
//         {
//             // 找绝对值最大的元素swap到主元位置
//             if (fabs(augmentedMatrix.data[k][i]) > fabs(augmentedMatrix.data[i][i]))
//             {
//                 // 交换行
//                 double* temp = augmentedMatrix.data[i];
//                 augmentedMatrix.data[i] = augmentedMatrix.data[k];
//                 augmentedMatrix.data[k] = temp;
//             }
//         }

//         // 确保主元不为零
//         if (augmentedMatrix.data[i][i] == 0)
//         {
//             FreeMatXD(&augmentedMatrix);
//             return CMU_STATUS_FAILURE;
//         }

//         // 将主元化为1
//         double pivot = augmentedMatrix.data[i][i];
//         for (int j = 0; j < 2 * n; j++) {
//             augmentedMatrix.data[i][j] /= pivot;
//         }

//         // 消元
//         for (int k = 0; k < n; k++) {
//             if (k != i) {
//                 double factor = augmentedMatrix.data[k][i];
//                 for (int j = 0; j < 2 * n; j++) {
//                     augmentedMatrix.data[k][j] -= factor * augmentedMatrix.data[i][j];
//                 }
//             }
//         }
//     }

//     // 输出结果
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             result->data[i][j] = augmentedMatrix.data[i][j + n];
//         }
//     }

//     FreeMatXD(&augmentedMatrix);
//     return CMU_STATUS_SUCCESS;
// }

Status InverseMatXD(const MatrixXD *matrix, MatrixXD *result)
{
    if (matrix->rows != matrix->cols || matrix->rows == 0 ||
        result->rows != result->cols || result->rows != matrix->rows)
    {
        return CMU_ERROR_MATRIX_OPERATION_INVALID;
    }
    MatrixXD L, U;
    VectorXD Y;
    int n = matrix->rows;
    CreateMatXD(n, n, &L);
    CreateMatXD(n, n, &U);
    CreateVecXD(n, &Y);

    // LU分解过程
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            U.data[i][j] = matrix->data[i][j]; // 上三角矩阵
            for (int k = 0; k < i; k++)
            {
                U.data[i][j] -= L.data[i][k] * U.data[k][j];
            }
        }
        for (int j = i; j < n; j++)
        {
            if (i == j)
            {
                L.data[i][i] = 1; // 对角线元素为1
            }
            else
            {
                L.data[j][i] = matrix->data[j][i]; // 下三角矩阵
                for (int k = 0; k < i; k++)
                {
                    L.data[j][i] -= L.data[j][k] * U.data[k][i];
                }
                L.data[j][i] /= U.data[i][i]; // 除以主元
            }
        }
    }

    // 依次求逆矩阵的每一列：首先求解 L*Y_i = e_i，然后 U*X_i = Y_i,X_i即是逆矩阵的第i列
    for (int i = 0; i < n; i++)
    {
        // 求解 L*Y = I，初始化Y为e_i
        InitZeroVecXD(&Y);
        Y.data[i] = 1;

        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < j; k++)
            {
                Y.data[j] -= L.data[j][k] * Y.data[k];
            }
        }
        for (int j = 0; j < n; j++)
        {
            result->data[j][i] = Y.data[j];
        }

        // 求解 U*X = Y；k为U的行，j为U的列
        for (int j = n - 1; j >= 0; j--)
        {
            result->data[j][i] /= U.data[j][j]; // 除以对角线元素
            for (int k = j - 1; k >= 0; k--)
            {
                result->data[k][i] -= U.data[k][j] * result->data[j][i];
            }
        }
    }

    // 释放内存
    FreeMatXD(&L);
    FreeMatXD(&U);
    FreeVecXD(&Y);
    return CMU_STATUS_SUCCESS;
}