/*
 * 文件名: cmu_dynamic_matrix.c
 * 简介: 矩阵向量运算头文件
 *
 * 本文件定义了动态矩阵和动态向量，并提供了一系列矩阵和向量的基本运算接口，
 * 包括加法、乘法、向量数乘、点乘、 向量的模长等计算等。
 *
 * @author [sanghao_nwpu]
 * @version 1.0.0
 * @date 2024年10月18日
 *
 * @section change_log 版本变更说明
 * - 1.0.0: 初始版本，添加了基本的向量和矩阵数据结构及运算接口。
 */

#ifndef _CMU_DYNAMIC_MATRIX_H_
#define _CMU_DYNAMIC_MATRIX_H_

#include "cmu_defs.h"

/**
 * @brief 创建动态矩阵
 * @param rows 矩阵的行数
 * @param cols 矩阵的列数
 * @param matrixPtr 指向动态矩阵的指针
 * @return 状态码
 */
Status CreateMatXD(const int rows, const int cols, MatrixXD *matrixPtr);

/**
 * @brief 创建动态向量
 * @param size 向量的长度
 * @param vectorPtr 指向动态向量的指针
 * @return 状态码
 */
Status CreateVecXD(const int size, VectorXD *vectorPtr);

/**
 * @brief 释放动态矩阵
 * @param matrix 指向待释放的动态矩阵的指针
 */
void FreeMatXD(MatrixXD *matrix);

/**
 * @brief 释放动态向量
 * @param vector 指向待释放的向量矩阵的指针
 */
void FreeVecXD(VectorXD *vector);

/**
 * @brief 初始化矩阵为单位矩阵
 * @param matrix 待初始化的矩阵，必须为已分配好内存的方阵
 * @return Status
 */
Status InitIdentityMatXD(MatrixXD *matrix);

/**
 * @brief 初始化矩阵为零矩阵
 * @param matrix 待初始化的矩阵，必须为已分配好内存的方阵
 * @return Status
 */
Status InitZeroMatXD(MatrixXD *matrix);

/**
 * @brief 初始化向量为零向量
 * @param vector 待初始化的向量，必须为已分配好内存的方阵
 * @return Status
 */
Status InitZeroVecXD(VectorXD *vector);

/**
 * @brief 矩阵求和
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status AddMatXD(const MatrixXD *A, const MatrixXD *B, MatrixXD *result);

/**
 * @brief 矩阵求差
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status SubMatXD(const MatrixXD *A, const MatrixXD *B, MatrixXD *result);

/**
 * @brief 矩阵乘法
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status DotMatXD(const MatrixXD *A, const MatrixXD *B, MatrixXD *result);

/**
 * @brief 矩阵数乘
 * @param scalar 乘以的数值
 * @param matrix 指向要进行数乘的矩阵的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status DotScalarMatXD(const double scalar, const MatrixXD *matrix, MatrixXD *result);

/**
 * @brief 矩阵转置
 * @param matrix 指向要转置的矩阵的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status TransposeMatXD(const MatrixXD *matrix, MatrixXD *result);

/**
 * @brief 计算矩阵的迹
 * @param matrix 指向要计算迹的矩阵的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status TraceMatXD(const MatrixXD *matrix, double *result);

/**
 * @brief 矩阵与向量相乘
 * @param matrix 指向要进行相乘的矩阵的指针
 * @param vector 指向要进行相乘的向量的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status DotMatVecXD(const MatrixXD *matrix, const VectorXD *vector, VectorXD *result);

/**
 * @brief 向量数乘
 * @param scalar 乘以的数值
 * @param vector 指向要进行数乘的向量的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status DotScalarVecXD(double scalar, const VectorXD *vector, VectorXD *result);

/**
 * @brief 向量点乘
 * @param vectorA 指向第一个向量的指针
 * @param vectorB 指向第二个向量的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status DotVecXD(const VectorXD *vectorA, const VectorXD *vectorB, double *result);

/**
 * @brief 向量转换为对角矩阵
 * @param vector 指向要转换的向量的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status DiagVecToMatXD(const VectorXD *vector, MatrixXD *result);

/**
 * @brief 向量模长计算
 * @param vector 指向要进行操作的向量的指针
 * @param result 指向结果的指针
 * @return 状态码
 */
Status NormVecXD(const VectorXD *vector, double *result);

/**
 * @brief 矩阵求逆
 * @param matrix 指向要进行求逆的矩阵的指针
 * @param result 指向结果的指针
 * @return Status
 */
Status InverseMatXD(const MatrixXD *matrix, MatrixXD *result);

#endif // _CMU_DYNAMIC_MATRIX_H_