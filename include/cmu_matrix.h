#ifndef _CMATHUTILS_MATRIX_H_
#define _CMATHUTILS_MATRIX_H_

/**
 * @file cmu_matrix.h
 * @brief 固定大小的矩阵向量运算头文件
 *
 * 本文件定义了二维、三维和四维的向量及对应的矩阵，并提供了一系列
 * 矩阵和向量的基本运算接口，包括加减法、乘法、向量数乘、叉乘、点乘、
 * 模长等计算等。
 *
 * @author [sanghao_nwpu]
 * @version 1.0.0
 * @date 2024年10月18日
 *
 * @section change_log 版本变更说明
 * - 1.0.0: 初始版本，添加了基本的向量和矩阵数据结构及运算接口。
 */

#include "cmu_defs.h"

/*******************************************************************************
 * @subsection 矩阵转置
*******************************************************************************/
/**
 * @brief 矩阵转置(2维)
 * @param A 指向矩阵的指针
 * @return 转置后的矩阵
 */
Matrix2D MatTranspose2D(const Matrix2D *A);

/**
 * @brief 矩阵转置(3维)
 * @param A 指向矩阵的指针
 * @return 转置后的矩阵
 */
Matrix3D MatTranspose3D(const Matrix3D *A);

/**
 * @brief 矩阵转置(4维)
 * @param A 指向矩阵的指针
 * @return 转置后的矩阵
 */
Matrix4D MatTranspose4D(const Matrix4D *A);

/*******************************************************************************
 * @subsection 矩阵加法
*******************************************************************************/
/**
 * @brief 矩阵加法(2维)
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @return 矩阵的和
 */
Matrix2D MatAdd2D(const Matrix2D *A, const Matrix2D *B);

/**
 * @brief 矩阵加法(3维)
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @return 矩阵的和
 */
Matrix3D MatAdd3D(const Matrix3D *A, const Matrix3D *B);

/**
 * @brief 矩阵加法(3维)
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @return 矩阵的和
 */
Matrix4D MatAdd4D(const Matrix4D *A, const Matrix4D *B);


/*******************************************************************************
 * @subsection 矩阵减法
*******************************************************************************/
/**
 * @brief 矩阵减法(2维)
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @return 矩阵的差
 */
Matrix2D MatSub2D(const Matrix2D *A, const Matrix2D *B);

/**
 * @brief 矩阵减法(3维)
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @return 矩阵的差
 */
Matrix3D MatSub3D(const Matrix3D *A, const Matrix3D *B);

/**
 * @brief 矩阵减法(4维)
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @return 矩阵的差
 */
Matrix4D MatSub4D(const Matrix4D *A, const Matrix4D *B);


/*******************************************************************************
 * @subsection 矩阵与矩阵相乘
*******************************************************************************/
/**
 * @brief 矩阵乘法(2维)
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @return 矩阵的乘积
 */
Matrix2D MatMul2D(const Matrix2D *A, const Matrix2D *B);

/**
 * @brief 矩阵乘法(3维)
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @return 矩阵的乘积
 */
Matrix3D MatMul3D(const Matrix3D *A, const Matrix3D *B);

/**
 * @brief 矩阵乘法(4维)
 * @param A 指向第一个矩阵的指针
 * @param B 指向第二个矩阵的指针
 * @return 矩阵的乘积
 */
Matrix4D MatMul4D(const Matrix4D *A, const Matrix4D *B);


/*******************************************************************************
 * @subsection 标量与矩阵相乘
*******************************************************************************/
/**
 * @brief 标量与矩阵相乘(2维)
 * @param A 指向矩阵的指针
 * @param scalar 要乘的标量
 * @return 标量和矩阵的乘积
 */
Matrix2D ScalarMulMat2D(const double scalar, const Matrix2D *A);

/**
 * @brief 标量与矩阵相乘(3维)
 * @param A 指向矩阵的指针
 * @param scalar 要乘的标量
 * @return 标量和矩阵的乘积
 */
Matrix3D ScalarMulMat3D(const double scalar, const Matrix3D *A);

/**
 * @brief 标量与矩阵相乘(4维)
 * @param A 指向矩阵的指针
 * @param scalar 要乘的标量
 * @return 标量和矩阵的乘积
 */
Matrix4D ScalarMulMat4D(const double scalar, const Matrix4D *A);


/*******************************************************************************
 * @subsection 矩阵与向量相乘
*******************************************************************************/
/**
 * @brief 矩阵与向量相乘(2维)
 * @param A 指向矩阵的指针
 * @param v 指向向量的指针
 * @return 矩阵和向量的乘积
 */
Vector2D MatMulVec2D(const Matrix2D *A, const Vector2D *v);

/**
 * @brief 矩阵与向量相乘(3维)
 * @param A 指向矩阵的指针
 * @param v 指向向量的指针
 * @return 矩阵和向量的乘积
 */
Vector3D MatMulVec3D(const Matrix3D *A, const Vector3D *v);

/**
 * @brief 矩阵与向量相乘(4维)
 * @param A 指向矩阵的指针
 * @param v 指向向量的指针
 * @return 矩阵和向量的乘积
 */
Vector4D MatMulVec4D(const Matrix4D *A, const Vector4D *v);


/*******************************************************************************
 * @subsection 向量加法
*******************************************************************************/
/**
 * @brief 向量加法(2维)
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量的和
 */
Vector2D VecAdd2D(const Vector2D *u, const Vector2D *v);

/**
 * @brief 向量加法(3维)
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量的和
 */
Vector3D VecAdd3D(const Vector3D *u, const Vector3D *v);

/**
 * @brief 向量加法(4维)
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量的和
 */
Vector4D VecAdd4D(const Vector4D *u, const Vector4D *v);

/*******************************************************************************
 * @subsection 向量减法
*******************************************************************************/
/**
 * @brief 向量减法(2维)
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量的差
 */
Vector2D VecSub2D(const Vector2D *u, const Vector2D *v);

/**
 * @brief 向量减法(3维)
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量的差
 */
Vector3D VecSub3D(const Vector3D *u, const Vector3D *v);

/**
 * @brief 向量减法(4维)
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量的差
 */
Vector4D VecSub4D(const Vector4D *u, const Vector4D *v);


/*******************************************************************************
 * @subsection 向量数乘
*******************************************************************************/
/**
 * @brief 向量数乘(2维)
 * @param v 指向向量的指针
 * @param scalar 要乘的标量
 * @return 向量数乘结果
 */
Vector2D ScalarMulVec2D(const double scalar, const Vector2D *v);

/**
 * @brief 向量数乘(3维)
 * @param v 指向向量的指针
 * @param scalar 要乘的标量
 * @return 向量数乘结果
 */
Vector3D ScalarMulVec3D(const double scalar, const Vector3D *v);

/**
 * @brief 向量数乘(4维)
 * @param v 指向向量的指针
 * @param scalar 要乘的标量
 * @return 向量数乘结果
 */
Vector4D ScalarMulVec4D(const double scalar, const Vector4D *v);


/*******************************************************************************
 * @subsection 向量点乘
*******************************************************************************/
/**
 * @brief 向量点乘
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量点乘结果
 */
double VecDot2D(const Vector2D *u, const Vector2D *v);

/**
 * @brief 向量点乘
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量点乘结果
 */
double VecDot3D(const Vector3D *u, const Vector3D *v);

/**
 * @brief 向量点乘
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量点乘结果
 */
double VecDot4D(const Vector4D *u, const Vector4D *v);


/*******************************************************************************
 * @subsection 向量叉乘(仅适用于三维向量)
*******************************************************************************/
/**
 * @brief 向量叉乘(仅适用于3维向量)
 * 
 * @param u 指向第一个向量的指针
 * @param v 指向第二个向量的指针
 * @return 向量叉乘结果
 */
Vector3D VecCross3D(const Vector3D *u, const Vector3D *v);


/*******************************************************************************
 * @subsection 向量和反对称矩阵互转(仅适用于3维向量)
*******************************************************************************/
/**
 * @brief 将3维向量转换为反对称矩阵
 * @param v 指向向量的指针
 * @return 该向量的反对称矩阵
 */
Matrix3D VecToSkew3D(const Vector3D *v);

/**
 * @brief 将3维反对称矩阵转换为向量
 * @param v 指向反对称矩阵的指针
 * @return 该反对称矩阵对应的向量
 */
Vector3D SkewToVec3D(const Matrix3D *skew);


/*******************************************************************************
 * @subsection 向量和对角矩阵互转
*******************************************************************************/
/**
 * @brief 将向量转换为对角矩阵(2维)
 * @param v 指向向量的指针
 * @return 向量对应的对角矩阵
 */
Matrix2D VecToDiag2D(const Vector2D *v);

/**
 * @brief 将向量转换为对角矩阵(3维)
 * @param v 指向向量的指针
 * @return 向量对应的对角矩阵
 */
Matrix3D VecToDiag3D(const Vector3D *v);

/**
 * @brief 将向量转换为对角矩阵(4维)
 * @param v 指向向量的指针
 * @return 向量对应的对角矩阵
 */
Matrix4D VecToDiag4D(const Vector4D *v);

/**
 * @brief 将对角矩阵转换为向量(2维)
 * @param v 指向对角的指针
 * @return 对角矩阵对应的向量
 */
Vector2D DiagToVec2D(const Matrix2D *diag);

/**
 * @brief 将对角矩阵转换为向量(3维)
 * @param v 指向对角的指针
 * @return 对角矩阵对应的向量
 */
Vector3D DiagToVec3D(const Matrix3D *diag);

/**
 * @brief 将对角矩阵转换为向量(2维)
 * @param v 指向对角的指针
 * @return 对角矩阵对应的向量
 */
Vector4D DiagToVec4D(const Matrix4D *diag);


/*******************************************************************************
 * @subsection 向量按位取平方
*******************************************************************************/
/**
 * @brief 向量按位取平方(2维)
 * @param v 指向向量的指针
 * @return 按位取平方后的向量
 */
Vector2D VecSquare2D(const Vector2D *v);

/**
 * @brief 向量按位取平方(3维)
 * @param v 指向向量的指针
 * @return 按位取平方后的向量
 */
Vector3D VecSquare3D(const Vector3D *v);

/**
 * @brief 向量按位取平方(4维)
 * @param v 指向向量的指针
 * @return 按位取平方后的向量
 */
Vector4D VecSquare4D(const Vector4D *v);


/*******************************************************************************
 * @subsection 向量按位取开方
*******************************************************************************/
/**
 * @brief 向量按位取开方(2维)
 * @param v 指向向量的指针
 * @return 按位取开方后的向量
 */
Vector2D VecSqrt2D(const Vector2D *v);

/**
 * @brief 向量按位取开方(3维)
 * @param v 指向向量的指针
 * @return 按位取开方后的向量
 */
Vector3D VecSqrt3D(const Vector3D *v);

/**
 * @brief 向量按位取开方(4维)
 * @param v 指向向量的指针
 * @return 按位取开方后的向量
 */
Vector4D VecSqrt4D(const Vector4D *v);


/*******************************************************************************
 * @subsection 计算向量模长
*******************************************************************************/
/**
 * @brief 计算向量的模长(2维)
 * @param v 指向向量的指针
 * @return 返回向量的模长
 */
double VecNorm2D(const Vector2D *v);

/**
 * @brief 计算向量的模长(3维)
 * @param v 指向向量的指针
 * @return 返回向量的模长
 */
double VecNorm3D(const Vector3D *v);

/**
 * @brief 计算向量的模长(4维)
 * @param v 指向向量的指针
 * @return 返回向量的模长
 */
double VecNorm4D(const Vector4D *v);


#endif // _CMATHUTILS_FIXED_MATRIX_H_
