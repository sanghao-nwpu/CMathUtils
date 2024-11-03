#ifndef _CMATHUTILS_ROTATION_H_
#define _CMATHUTILS_ROTATION_H_

/** 
 * 文件名: cmu_rotation.h
 * 简介: 旋转表示与计算的头文件
 *
 * 本文件定义了用于表示旋转的不同数据类型，包括欧拉角、角轴、
 * 四元数和旋转矩阵，并提供了这些旋转表示之间的互转接口、求逆运算等。
 *
 * @author [sanghao_nwpu]
 * @version 1.0.0
 * @date 2024年10月18日
 * @note change_log 版本变更说明
 * - 1.0.0: 初始版本，添加了旋转数据结构及互转接口的定义。
 */

#include "cmu_defs.h"

/*******************************************************************************
 * @subsection 旋转表示互转接口
*******************************************************************************/
/** 
 * @brief 将欧拉角转换为旋转矩阵 
 * @param euler 欧拉角结构体指针
 * @return RotationMatrixd 旋转矩阵
 */
RotationMatrixd EulerToRotMat(const EulerAngled *euler);

/** 
 * @brief 将旋转矩阵转换为欧拉角  
 * @param R 旋转矩阵结构体指针
 * @return EulerAngled 欧拉角
 */
EulerAngled RotMatToEuler(const RotationMatrixd *R);

/** 
 * @brief 将角轴转换为旋转矩阵 
 * @param aa 角轴结构体指针
 * @return RotationMatrixd 旋转矩阵
 */
RotationMatrixd AngleAxisToRotMat(const AngleAxisd *aa);

/** 
 * @brief 将旋转矩阵转换为角轴 
 * @param R 旋转矩阵结构体指针
 * @return AngleAxisd 角轴
 */
AngleAxisd RotMatToAngleAxis(const RotationMatrixd *R);

/** 
 * @brief 将四元数转换为旋转矩阵 
 * @param q 四元数结构体指针
 * @return RotationMatrixd 旋转矩阵
 */
RotationMatrixd QuatToRotMat(const Quaterniod *q);

/** 
 * @brief 将旋转矩阵转换为四元数 
 * @param R 旋转矩阵结构体指针
 * @return Quaterniod 四元数
 */
Quaterniod RotMatToQuat(const RotationMatrixd *R);

/**
 * @brief 将欧拉角转换为四元数 
 * @param euler 欧拉角结构体指针
 * @return Quaterniod 四元数
 */
Quaterniod EulerToQuat(const EulerAngled *euler);

/** 
 * @brief 将四元数转换为欧拉角 
 * @param q 四元数结构体指针
 * @return EulerAngled 欧拉角
 */
EulerAngled QuatToEuler(const Quaterniod *q);

/** 
 * @brief 将角轴转换为四元数 
 * @param aa 角轴结构体指针
 * @return Quaterniod 四元数
 */
Quaterniod AngleAxisToQuat(const AngleAxisd *aa);

/** 
 * @brief 将四元数转换为角轴 
 * @param q 四元数结构体指针
 * @return AngleAxisd 角轴
 */
AngleAxisd QuatToAngleAxis(const Quaterniod *q);

/** 
 * @brief 将欧拉角转换为角轴 
 * @param euler 欧拉角结构体指针
 * @return AngleAxisd 角轴
 */
AngleAxisd EulerToAngleAxis(const EulerAngled *euler);

/** 
 * @brief 将角轴转换为欧拉角 
 * @param aa 角轴结构体指针
 * @return EulerAngled 欧拉角
 */
EulerAngled AngleAxisToEuler(const AngleAxisd *aa);


/*******************************************************************************
 * @subsection 旋转的逆运算
*******************************************************************************/
/**
 * @brief 计算四元数的逆
 * @param q 四元数结构体指针
 * @return Quaterniod 逆四元数
 */
Quaterniod InverseQuat(const Quaterniod *q);

/**
 * @brief 计算旋转矩阵的逆
 * @param R 旋转矩阵结构体指针
 * @return RotationMatrixd 逆旋转矩阵
 */
RotationMatrixd InverseRotMat(const RotationMatrixd *R);


/*******************************************************************************
 * @subsection 旋转的乘法运算
*******************************************************************************/
/**
 * @brief 计算两个旋转矩阵的乘积
 * @param R1 旋转矩阵结构体指针
 * @param R2 旋转矩阵结构体指针
 * @return RotationMatrixd 乘积旋转矩阵
 */
RotationMatrixd RotMatMulRotMat(const RotationMatrixd *R1, 
                                const RotationMatrixd *R2);

/**
 * @brief 计算两个四元数的乘积
 * @param q1 四元数结构体指针
 * @param q2 四元数结构体指针
 * @return Quaterniod 乘积四元数
 */
Quaterniod QuatMulQuat(const Quaterniod *q1, const Quaterniod *q2);

#endif // _CMATHUTILS_ROTATION_H_
