#ifndef _CMUCMATHUTILS_EARTH_H_
#define _CMUCMATHUTILS_EARTH_H_

/**
 * @file cmu_earth.h
 * @brief 地球相关运算
 *
 * 本文件定义了地球计算的相关函数接口，包括经纬高和东北天坐标系的互转。
 *
 * @author [sanghao_nwpu]
 * @version 1.0.0
 * @date 2024年10月28日
 *
 * @section change_log 版本变更说明
 * - 1.0.0: 初始版本，添加了基本的向量和矩阵数据结构及运算接口。
 */

#include "cmu_defs.h"

/*******************************************************************************
 * @subsection 相关函数声明
*******************************************************************************/
/**
 * @brief 计算地球上某点的重力加速度
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 重力加速度(m/s^2)
 */
double calculate_gravity(const Vector3D *blh);

/**
 * @brief 计算子午圈半径
 * @param lat 纬度坐标(latitude{rad})
 * @return 子午圈半径(m)
 */
double calculate_meridian_radius(const double lat);

/**
 * @brief 计算卯酉圈半径
 * @param lat 纬度坐标(latitude{rad})
 * @return 卯酉圈半径(m)
 */
double calculate_prime_vertical_radius(const double lat);

/**
 * @brief 计算东北天系(enu)到地心地固系(ecef)的转换矩阵
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 转换矩阵 
 */
Matrix3D calculate_enu_to_ecef_matrix(const Vector3D *blh);

/**
 * @brief 地心地固系(ecef)坐标转换为纬经高坐标(blh)
 * @param ecef 地心地固系坐标(x{m}, y{m}, z{m})
 * @return 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 */
Vector3D ecef_to_blh(const Vector3D *ecef);

/**
 * @brief 纬经高坐标(blh)转换为地心地固系(ecef)坐标
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 地心地固系坐标(x{m}, y{m}, z{m})
 */
Vector3D blh_to_ecef(const Vector3D *blh);

/**
 * @brief 纬经高坐标(blh)转换为东北天系(enu)坐标
 * @param ref_blh 参考点的纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @param blh 待转换点的纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 东北天系坐标(east{m}, north{m}, up{m})
 */
Vector3D blh_to_enu(const Vector3D *ref_blh, const Vector3D *blh);

/**
 * @brief 东北天系(enu)坐标转换为纬经高坐标(blh)
 * @param ref_blh 参考点的纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @param enu 待转换点的东北天系坐标(east{m}, north{m}, up{m})
 * @return 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 */
Vector3D enu_to_blh(const Vector3D *ref_blh, const Vector3D *enu);


#endif // _CMUCMATHUTILS_EARTH_H_