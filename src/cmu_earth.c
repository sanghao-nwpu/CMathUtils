/**
 * @file cmu_earth.c
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

#include "cmu_earth.h"
#include <math.h>

/**
 * @brief 计算地球上某点的重力加速度
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 重力加速度(m/s^2)
 */
double calculate_gravity(const Vector3D *blh)
{
    double squre_sin_lat = 0.0;
    double gravity = 0.0;
    squre_sin_lat = sin(blh->x);
    squre_sin_lat *= squre_sin_lat;
    /** g = 9.7803267715 * (1 + 0.0052790414 * sin^2(lat) + 0.0000232718 * sin^2(lat)) + 
     * blh.z * (0.0000000043977311 * sin^2(lat) - 0.0000030876910891) 
     * + 0.0000000000007211 * z^2 */
    gravity = 1 + 0.0052790414 * squre_sin_lat;
    gravity += 0.0000232718 * squre_sin_lat * squre_sin_lat;
    gravity = 9.7803267715 * gravity;
    gravity += blh->z * (0.0000000043977311 * squre_sin_lat - 0.0000030876910891);
    gravity += 0.0000000000007211 * blh->z * blh->z;
    return gravity;
}


/**
 * @brief 计算子午圈半径
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 子午圈半径(m)
 */
double calculate_meridian_radius(const Vector3D *blh)
{
    double meridian_radius = 0.0;
}