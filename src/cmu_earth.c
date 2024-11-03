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
#include "cmu_matrix.h"
#include <math.h>
#include <stdio.h>

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
 * @param lat 纬度坐标(latitude{rad})
 * @return 子午圈半径(m)
 */
double calculate_meridian_radius(const double lat)
{
    double meridian_radius = 0.0;
    double squre_sin_lat = 0.0;
    squre_sin_lat = sin(lat);
    squre_sin_lat *= squre_sin_lat;
    meridian_radius = WGS84_RA * (1.0 - WGS84_SQURE_E1) / pow(1 - WGS84_SQURE_E1 * squre_sin_lat, 1.5);
    return meridian_radius;
}

/**
 * @brief 计算卯酉圈半径
 * @param lat 纬度坐标(latitude{rad})
 * @return 卯酉圈半径(m)
 */
double calculate_prime_vertical_radius(const double lat)
{
    double prime_vertical_radius = 0.0;
    double sinlat = sin(lat);
    prime_vertical_radius = WGS84_RA / sqrt(1.0 - WGS84_SQURE_E1 * sinlat * sinlat);
    return prime_vertical_radius;
}

/**
 * @brief 计算东北天系(enu)到地心地固系(ecef)的转换矩阵
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 转换矩阵 
 */
Matrix3D calculate_enu_to_ecef_matrix(const Vector3D *blh)
{
    Matrix3D enu_to_ecef_matrix = { 0 };
    double sinlat = sin(blh->x);
    double coslat = cos(blh->x);
    double sinlon = sin(blh->y);
    double coslon = cos(blh->y);

    enu_to_ecef_matrix.data[0][0] = -sinlon;
    enu_to_ecef_matrix.data[0][1] = -coslon * sinlat;
    enu_to_ecef_matrix.data[0][2] = coslat * coslon;

    enu_to_ecef_matrix.data[1][0] = coslon; 
    enu_to_ecef_matrix.data[1][1] = -sinlon * sinlat;
    enu_to_ecef_matrix.data[1][2] = coslat * sinlon;

    enu_to_ecef_matrix.data[2][0] = 0.0;
    enu_to_ecef_matrix.data[2][1] = coslat;
    enu_to_ecef_matrix.data[2][2] = sinlat;

    return enu_to_ecef_matrix;
}

/**
 * @brief 地心地固系(ecef)坐标转换为纬经高坐标(blh)
 * @param ecef 地心地固系坐标(x{m}, y{m}, z{m})
 * @return 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 */
Vector3D ecef_to_blh(const Vector3D *ecef)
{
    Vector3D blh = { 0 };
    double p = sqrt(ecef->x * ecef->x + ecef->y * ecef->y);
    double rn;
    double lat, lon;
    double h = 0.0, h2 = 0.0;

    // 初始状态
    lat = atan2(ecef->z , (p * (1.0 - WGS84_SQURE_E1)));
    lon = atan2(ecef->y, ecef->x);

    while (TRUE) {
        h2 = h;
        rn = calculate_prime_vertical_radius(lat);
        h = p / cos(lat) - rn;
        // lat = atan(ecef->z / (p * (1.0 - WGS84_SQURE_E1 * rn / (rn + h))));
        lat = atan2(ecef->z + rn * WGS84_SQURE_E1 * sin(lat), p);
        if (fabs(h - h2) <= 1.0e-4) {
            break; // 当条件满足时退出循环
        }
    }

    blh.x = lat;
    blh.y = lon;
    blh.z = h;

    return blh;
}

/**
 * @brief 纬经高坐标(blh)转换为地心地固系(ecef)坐标
 * @param blh 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 地心地固系坐标(x{m}, y{m}, z{m})
 */
Vector3D blh_to_ecef(const Vector3D *blh)
{
    Vector3D ecef = { 0 };
    double coslat, sinlat, coslon, sinlon;
    double rnh, rn;

    coslat = cos(blh->x);
    sinlat = sin(blh->x);
    coslon = cos(blh->y);
    sinlon = sin(blh->y);

    rn  = calculate_prime_vertical_radius(blh->x);
    rnh = rn + blh->z;
    ecef.x = rnh * coslat * coslon;
    ecef.y = rnh * coslat * sinlon;
    ecef.z = (rnh - rn * WGS84_SQURE_E1) * sinlat;
    return ecef;
}

/**
 * @brief 纬经高坐标(blh)转换为东北天系(enu)坐标
 * @param ref_blh 参考点的纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @param blh 待转换点的纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @return 东北天系坐标(east{m}, north{m}, up{m})
 */
Vector3D blh_to_enu(const Vector3D *ref_blh, const Vector3D *blh)
{
    Vector3D enu = { 0 };
    Vector3D ref_ecef = { 0 };
    Vector3D ecef = { 0 };
    Matrix3D enu_to_ecef_matrix = { 0 };
    Matrix3D ecef_to_enu_matrix = { 0 };
    Vector3D delta_ecef = { 0 };

    ref_ecef = blh_to_ecef(ref_blh);
    ecef = blh_to_ecef(blh);
    enu_to_ecef_matrix = calculate_enu_to_ecef_matrix(ref_blh);
    ecef_to_enu_matrix = MatTranspose3D(&enu_to_ecef_matrix);

    delta_ecef = VecSub3D(&ecef, &ref_ecef);
    enu = MatMulVec3D(&ecef_to_enu_matrix, &delta_ecef);

    return enu;
}

/**
 * @brief 东北天系(enu)坐标转换为纬经高坐标(blh)
 * @param ref_blh 参考点的纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 * @param enu 待转换点的东北天系坐标(east{m}, north{m}, up{m})
 * @return 纬经高坐标(latitude{rad}, longitude{rad}, height{m})
 */
Vector3D enu_to_blh(const Vector3D *ref_blh, const Vector3D *enu)
{
    Vector3D blh = { 0 };
    Vector3D ref_ecef = { 0 };
    Vector3D ecef = { 0 };
    Matrix3D enu_to_ecef_matrix = { 0 };
    Matrix3D ecef_to_enu_matrix = { 0 };
    Vector3D delta_ecef = { 0 };

    ref_ecef = blh_to_ecef(ref_blh);
    enu_to_ecef_matrix = calculate_enu_to_ecef_matrix(ref_blh);

    delta_ecef = MatMulVec3D(&enu_to_ecef_matrix, enu);
    ecef = VecAdd3D(&ref_ecef, &delta_ecef);

    blh = ecef_to_blh(&ecef);

    return blh;
}