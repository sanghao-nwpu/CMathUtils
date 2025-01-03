#include <stdio.h>
#include <windows.h>
#include "cmu_earth.h"

void print_vector3D(const char *label, const Vector3D *vec) {
    printf("%s: (x: %.8f, y: %.8f, z: %.8f)\n", label, vec->x, vec->y, vec->z);
}

int main() {

    // 设置控制台输出为UTF-8编码
    SetConsoleOutputCP(CP_UTF8);

    // 测试参数
    Vector3D ref_blh = {D2R(30.0), D2R(120.0), 0.0};
    Vector3D blh = {D2R(30.001), D2R(120.001), 10.0};
    Vector3D ecef;
    Vector3D enu = {10.0, 20.0, 30.0};
    Vector3D result = { 0.0 };
    
    // 测试 CalculateGravity
    double gravity = CalculateGravity(&blh);
    printf("重力加速度: %f m/s^2\n", gravity);
    
    // 测试 CalculateMeridianRadius
    double meridian_radius = CalculateMeridianRadius(blh.x);
    printf("子午圈半径: %f m\n", meridian_radius);
    
    // 测试 CalculatePrimeVerticalRadius
    double prime_vertical_radius = CalculatePrimeVerticalRadius(blh.x);
    printf("卯酉圈半径: %f m\n", prime_vertical_radius);
    
    // 测试 CalculateEnuToEcefMatrix
    Matrix3D enu_to_ecef_matrix = CalculateEnuToEcefMatrix(&blh);
    printf("ENU到ECEF转换矩阵计算完成。\n");
    
    // 测试 EcefToBlh
    ecef = BlhToEcef(&blh);
    printf("地心地固系坐标转换为纬经高坐标:\n");
    print_vector3D("ECEF坐标", &ecef);
    
    // 测试 BlhToEcef
    Vector3D converted_blh = EcefToBlh(&ecef);
    printf("纬经高坐标转换为地心地固系坐标:\n");
    print_vector3D("BLH坐标", &converted_blh);
    
    // 测试 BlhToEnu
    result = BlhToEnu(&ref_blh, &blh);
    printf("BLH转换为ENU坐标:\n");
    print_vector3D("ENU坐标", &result);
    
    // 测试 EnuToBlh
    result = EnuToBlh(&ref_blh, &enu);
    result.x = R2D(result.x);
    result.y = R2D(result.y);
    printf("ENU坐标转换为BLH坐标:\n");
    print_vector3D("转换后的BLH坐标", &result);

    printf("Press Enter to continue...");
    while(getchar() != '\n');
    return 0;
}
