#include <stdio.h>
#include "cmu_rotation.h"

void testRotation() {
    EulerAngled euler = {M_PI / 4, M_PI / 4, M_PI / 4};  // 45 degrees around each axis
    AngleAxisd angle_axis = {1 / (M_PI / 2), 0, 0}; // 90 degrees around x-axis
    Quaterniod quat = {0.7071, 0.7071, 0, 0}; // Example quaternion

    // 测试欧拉角到旋转矩阵
    RotationMatrixd R1 = EulerToRotMat(&euler);
    printf("EulerToRotMat:\n");
    printf("[[%f, %f, %f], [%f, %f, %f], [%f, %f, %f]]\n",
           R1.data[0][0], R1.data[0][1], R1.data[0][2],
           R1.data[1][0], R1.data[1][1], R1.data[1][2],
           R1.data[2][0], R1.data[2][1], R1.data[2][2]);

    // 测试旋转矩阵到欧拉角
    EulerAngled eulerBack = RotMatToEuler(&R1);
    printf("RotMatToEuler: roll = %f, pitch = %f, yaw = %f\n", eulerBack.roll, eulerBack.pitch, eulerBack.yaw);

    // 测试角轴到旋转矩阵
    RotationMatrixd R2 = AngleAxisToMat(&angle_axis);
    printf("AngleAxisToMat:\n");
    printf("[[%f, %f, %f], [%f, %f, %f], [%f, %f, %f]]\n",
           R2.data[0][0], R2.data[0][1], R2.data[0][2],
           R2.data[1][0], R2.data[1][1], R2.data[1][2],
           R2.data[2][0], R2.data[2][1], R2.data[2][2]);

    // 测试旋转矩阵到角轴
    AngleAxisd axisBack = RotMatToAngleAxis(&R2);
    printf("RotMatToAngleAxis: axis = [%f, %f, %f], angle = %f\n", axisBack.axis[0], axisBack.axis[1], axisBack.axis[2], axisBack.angle);

    // 测试四元数到旋转矩阵
    RotationMatrixd R3 = QuatToRotMat(&quat);
    printf("QuatToRotMat:\n");
    printf("[[%f, %f, %f], [%f, %f, %f], [%f, %f, %f]]\n",
           R3.data[0][0], R3.data[0][1], R3.data[0][2],
           R3.data[1][0], R3.data[1][1], R3.data[1][2],
           R3.data[2][0], R3.data[2][1], R3.data[2][2]);

    // 测试旋转矩阵到四元数
    Quaterniod quatBack = RotMatToQuat(&R3);
    printf("RotMatToQuat: w = %f, x = %f, y = %f, z = %f\n", quatBack.w, quatBack.x, quatBack.y, quatBack.z);

    // 测试欧拉角到四元数
    Quaterniod quatFromEuler = EulerToQuat(&euler);
    printf("EulerToQuat: w = %f, x = %f, y = %f, z = %f\n", quatFromEuler.w, quatFromEuler.x, quatFromEuler.y, quatFromEuler.z);

    // 测试四元数到欧拉角
    EulerAngled eulerFromQuat = QuatToEuler(&quat);
    printf("QuatToEuler: roll = %f, pitch = %f, yaw = %f\n", eulerFromQuat.roll, eulerFromQuat.pitch, eulerFromQuat.yaw);
}

int main() {
    testRotation();

    printf("Press Enter to continue...");
    while(getchar() != '\n');

    return 0;
}
