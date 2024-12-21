#include <math.h>
#include "cmu_rotation.h"

/**
 * @brief 将欧拉角转换为旋转矩阵
 * @param euler 欧拉角结构体指针
 * @return RotationMatrixd 旋转矩阵
 */
RotationMatrixd EulerToRotMat(const EulerAngled *euler)
{
    RotationMatrixd R;
    double cr = cos(euler->roll);
    double sr = sin(euler->roll);
    double cp = cos(euler->pitch);
    double sp = sin(euler->pitch);
    double cy = cos(euler->yaw);
    double sy = sin(euler->yaw);

    R.data[0][0] = cp * cy;
    R.data[0][1] = cy * sp * sr - cr * sy;
    R.data[0][2] = sr * sy + cr * sp * cy;

    R.data[1][0] = cp * sy;
    R.data[1][1] = cr * cy + sp * sr * sy;
    R.data[1][2] = cr * sp * sy - sr * cy;

    R.data[2][0] = -sp;
    R.data[2][1] = cp * sr;
    R.data[2][2] = cr * cp;

    return R;
}

/**
 * @brief 将旋转矩阵转换为欧拉角
 * @param R 旋转矩阵结构体指针
 * @return EulerAngled 欧拉角
 * @note 这里的欧拉角顺序为roll-pitch-yaw
 * @note 在万向锁情况下，只能计算roll和pitch的和或差，此时定义roll为0
 */
EulerAngled RotMatToEuler(const RotationMatrixd *R)
{
    EulerAngled euler;
    euler.pitch = asin(-R->data[2][0]);

    if (cos(euler.pitch) > 1e-6)
    {
        euler.yaw = atan2(R->data[1][0], R->data[0][0]);
        euler.roll = atan2(R->data[2][1], R->data[2][2]);
    }
    else
    {
        euler.yaw = atan2(-R->data[0][1], R->data[1][1]);
        euler.roll = 0;
    }

    return euler;
}

/**
 * @brief 将角轴转换为旋转矩阵
 * @param aa 角轴结构体指针
 * @return RotationMatrixd 旋转矩阵
 */
RotationMatrixd AngleAxisToRotMat(const AngleAxisd *aa)
{
    RotationMatrixd R;
    double angle = aa->angle;
    double x = aa->axis[0];
    double y = aa->axis[1];
    double z = aa->axis[2];

    double c = cos(angle);
    double s = sin(angle);
    double t = 1 - c;

    R.data[0][0] = t * x * x + c;
    R.data[0][1] = t * x * y - z * s;
    R.data[0][2] = t * x * z + y * s;

    R.data[1][0] = t * x * y + z * s;
    R.data[1][1] = t * y * y + c;
    R.data[1][2] = t * y * z - x * s;

    R.data[2][0] = t * x * z - y * s;
    R.data[2][1] = t * y * z + x * s;
    R.data[2][2] = t * z * z + c;

    return R;
}

/**
 * @brief 将旋转矩阵转换为角轴
 * @param R 旋转矩阵结构体指针
 * @return AngleAxisd 角轴
 * @note 通过quat作为中间量来计算
 */
AngleAxisd RotMatToAngleAxis(const RotationMatrixd *R)
{
    AngleAxisd aa;
    Quaterniod quat;

    quat = RotMatToQuat(R);
    aa = QuatToAngleAxis(&quat);

    return aa;
}

/**
 * @brief 将四元数转换为旋转矩阵
 * @param q 四元数结构体指针
 * @return RotationMatrixd 旋转矩阵
 */
RotationMatrixd QuatToRotMat(const Quaterniod *q)
{
    RotationMatrixd R;
    double w = q->w, x = q->x, y = q->y, z = q->z;

    R.data[0][0] = 1 - 2 * (y * y + z * z);
    R.data[0][1] = 2 * (x * y - z * w);
    R.data[0][2] = 2 * (x * z + y * w);

    R.data[1][0] = 2 * (x * y + z * w);
    R.data[1][1] = 1 - 2 * (x * x + z * z);
    R.data[1][2] = 2 * (y * z - x * w);

    R.data[2][0] = 2 * (x * z - y * w);
    R.data[2][1] = 2 * (y * z + x * w);
    R.data[2][2] = 1 - 2 * (x * x + y * y);

    return R;
}

/**
 * @brief 将旋转矩阵转换为四元数
 * @param R 旋转矩阵结构体指针
 * @return Quaterniod 四元数
 */
Quaterniod RotMatToQuat(const RotationMatrixd *R)
{
    Quaterniod q;
    double tr = R->data[0][0] + R->data[1][1] + R->data[2][2];
    if (tr > 0)
    {
        double S = sqrt(tr + 1.0) * 2; // S = 4 * q.w
        q.w = 0.25 * S;
        q.x = (R->data[2][1] - R->data[1][2]) / S;
        q.y = (R->data[0][2] - R->data[2][0]) / S;
        q.z = (R->data[1][0] - R->data[0][1]) / S;
    }
    else if ((R->data[0][0] > R->data[1][1]) && (R->data[0][0] > R->data[2][2]))
    {
        double S = sqrt(1.0 + R->data[0][0] - R->data[1][1] - R->data[2][2]) * 2; // S = 4 * q.x
        q.w = (R->data[2][1] - R->data[1][2]) / S;
        q.x = 0.25 * S;
        q.y = (R->data[0][1] + R->data[1][0]) / S;
        q.z = (R->data[0][2] + R->data[2][0]) / S;
    }
    else if (R->data[1][1] > R->data[2][2])
    {
        double S = sqrt(1.0 + R->data[1][1] - R->data[0][0] - R->data[2][2]) * 2; // S = 4 * q.y
        q.w = (R->data[0][2] - R->data[2][0]) / S;
        q.x = (R->data[0][1] + R->data[1][0]) / S;
        q.y = 0.25 * S;
        q.z = (R->data[1][2] + R->data[2][1]) / S;
    }
    else
    {
        double S = sqrt(1.0 + R->data[2][2] - R->data[0][0] - R->data[1][1]) * 2; // S = 4 * q.z
        q.w = (R->data[1][0] - R->data[0][1]) / S;
        q.x = (R->data[0][2] + R->data[2][0]) / S;
        q.y = (R->data[1][2] + R->data[2][1]) / S;
        q.z = 0.25 * S;
    }

    return q;
}

/**
 * @brief 将欧拉角转换为四元数
 * @param euler 欧拉角结构体指针
 * @return Quaterniod 四元数
 */
Quaterniod EulerToQuat(const EulerAngled *euler)
{
    Quaterniod q;

    double cr = cos(euler->roll * 0.5);
    double sr = sin(euler->roll * 0.5);
    double cp = cos(euler->pitch * 0.5);
    double sp = sin(euler->pitch * 0.5);
    double cy = cos(euler->yaw * 0.5);
    double sy = sin(euler->yaw * 0.5);

    q.w = cr * cp * cy + sr * sp * sy;
    q.x = sr * cp * cy - cr * sp * sy;
    q.y = cr * sp * cy + sr * cp * sy;
    q.z = cr * cp * sy - sr * sp * cy;

    return q;
}

/**
 * @brief 将四元数转换为欧拉角
 * @param q 四元数结构体指针
 * @return EulerAngled 欧拉角
 */
EulerAngled QuatToEuler(const Quaterniod *q)
{
    EulerAngled euler;
    double sinr_cosp = 2 * (q->w * q->x + q->y * q->z);
    double cosr_cosp = 1 - 2 * (q->x * q->x + q->y * q->y);
    euler.roll = atan2(sinr_cosp, cosr_cosp);

    double sinp = 2 * (q->w * q->y - q->z * q->x);
    if (fabs(sinp) >= 1)
        euler.pitch = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    else
        euler.pitch = asin(sinp);

    double siny_cosp = 2 * (q->w * q->z + q->x * q->y);
    double cosy_cosp = 1 - 2 * (q->y * q->y + q->z * q->z);
    euler.yaw = atan2(siny_cosp, cosy_cosp);

    return euler;
}

/**
 * @brief 将角轴转换为四元数
 * @param aa 角轴结构体指针
 * @return Quaterniod 四元数
 */
Quaterniod AngleAxisToQuat(const AngleAxisd *aa)
{
    Quaterniod q;
    double halfAngle = aa->angle * 0.5;
    double s = sin(halfAngle);

    q.w = cos(halfAngle);
    q.x = aa->axis[0] * s;
    q.y = aa->axis[1] * s;
    q.z = aa->axis[2] * s;

    return q;
}

/**
 * @brief 将四元数转换为角轴
 * @param q 四元数结构体指针
 * @return AngleAxisd 角轴
 */
AngleAxisd QuatToAngleAxis(const Quaterniod *q)
{
    AngleAxisd aa;
    aa.angle = 2 * acos(q->w);
    double s = sqrt(1 - q->w * q->w);

    if (s < 1e-6)
    {
        aa.axis[0] = 0;
        aa.axis[1] = 0;
        aa.axis[2] = 0;
    }
    else
    {
        aa.axis[0] = q->x / s;
        aa.axis[1] = q->y / s;
        aa.axis[2] = q->z / s;
    }

    return aa;
}

/**
 * @brief 将欧拉角转换为角轴
 * @param euler 欧拉角结构体指针
 * @return AngleAxisd 角轴
 * @note 通过旋转矩阵作为中间量来计算
 */
AngleAxisd EulerToAngleAxis(const EulerAngled *euler)
{
    AngleAxisd aa;
    RotationMatrixd R = {0};
    R = EulerToRotMat(euler);
    aa = RotMatToAngleAxis(&R);

    return aa;
}

/**
 * @brief 将角轴转换为欧拉角
 * @param aa 角轴结构体指针
 * @return EulerAngled 欧拉角
 * @note 通过旋转矩阵作为中间量来计算
 */
EulerAngled AngleAxisToEuler(const AngleAxisd *aa)
{
    // 这里需要一个根据实际角轴计算欧拉角的实现
    // 示例返回值
    EulerAngled euler;
    RotationMatrixd R = {0};
    R = AngleAxisToRotMat(aa);
    euler = RotMatToEuler(&R);

    return euler;
}

/**
 * @brief 计算四元数的逆
 * @param q 四元数结构体指针
 * @return Quaterniod 逆四元数
 */
Quaterniod InverseQuat(const Quaterniod *q)
{
    Quaterniod inverse;
    double norm = q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z;

    inverse.w = q->w / norm;
    inverse.x = -q->x / norm;
    inverse.y = -q->y / norm;
    inverse.z = -q->z / norm;

    return inverse;
}

/**
 * @brief 计算旋转矩阵的逆
 * @param R 旋转矩阵结构体指针
 * @return RotationMatrixd 逆旋转矩阵
 */
RotationMatrixd InverseRotMat(const RotationMatrixd *R)
{
    RotationMatrixd inverse;
    // 旋转矩阵的逆等于其转置
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            inverse.data[j][i] = R->data[i][j];
        }
    }
    return inverse;
}

/**
 * @brief 计算两个旋转矩阵的乘积
 * @param R1 旋转矩阵结构体指针
 * @param R2 旋转矩阵结构体指针
 * @return RotationMatrixd 乘积旋转矩阵
 */
RotationMatrixd RotMatMulRotMat(const RotationMatrixd *R1, const RotationMatrixd *R2)
{
    RotationMatrixd result;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result.data[i][j] = 0;
            for (int k = 0; k < 3; k++)
            {
                result.data[i][j] += R1->data[i][k] * R2->data[k][j];
            }
        }
    }
    return result;
}

/**
 * @brief 计算两个四元数的乘积
 * @param q1 四元数结构体指针
 * @param q2 四元数结构体指针
 * @return Quaterniod 乘积四元数
 */
Quaterniod QuatMulQuat(const Quaterniod *q1, const Quaterniod *q2)
{
    Quaterniod result;
    result.w = q1->w * q2->w - q1->x * q2->x - q1->y * q2->y - q1->z * q2->z;
    result.x = q1->w * q2->x + q1->x * q2->w + q1->y * q2->z - q1->z * q2->y;
    result.y = q1->w * q2->y - q1->x * q2->z + q1->y * q2->w + q1->z * q2->x;
    result.z = q1->w * q2->z + q1->x * q2->y - q1->y * q2->x + q1->z * q2->w;

    return result;
}
