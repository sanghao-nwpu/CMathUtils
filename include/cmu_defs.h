#ifndef _CMATHUTILS_DEFINE_H_
#define _CMATHUTILS_DEFINE_H_

/**
 * @file cmu_defs.h
 * @brief CMathUtils的基本数据类型定义
 *
 * 本文件定义了CMathUtils的基本数据类型定义，包括向量、矩阵、旋转的四种定义。
 *
 * @author [sanghao_nwpu]
 * @version 1.0.0
 * @date 2024年10月28日
 *
 * @section change_log 版本变更说明
 * - 1.0.0: 初始版本。
 */

/*******************************************************************************
 * @subsection 基本数学定义
 *******************************************************************************/
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define D2R(x) ((x) * 1.7453292519943295769e-2)
#define R2D(x) ((x) * 5.7295779513082320877e+1)

#define EPSILON (1e-6)

/*******************************************************************************
 * @subsection WGS84椭球模型参数
 * @note 如果使用其他椭球模型需要修改椭球参数
 *******************************************************************************/
#define WGS84_WIE (7.2921151467E-5)            /** 地球自转角速度*/
#define WGS84_F (0.0033528106647474805)        /** 地球扁率 */
#define WGS84_RA (6378137.0000000000)          /** 长半轴a */
#define WGS84_RB (6356752.3142451793)          /** 短半轴b */
#define WGS84_GM0 (398600441800000.00)         /** 地球引力常数 */
#define WGS84_SQURE_E1 (0.0066943799901413156) /** 第一偏心率平方 */
#define WGS84_SQURE_E2 (0.0067394967422764341) /** 第二偏心率平方 */

/*******************************************************************************
 * @subsection 状态类型定义
 *******************************************************************************/
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

#define CMU_STATUS_SUCCESS (0)
#define CMU_STATUS_FAILURE (-1)
#define CMU_ERROR_MATRIX_OPERATION_INVALID (-2)

typedef int Status;

/**
 * @brief 状态枚举类型
 * @param STATUS_SUCCESS 成功
 * @param STATUS_FAILURE 失败
 */
typedef enum
{
    STATUS_SUCCESS = 0,
    STATUS_FAILURE = 1
} StatusEnum;

/*******************************************************************************
 * @subsection 向量和矩阵类型定义
 *******************************************************************************/
/**
 * @brief 2维向量定义(联合体)
 * @param data[2] 存储向量的数组
 * @param x 向量的 x 分量
 * @param y 向量的 y 分量
 */
typedef union
{
    double data[2];
    struct
    {
        double x;
        double y;
    };
} Vector2D;

/**
 * @brief 3维向量定义(联合体)
 * @param data[3] 存储向量的数组
 * @param x 向量的 x 分量
 * @param y 向量的 y 分量
 * @param z 向量的 z 分量
 */
typedef union
{
    double data[3];
    struct
    {
        double x;
        double y;
        double z;
    };
} Vector3D;

/**
 * @brief 4维向量定义(联合体)
 * @param data[4] 存储向量的数组
 * @param w 向量的 w 分量(no.1)
 * @param x 向量的 x 分量(no.2)
 * @param y 向量的 y 分量(no.3)
 * @param z 向量的 z 分量(no.4)
 */
typedef union
{
    double data[4];
    struct
    {
        double w;
        double x;
        double y;
        double z;
    };
} Vector4D;

/**
 * @brief 2维矩阵定义
 * @param data[2][2] 存储矩阵的数组
 */
typedef struct
{
    double data[2][2];
} Matrix2D;

/**
 * @brief 3维矩阵定义
 * @param data[3][3] 存储矩阵的数组
 */
typedef struct
{
    double data[3][3];
} Matrix3D;

/**
 * @brief 4维矩阵定义
 * @param data[4][4] 存储矩阵的数组
 */
typedef struct
{
    double data[4][4];
} Matrix4D;

/**
 * @brief 动态向量数据结构(结构体)
 * @param rows 矩阵的行数
 * @param cols 矩阵的列数
 * @param data 二维数组，用于存储矩阵数据
 * @note rows为0表示矩阵并未初始化，data指针不可用
 */
typedef struct
{
    int rows;
    int cols;
    double **data;
} MatrixXD;

/**
 * @brief 动态向量数据结构(结构体)
 * @param size 向量的大小
 * @param data 一维数组，用于存储向量数据
 * @note size为0表示内存并未初始化，data指针不可用
 */
typedef struct
{
    int size;
    double *data;
} VectorXD;

/*******************************************************************************
 * @subsection 旋转的四种定义：欧拉角、角轴、四元数、旋转矩阵
 *******************************************************************************/
/**
 * @brief 欧拉角数据结构(联合体)
 * @param roll 滚转角，绕 x 轴旋转角度 (rad)
 * @param pitch 俯仰角，绕 y 轴旋转角度 (rad)
 * @param yaw 偏航角，绕 z 轴旋转角度 (rad)
 * @param angles[3] 使用数组表示角度(rad)
 */
typedef union
{
    struct
    {
        double roll;
        double pitch;
        double yaw;
    };
    double angles[3];
} EulerAngled;

/**
 * @brief 角轴数据结构(联合体)
 * @param axis 旋转轴的单位向量
 * @param angle 旋转角度 (rad)
 * @param ax 旋转轴的 x 分量
 * @param ay 旋转轴的 y 分量
 * @param az 旋转轴的 z 分量
 * @param theta 作为单独角度访问 (rad)
 */
typedef union
{
    struct
    {
        double axis[3];
        double angle;
    };
    struct
    {
        double ax;
        double ay;
        double az;
        double theta;
    };
} AngleAxisd;

/**
 * @brief 四元数数据结构(联合体)
 * @param w 标量部分
 * @param x x 分量
 * @param y y 分量
 * @param z z 分量
 * @param data[4] 使用数组表示四元数的四个分量
 */
typedef union
{
    struct
    {
        double w;
        double x;
        double y;
        double z;
    };
    double data[4];
} Quaterniod;

/**
 * @brief 旋转矩阵数据结构(联合体)
 * @param data[3][3] 3x3 旋转矩阵
 * @param flat[9] 使用一维数组访问旋转矩阵
 */
typedef union
{
    struct
    {
        double data[3][3];
    };
    double flat[9];
} RotationMatrixd;

/**
 * @brief 位姿矩阵
 * @param data 4*4矩阵, [R, t; 0, 1]
 */
typedef struct
{
    double data[4][4];
} Pose;

#endif // _CMATHUTILS_DEFINE_H_