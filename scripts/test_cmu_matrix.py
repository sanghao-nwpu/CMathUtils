# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 22:20:51 2024

@author: sanghao
"""

import numpy as np
from scipy.linalg import lu

# 测试矩阵
test_matrix = np.array([[2, 1, 3],
                         [1, 2, 1],
                         [3, 0, 2]])

# 使用 scipy 的 LU 分解
P, L, U = lu(test_matrix)

# 输出结果
print("Test Matrix:")
print(test_matrix)

print("\nP Matrix (Permutation Matrix):")
print(P)

print("\nL Matrix (Lower Triangular Matrix):")
print(L)

print("\nU Matrix (Upper Triangular Matrix):")
print(U)

# 使用 NumPy 求逆
inverse_matrix = np.linalg.inv(test_matrix)

# 输出逆矩阵
print("\nInverse Matrix:")
print(inverse_matrix)

# 验证 A * A_inv = I
identity_check = np.dot(test_matrix, inverse_matrix)
print("\nIdentity Check (A * A_inv):")
print(identity_check)
