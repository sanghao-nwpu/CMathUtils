# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 22:20:51 2024

@author: sanghao
"""

import numpy as np


# 示例
A = np.array([[1, 0, 1],
              [0, 2, 0],
              [1, 0, 3]], dtype=float)

B = np.array([[1, 2, 0, 0],
              [0, 2, 0, 0],
              [0, 0, 3, 0],
              [2, 0, 0, 4]], dtype=float)

print(np.linalg.inv(A))
print(np.linalg.inv(B))