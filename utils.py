import numpy as np


def GetPrevAndNext(x, array):
    cur_idx = array.index(x)
    assert cur_idx != 0 and cur_idx != len(array) - 1
    prev = array[cur_idx - 1]
    next = array[cur_idx + 1]
    return prev, next


def NormalToPlane(x1, x2, x3):
    v1 = x2 - x1
    v2 = x3 - x1
    return np.cross(v1, v2)


def Normalized(v: np.ndarray):
    return v / np.linalg.norm(v)
