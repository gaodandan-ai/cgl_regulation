import os
import sys
import math

if not hasattr(math, 'comb'):
    def math_comb(n, k):
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        k = min(k, n - k)
        numerator = 1
        denominator = 1
        for i in range(1, k + 1):
            numerator *= n - i + 1
            denominator *= i
        return numerator // denominator
    math.comb = math_comb

# Add the parent directory to Python path so we can import backend.app
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)

from backend.app import app
