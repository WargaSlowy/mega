# MIT License

# Copyright (c) 2025 WargaSlowy

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# cython: language_level=3
# mega/utils/constant.pyx

cdef double PI_NUMBER = 3.141592653589793238462643383279502884197
cdef double EULER_NUMBER = 2.718281828459045235360287471352
cdef double GAUSS_NUMBER = 0.8346268
cdef double SQRT_PI = 1.772453850905516027298167483341145182

cpdef unsigned long lucas_number(int n) except 0:
    """
    compute the nth lucas number using iterative

    Parameter:
        n (int): index in the lucas squence, must be non-negative integers

    Return:
        (unsigned long): the nth lucas number

    Example:
    >>> lucas_number(0)
    2
    >>> lucas_number(40)
    228826127
    """
    if n < 0:
        raise ValueError("index must be non negative number")

    if n == 0:
        return 2
    elif n == 1:
        return 1

    cdef unsigned long prev = 2, curr = 1, next_val
    cdef int i
    
    for i in range(2, n + 1):
        next_val = prev + curr
        # shifting previous value forward
        prev = curr
        # update current value for next iteration
        curr = next_val
    return curr
