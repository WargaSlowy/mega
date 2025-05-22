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

cdef void mat_mult(unsigned long A[2][2], unsigned long B[2][2], unsigned long result[2][2]) except *:
    """
    multiply 2x2 matrices A B, and storing in result variable

    formula:
        [ a00 a01 ]   [ b00 b01 ]   [ a00*b00 + a01*b10   a00*b01 + a01*b11 ]
        [ a10 a11 ] × [ b10 b11 ] = [ a10*b00 + a11*b10   a10*b01 + a11*b11 ]

    Parameter:
        A (unsigned long[2][2]): first input 2x2 matrices
        B (unsigned long[2][2]): second input 2x2 matrices
        result (unsigned long[2][2]): output matrices to storing the product
    """
    # top left element result matrices
    result[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0]
    # top right element result matrices
    result[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1]
    # bottom left element of result matrices
    result[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0]
    # bottom right element of result matrices
    result[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1]

cdef void mat_pow(unsigned long M[2][2], int power, unsigned long result[2][2]) except *:
    """
    raise 2x2 matrice to given power using binary exponent
    
    Parameter:
        M (unsigned long[2][2]): input 2x2 matrices to be raise to power
        power (int): non-negative integer exponent
        result (unsigned long[2][2]): output matrix where the result will stored
    """
    result[0][0] = 1
    result[0][1] = 0
    result[1][0] = 0
    result[1][1] = 1

    cdef unsigned long temp[2][2]

    while power > 0:
        if power % 2 == 1:
            mat_mult(result, M, temp)
            # inline copy back to the res to avoid extra function call
            result[0][0] = temp[0][0]
            result[0][1] = temp[0][1]
            result[1][0] = temp[1][0]
            result[1][1] = temp[1][1]
        # square current matrix M
        mat_mult(M, M, temp)
        M[0][0] = temp[0][0]
        M[0][1] = temp[0][1]
        M[1][0] = temp[1][0]
        M[1][1] = temp[1][1]
        # halve power
        power //= 2


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
    cdef unsigned long prev = 2, curr = 1, next_val
    cdef int i
    cdef unsigned long res[2][2]
    cdef int power
    cdef unsigned long F[2][2]
    F[0][0] = 1
    F[0][1] = 1
    F[1][0] = 1
    F[1][1] = 0

    if n < 0:
        raise ValueError("index must be non negative number")

    if n == 0:
        return 2
    elif n == 1:
        return 1

    if n <= 20:
        for i in range(2, n + 1):
            next_val = prev + curr
            # shifting previous value forward
            prev = curr
            # update current value for next iteration
            curr = next_val
        return curr

    else:        
        power = n - 1
        mat_pow(F, power, res)
        return 2 * res[0][1] + res[0][0]
