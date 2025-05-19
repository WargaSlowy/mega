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
# mega/op/gamma.pyx

from mega.utils.constant cimport PI_NUMBER, SQRT_PI
from libc.math cimport sqrt, exp, log, sin, pow

cpdef double gamma(double point):
    """
    computing gamma function Γ(point) for real-valued `point > 0`

    implement using lanczos approximation which providing an efficient and
    accurate way to compute Γ(z) for non-integer values

    Parameter:
        point (double): real number greater than zero at which to evaluate  gamma function

    Return:
        (double): value of the gamma function at `gamma`

    Note:
        - Γ(0.5) = √π ≈ 1.77245385091
        - Γ(1.0) = 1
        - Γ(2.0) = 1

    Example:
    >>> gamma(0.5)
    1.77245385091
    
    >>> gamma(5.0)
    24.0

    Reference:
     - https://en.wikipedia.org/wiki/Gamma_function
     - https://en.wikipedia.org/wiki/Lanczos_approximation
    """
    cdef double z = point
    cdef double x, y, t, tmp, fact
    cdef int i, n
    cdef int neg = 0

    if point <= 0:
        raise ValueError("point must be bigger than zero")
    if point > 175.5:
        raise OverflowError("input value too large, will be OverflowError")

    # special case
    if point == 0.5:
        return SQRT_PI

    if point == 1.0 or point == 2.0:
        return 1.0

    # factorial shorcut for exact integer input
    if z == <int>z:
        n  = <int>z;
        fact = 1.0
        for i in range(2, n):
            fact *= i
        return fact

    # duplicate check to avoid confusion
    if point == int(point):
        n = <int>point
        result = 1.0
        for i in range(2, n):
            result *= i
        return result

    # reset working variable
    z = point
    y = z
    
    # apply reflection formula
    # this allow to compute (for small value)
    neg = 0
    if z < 0.5:
        y = 1.0 - z
        neg = 1
    else:
        y = z

    x = 0.99999999999980993

    # coefficient from lanczos approximation with g=7 and n=8 terms
    # theres are precomputed to balance accuracy and performance
    cdef double[8] LANCZOS_COEFF = [
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    ]
    
    # accumulate terms of the series
    for i in range(8):
        x += LANCZOS_COEFF[i] / (y + i)
    # compute intermediate value for the final formula
    t = y + 7.5
    tmp = sqrt(2.0 * PI_NUMBER) * x * pow(t, y - 0.5) * exp(-t)
    # apply reflection formula if needed
    if neg:
        return PI_NUMBER / (sin(PI_NUMBER * z) * tmp)
    return tmp

cpdef int euler_phi(int n) except -1:
    """
    compute euler totients function which couynt the number of integers
    less than or equal to `n` that are coprime to `n`

    formula:
    φ(n) = n × ∏(p|n) (1 - 1/p)

    where the product is over all distinct prime factors p of n

    Parameter:
        n (int): positive integer greater than 0

    Return:
        (int): value of euler totient function

    Example:
    >>> euler_phi(10)
    4
    >>> euler_phi(100)
    40
    """
    if n <= 0:
        raise ValueError("only positive number will accept")

    cdef int result = n
    cdef int i = 2 # start check from smallest prime factor

    if n % 2 == 0:
        # apply the formula
        # result = result * (1 - 1 / 2) = result - result // 2
        result -= result // 2
        # remove all occurrences of 2 from n
        while n % 2 == 0:
            n //= 2

    i = 3
    while i * i <= n:
        if n % i == 0:
            # apply the formula
            # result = result * (1 - 1 / i) = result - result // i
            result -= result // i
            # remove all occurence of current prime factor i
            while n % i == 0:
                n //= i
        i += 2

    # if remaining n is a prime > 2, apply last adjustment
    if n > 1:
        result -= result // n

    return result
