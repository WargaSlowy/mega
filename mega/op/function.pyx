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
from libc.math cimport sqrt, exp, sin, pow, cos

def prime_factors(int n, bint unique=False):
    """"
    return prime factor of positive integer n

    Parameter:
        n (int): integer to be factorized, must greater than zero
        unique (bool): if True, return only distinct prime factors
        
    Return:
        list[int]: list of prime factor of `n`, sorted in increasing order

    Example:
    >>> prime_factors(60)
    [2, 2, 3, 5]
    >>> prime_factors(13, unique=True)
    [13]
    """
    if n <= 0:
        raise ValueError("only positive integer are acc")
    cdef int i = 2
    cdef list factors = []
    
    # trie division algorithm
    while i * i <= n:
        while n % i == 0:
            # append current factor
            factors.append(i)
            # divide n by the factor
            n //= i
        # move to next candidate divisor
        i += 1

    # remaining prime factors
    if n > 1:
        factors.append(n)
    # remove duplicate if unique flag is set
    if unique:
        seen = set()
        # using list comprehension to preserve order while removing duplicate
        return [x for x in factors if not (x in seen or seen.add(x))]
    else:
        return factors

cdef class Haversine:
    """
    compute the haversine function of an angle θ

    formula haversine function:
    haversine(θ) = sin^2(θ/2) = (1 - cos(θ)) / 2

    Attribute:
        theta (double): angle in radian passed during initialization

    Methods
        compute(self) (double): compute haversine(theta)
        get_theta(self) (double): get current theta value
        set_theta(self, double new_theta): set new angle in raadians

    Example:
    >>> h = Haversine(0)
    >>> h.compute()
    0.0
    >>> import math
    >>> h.set_theta(math.pi / 2)
    >>> h.compute()
    0.5
    """
    cdef double theta

    def __cinit__(self, double theta):
        """
        constructor for haversine class

        called once when object is created

        Parameter:
            theta (double): angle in radians - must be numeric float or double
        """
        self.theta = theta

    def __dealloc__(self):
        """
        currently empty because no dynamic memory allocation
        is used
        """
        pass

    cpdef double compute(self):
        """
        compute haversine function
        
        Return:
            (double): value of haversine(theta) = (1 - cos(theta)) / 2
        """
        # cosine based definition for avoiding compute square root or power
        return (1.0 - cos(self.theta)) / 2.0

    def get_theta(self):
        """
        get the current angle stored in the haversine instance
        usefull for debug or chaining operations
        
        Return:
            (double): current theta (in radians)
        """
        return self.theta

    def set_theta(self, double new_theta):
        """
        update internal angle theta

        allow reuse of same object for multiple computation

        Parameter:
            new_theta (double): new angle in radians
        """
        self.theta = new_theta

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

cpdef long jordan_totient(int n, int k) except -1:
    """
    compute jordan totient function which generalizes euler totient function

    Parameter:
        n (int): positive integer 
        k (int): non-negative integer exponent

    Return:
        (long): value of totient function

    Example:
    >>> jordan_totient(6, 1)
    2
    >>> jordan_totient(10, 2)
    80
    """
    if n <= 0:
        raise ValueError("input `n` must positive integers")
    if k < 0:
        raise ValueError("exponent `k` must non-negative integer")

    if k == 0:
        return 0

    if n == 1:
        return 1

    # acompute n^k
    cdef long res = <long>(pow(n, k))
    # get unique prime factor of n
    cdef list primes = prime_factors(n)
    # apply the formula for each prime factor
    cdef int p
    cdef long pk, numerator, denominator

    for p in primes:
        # compute p^k
        pk = <long>(pow(p, k))
        # compute (p^k - 1) / p^k as two separate integer operations
        numerator = pk - 1
        denominator = pk
        # multiply result by numerator / denominator using integer
        # arithmetic to make sure division happens before multiplication
        # to preventing overflow
        res = res // denominator * numerator
    return res
