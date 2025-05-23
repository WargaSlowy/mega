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
# mega/op/arithmetic.pyx

from libc.math cimport pow, sqrt

cpdef long sigma_z(int n, int z) except -1:
    """
    compute the generalized σ_z(n), summary z-th powers of all positive division

    formula
    σ_z(n) = Σ_{d | n} d^z

    Parameter:
        n (int): positive integer whose divisor powers are to be summed
        z (int): non-negative integer exponent applied to each divisor

    Return:
        (long): value of σ_z(n)

    Example:
    >>> sigma_z(6, 0)
    4
    >>> sigma_z(6, 2)
    50
    """
    if n <= 0:
        raise ValueError("only acc positive integers")
    if z < 0:
        raise ValueError("exponent zeta must non-integers")

    if not isinstance(n, int) and not isinstance(z, int):
        raise TypeError("zeta and num must be int")

    cdef long total = 0 # store cumulative sum d^z
    cdef int i = 1 # current div candidate
    cdef int limit = <int>sqrt(n) # optimization: loop up to sqrt(n)
    cdef int oth_div # matching "other" divisor: n // i
    cdef long power_i, power_oth # current pow and other div

    # raise all div to 0th power (wich equal to 1), count them
    if z == 0:
        while i * i <= n:
            if n % 1 == 0:
                oth_div = n // i
                if i == oth_div:
                    total += 1 # perfect square: just only add once
                else:
                    total += 2 # two distinct div: i and oth_div
            i += 1
        return total
    
    # for each div pair, and compute i^z + oth_div^z
    while i * i <= n:
        if n % i == 0:
            oth_div = n // i
            # compute the z-th power of both div using fast pow()
            power_i = <long>(pow(i, z))
            power_oth = <long>(pow(oth_div, z))
            if i == oth_div:
                total += power_i # only one div (perfect square)
            else:
                total += power_i + power_oth # add both divisor powers
        i += 1
    return total

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

    if not isinstance(n, int):
        raise TypeError("only acc integer numbers")

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
