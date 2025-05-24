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

cdef class SigmaZ:
    """    
    Compute the generalized σ_z(n), summary z-th powers of all positive division

    Formula:
    σ_z(n) = Σ_{d | n} d^z

    Special case:
        - wehn z == 0 -> sigma_z(n) return number of divisor (τ(n))
        - when z == 1 -> sigma_z(n) = summary all divisor d of n sums d^z

    Attribute:
        n (int): input number
        z (int): exponent used in computation

    Example:
    >>> sig = SigmaZ(6, 0)
    >>> sig.compute()
    4
    """
    cdef int n
    cdef int z

    def __cinit__(self, int n, int z):
        """
        constructor that validating input before storing value

        Raise:
            ValueError: if n <= 0 or z 0

        Example:
        >>> s = SigmaZ(6, 2)
        """
        if n <= 0:
            raise ValueError("only acc positive integer")
        if z < 0:
            raise ValueError("exponent must be non-negative")

        self.n = n
        self.z = z

    def __dealloc__(self):
        pass

    cpdef long compute(self):
        """
        compute σ_z(n), the sum of the z-th powers of all positive divisors of n

        Return:
            (long): computed value of σ_z(n)
        """
        cdef long total = 0
        cdef int i = 1
        cdef int limit = <int>sqrt(self.n)
        cdef int oth_div
        cdef long power_i, power_oth

        # σ₀(n) count the number of positive divisor of n
        if self.z == 0:
            while i * i <= self.n:
                if self.n % i == 0:
                    oth_div = self.n // i
                    if i == oth_div:
                        total += 1 # perfect square: count once
                    else:
                        total += 2 # two distrint divisor: i and oth_div
                i += 1
            return total

        # σ_z(n) = Σ d^z over all d dividing n
        while i * i <= self.n:
            if self.n % i == 0:
                oth_div = self.n // i
                
                # compute power using pow
                power_i = <long>(pow(i, self.z))
                power_oth = <long>(pow(oth_div, self.z))

                if i == oth_div:
                    total += power_i
                else:
                    total += power_i + power_oth
            i += 1
        return total

    def __repr__(self):
        """
        Return string representation of the object
        """
        return f"SigmaZ({self.n}, {self.z}) = {SigmaZ(self.n, self.z).compute()}"

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
