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

from libc.math cimport pow, sqrt, log

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
    cdef object z

    def __cinit__(self, int n, object z):
        """
        constructor that validating input before storing value

        Raise:
            ValueError: if n <= 0 or z 0

        Example:
        >>> s = SigmaZ(6, 2)
        """
        if n <= 0:
            raise ValueError("only acc positive integer")
        if isinstance(z, (float, int)):
            self.z = <double>z
        elif isinstance(z, complex):
            self.z = z
        else:
            raise TypeError("exponent must be numeric")

        self.n = n
        self.z = z

    def __dealloc__(self):
        pass

    cpdef object compute(self):
        """
        compute σ_z(n), the sum of the z-th powers of all positive divisors of n

        Return:
            (long): computed value of σ_z(n)
        """
        cdef int i = 1
        cdef int oth_div
        cdef double power_i, power_oth
        cdef double z_real
        cdef complex z_complex
        cdef long total_real = 0
        
        # when z == 0, just counting number of divisor
        if self.z == 0:
            total_real = 0
            while i * i <= self.n:
                if self.n % i == 0:
                    oth_div = self.n // i
                    if i == oth_div:
                        total_real += 1 # perfect square
                    else:
                        total_real += 2 # two distinct divisor
                i += 1
            return total_real

        # try convert z to double for optimal compute
        try:
            z_real = <double>self.z
        except TypeError:
            z_real = -1.0 # mark as invalid for fallback logic

        # if z is numeric and >= 0
        if z_real >= 0:
            total_real = 0
            while i * i <= self.n:
                if self.n % i == 0:
                    oth_div = self.n // i
                    power_i_real = <long>(pow(i, z_real))
                    power_oth_real = <long>(pow(oth_div, z_real))

                    if i == oth_div:
                        total_real += power_i_real # add once for perfect square
                    else:
                        total_real += power_i_real + power_oth_real
                i += 1
            return total_real


        import numpy as np
        z_complex = self.z
        total_complex = 0j

        # same loop structure, using numpy for complex support
        while i * i <= self.n:
            if self.n % i == 0:
                oth_div = self.n // i
                total_complex += np.power(i, z_complex)
                if i != oth_div:
                    total_complex += np.power(oth_div, z_complex)
            i += 1

        return total_complex
       
    def __repr__(self):
        """
        Return string representation of the object
        """
        return f"SigmaZ({self.n}, {self.z})"

cdef class EulerPhi:
    """
    compute euler totient function which is count the number of integer <= n that
    are coprime to n

    Formula:
        φ(n) = n × ∏(p|n) (1 - 1/p)

    Attribute:
        n (int): positive integer at which to evaluate EulerPhi

    Methods:
        compute(): compute using optimized logic
        __repr__(): string representation of class
        __getitem__(key): retrieved cached or compute new totient value
        __setitem__(key, value): manually store precomptued totient value

    Example:
    >>> e = EulerPhi(10)
    >>> e.compute()
    4
    >>> e[6] = 2
    >>> e[6]
    2

    Info:
        - using floating point math to avoid rounding error
        - safe up to ~n = 10^8 on mostly system
    """
    cdef public int n
    cdef dict _cache

    def __cinit__(self, int n):
        """
        initialize an instance of EulerPhi with input validation

        Parameter:
            n (int): must be >= 1 - domain of phi(n)
        """
        if n < 1:
            raise ValueError("only positive integer are acc")
        self.n = n
        self._cache

    def __dealloc__(self):
        pass

    cpdef int compute(self):
        """
        compute euler phi using integer arithmetic

        Return:
            (int): value of phi(n)
        """
        if self.n == 1:
            return 1

        cdef int result = self.n
        cdef int i = 2

        # handle factor 2 separately
        if self.n % i == 0:
            result -= result // i
            while self.n % i == 0:
                self.n //= i

        i += 1

        # check odd factors up to sqrt(n)
        while i * i <= self.n:
            if self.n % i == 0:
                result -= result // i
                while self.n % i == 0:
                    self.n //= i
            i += 1

        # if remaining n is a prime > 2, apply one final adjustment
        if self.n > 1:
            result -= result // self.n
        return result

    def __getitem__(self, int key):
        """
        retrueve phi(key) from internal cache, or compute and stored it

        allow dictionary-style access:
            e[6] -> return phi(6) = 2

        Parameter:
            key (int): positive integer (> 0)

        Return:
            (int): phi(key)
        """
        if key < 1:
            raise ValueError("only positive integer are acc")

        if key in self._cache:
            return self._cache[key]

        # create temporary instance to compute
        temp = EulerPhi(key)
        result = temp.compute()
        self._cache[key] = result
        return result

    def __setitem__(self, int key, int value):
        """
        manually set cached value for phi(key)

        usefully for precomputing or testing

        Parameter:
            key (int): positive integer (> 0)
            value (int): precomptued value of phi(key)
        """
        if key < 1:
            raise ValueError("key must be positive integer")
        self._cache[key] = value

    def __repr__(self):
        return f"EulerPhi({self.n})"

cdef class Chebyshev:
    """
    compute chebyshev function ϑ(x), with formula:

    ϑ(x) = Σ_{p ≤ x} log(p)

    Attribute:
        x (double): upper bound of summation

    Methods:
        compute(): compute ϑ(x) using trial division
        __repr__(): representation of chebyshev function
        __getitem__(key): retrieve cache value or compute it
        __setitem__(key, value): manually cache a compute result

    Example:
    >>> cheb = Chebyshev(10.0)
    >>> cheb.compute()
    5.347107530637821
    """
    cdef public double x
    cdef dict _cache

    def __cinit__(self, double x):
        """
        initialize chebyshev instance with input validation

        Parameter:
            x (double): must be >= 2
        """
        if x < 2:
            raise ValueError("x must be >= 2")
        self.x = x
        self._cache = {}

    def __dealloc__(self):
        self._cache.clear()

    cpdef double compute(self):
        """
        compute the first chebyshev function

        this method using trial divison to check primality
        then accumulate log(p) for each prime <= floor(x)

        Return:
            (double): result of chebyshev functon
        """
        cdef double result = 0.0
        cdef int i = 2
        cdef int limit = <int>self.x
        cdef int j, is_prime

        while i <= limit:
            is_prime = 1
            for j in range(2, <int>sqrt(i) + 1):
                if i % j == 0:
                    is_prime = 0
                    break
            if is_prime:
                result += log(i)
            i += 1
        return result
 
    def __getitem__(self, double key):
        """
        retrieve chebyshev function from internal cache, or compute
        and store it

        Parameter:
            key (double): positive number >= 2

        Return:
            (double): Chebyshev(key)
        """
        if key in self._cache:
            return self._cache[key]
        temp = Chebyshev(key)
        result = temp.compute()
        self._cache[key] = result
        return result

    def __setitem__(self, double key, double value):
        """
        manually set cached value for chebyshev(key)

        Parameter:
            key(double): positive number >= 2
            value(doube): precomputed value of chebyshev(key)
        """
        self._cache[key] = value

    
    def __repr__(self):
        return f"Chebyshev(x={self.x})"
