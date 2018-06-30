
##Copyright (c) 2014 duncan g. smith
##
##Permission is hereby granted, free of charge, to any person obtaining a
##copy of this software and associated documentation files (the "Software"),
##to deal in the Software without restriction, including without limitation
##the rights to use, copy, modify, merge, publish, distribute, sublicense,
##and/or sell copies of the Software, and to permit persons to whom the
##Software is furnished to do so, subject to the following conditions:
##
##The above copyright notice and this permission notice shall be included
##in all copies or substantial portions of the Software.
##
##THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
##OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
##FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
##THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
##OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
##ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
##OTHER DEALINGS IN THE SOFTWARE.

from __future__ import division

from math import log

from bitstring import hamdist, digits, popcount
from tabhash import SimpleTabulation
from set_like import BloomFilter


class Minwise(object):
    def __init__(self, m, q=64, seed=None, klass=SimpleTabulation):
        # m is the length of the list of hashes returned by the hash method
        # if seed is not None, then its value will
        # be used to seed the PRNG
        # numbits denotes length of hash returned by the hash method
        if not m > 0:
            raise ValueError('Minwise hash must have length > 0')
        if not seed is None:
            np.random.seed(seed)
        self._hashers = [klass(q=q) for _ in range(m)]
        self._cache = {}

    @property
    def hashers(self):
        return self._hashers

    def hash(self, tokens):
        # returns a list of minwise hashes
        # ensure tokens is a frozenset
        if not isinstance(tokens, frozenset):
            tokens = frozenset(tokens)
        if not tokens in self._cache:
            self._cache[tokens] = [min(h.hash(s) for s in tokens) for h in self.hashers]
        return self._cache[tokens]


class B_bit(Minwise):
    def __init__(self, b, m, q=64, seed=None):
        super(B_bit, self).__init__(m, q, seed)
        self._mod = 2**b

    def hash(self, tokens):
        # returns a list of b-bit minwise hashes
        if not isinstance(tokens, frozenset):
            tokens = frozenset(tokens)
        if not tokens in self._cache:
            mod = self._mod
            self._cache[tokens] = [x % mod for x in super(B_bit, self).hash(tokens)]
        return self._cache[tokens]


class Concatenated(B_bit):
    def __init__(self, m, q=64, seed=None):
        super(Concatenated, self).__init__(1, m, q, seed)

    def hash(self, tokens):
        # returns a concatenated 1-bit hash
        if not isinstance(tokens, frozenset):
            tokens = frozenset(tokens)
        if not tokens in self._cache:
            h = 0
            for i, bit in enumerate(reversed(super(Concatenated, self).hash(tokens))):
                if bit:
                    h += (1 << i)
            self._cache[tokens] = C_hash(h, len(self.hashers))
        return self._cache[tokens]


class C_hash(long):
    # a 1-bit concatenated hash instance
    # i.e. a Python long with a couple of
    # of additional attributes
    def __new__(cls, bits, m, N=1):
        obj = long.__new__(cls, bits)
        obj._m = m # number of significant bits
        obj._N = N # compression factor for XOR compression
        return obj

    def hex(self):
        # return a hex representation
        return hex(self)

    def compressed(self, m):
        # returns a compressed hash of length m
        # by simply reducing number of significant bits
        if not m <= self._m:
            raise ValueError('Cannot compress to larger size')
        return self.__class__.__new__(self.__class__, self % 2**m, m)

    def XOR(self, N):
        # returns a compressed hash using XOR
        # N must be a positive power of 2
        # (and for most applications probably not more
        # than around 8)
        if not popcount(N) == 1:
            raise ValueError(' N must be a positive power of 2')
        chunksize = self.m // N
        if not N * chunksize == self.m:
            raise ValueError('Hash size %s not divisible by %s' % (self._m, N))
        h = 0
        x = self
        mod = 2**chunksize
        for _ in range(N):
            h ^= (x % mod)
            x = x >> chunksize
        return self.__class__.__new__(self.__class__, h, chunksize, N=self._N*N)

    @ property
    def digits(self):
        # returns the bitstring representation
        return digits(self, pad=self._m)

    @property
    def m(self):
        return self._m

    @property
    def N(self):
        return self._N


############## Functions based on token sets ##############


def Jaccard(tokens1, tokens2):
    # returns exact Jaccard
    # similarity measure for
    # two token sets
    tokens1 = set(tokens1)
    tokens2 = set(tokens2)
    return len(tokens1&tokens2) / len(tokens1|tokens2)

def dice_coefficient(tokens1, tokens2):
    # returns exact dice coefficient
    # similarity measure for
    # two token sets
    tokens1 = set(tokens1)
    tokens2 = set(tokens2)
    return 2 * len(tokens1&tokens2) / (len(tokens1) + len(tokens2))


############## Functions for b-bit hashes ##############


def J_hat(hashes1, hashes2, b):
    # returns estimated Jaccard score from lists of
    # minwise hashes 
    # hashes1 and hashes2 are lists for different
    # token sets using the same (minwise) hash functions
    # only the b most significant bits are used
    # in the evaluation of the measure
    # ensuring that b is 'much' less than the
    # number of bits in the full hash values
    # assumes the probability of collisions on the
    # full hash is negligible
    if not len(hashes1) == len(hashes2):
        raise ValueError('Hash lists must have equal length')
    mod = 2**b
    frac = sum([(a-b) % mod == 0 for (a,b) in zip(hashes1, hashes2)]) / len(hashes1)
    c = (1/2)**b
    return (frac-c)/(1-c)


############## Functions for 1-bit concatenated hashes ##############


def J_hat_from_conc(a, b, m, N=1, truncate=True):
    # returns estimated Jaccard
    # similarity measure for
    # a pair of comparable concatenated hashes
    # if truncate is True, then 0 is
    # returned rather than a negative value
    # N is the compression factor
    try:
        res = (1-2*hamdist(a,b)/m)**(1/N)
    except ValueError:
        # enforce truncation for N > 1
        return 0
    if truncate and res < 0:
        return 0
    return res

def var_J_hat(J, m, N=1):
    # returns variance of
    # above estimator (without truncation)
    # for true Jaccard score J, compressed hash length m
    # and XOR compression factor N
    # result is approximate for N > 1
    return (1/J**(2*N-2)-J**2)/m/N/N

def mean_J_hat(J, m, N=1):
    # returns expected value of
    # above estimator (without truncation)
    # for true Jaccard score J, compressed hash length m
    # and XOR compression factor N
    # result is approximate for N > 1
    return J - (N-1)*(1/J**(2*N-1)-J)/m/N/N/2

def MSE_J_hat(J, m, N=1):
    # returns MSE of
    # above estimator (without truncation)
    # for true Jaccard score J, compressed hash length m
    # and XOR compression factor N
    # result is approximate for N > 1
    if N == 1:
        return var_J_hat(J, m, N)
    return var_J_hat(J, m, N) + (J - mean_J_hat(J, m, N))**2


############## Bloom filter functions ##############


def J_hat_from_bf(a, b):
    # returns estimated Jaccard
    # similarity measure for
    # a pair of comparable Bloom filters
    return popcount(a&b) / popcount(a|b)

def D_hat_from_bf(a, b):
    # returns estimated Dice coefficient
    # similarity measure for
    # a pair of comparable Bloom filters
    return 2*popcount(a&b)/(popcount(a) + popcount(b))

def J_hat_from_bf_corrected(a, b, m):
    # bias corrected estimator due to Swamidass and Baldi (2007)
    A = -m*log(1-popcount(a)/m)
    B = -m*log(1-popcount(b)/m)
    AB = -m*log(1-popcount(a|b)/m)
    num = max(A+B-AB, 0)
    denom = min(AB, A+B)
    return num / denom
