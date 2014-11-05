
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
import functools

from bitstring import digits, popcount
from tabhash import SimpleTabulation


def k_hashes(k, m, seed=None):
    # returns a list of k hash functions
    # suitable for a Bloom filter of length m
    # the PRNG used to generate hash functions can be seeded
    cache = {}
    # generate required size for hashes
    for q in SimpleTabulation.sizes:
        if 2**q >= m:
            break
    else:
        raise ValueError('m is too large (> 2**%d)' % q)
    hash1 = SimpleTabulation(q, seed)
    hash2 = SimpleTabulation(q)
    def f(item, i):
        key = (item, i)
        if not key in cache:
            cache[key] = hash1.hash(item) + i * hash2.hash(item)
        return cache[key]
    return [functools.partial(f, i=i) for i in range(k)]


class BloomFilter(object):
    def __init__(self, m, funcs, items=None):
        # Bloom filter of length m
        # using hash functions in funcs
        self._m = m # number of bits
        self.funcs = funcs
        self.bits = 0
        if items is not None:
            for item in items:
                self.add(item)

    def __str__(self):
        return str(self.bits)

    def __repr__(self):
        return repr(self.bits)

    def hex(self):
        # return a hex representation of self.bits
        return hex(self.bits)

    def add(self, item):
        for func in self.funcs:
            index = func(item) % self._m
            self.bits = self.bits | (1 << index)

    def __contains__(self, item):
        # can generate false positives
        for func in self.funcs:
            index = func(item) % self._m
            if (self.bits >> index) & 1 == 0:
                return False
        return True

    def union(self, other):
        # returns a Bloom filter containing
        # the union of items in self and other
        # (which can be any container type)
        # Assumes Bloom filter union is valid
        # (same hash functions etc.)
        res = self.__class__(self._m, self.funcs)
        if isinstance(other, self.__class__):
            res.bits = self.bits | other.bits
        else:
            res.bits = self.bits
            for item in other:
                res.add(item)
        return res

    def intersection(self, other):
        # returns a Bloom filter containing
        # the intersection of bits in self and other,
        # or the intersection of bits in self and
        # a Bloom filter with items added from other
        res = self.__class__(self._m, self.funcs)
        if isinstance(other, self.__class__):
            res.bits = self.bits & other.bits
        else:
            for item in other:
                res.add(item)
            res.bits = self.bits & res.bits
        return res

    def estimated_size(self):
        # estimate of number of elements in the Bloom filter
        # note: estimate of size of intersection
        # can be estimated as the sum of
        # the estimated sizes of A and B
        # minus the estimated size of A.union(B)
        return -self._m * log(1-popcount(self.bits)/self._m) / len(self.funcs)

    @ property
    def digits(self):
        # returns the bitstring representation
        return digits(self.bits, pad=self._m)

    @property
    def m(self):
        return self._m


def get_k(m, n, p=0.5):
    # returns an estimate of the number of hash functions k s.t.
    # the expected proportion of set bits is p, given m bits in
    # Bloom filter and n added items
    return log(1-p) / log(1-1/m) / n

def opt_k(m, n):
    # returns an estimate of the number of hash functions k that
    # minimizes the false positive rate for given m and n
    return (m/n) * log(2)


"""
>>> import set_like
>>> bf = set_like.BloomFilter(100, set_like.k_hashes(3, 10, seed=6))
>>> bf.add('one')
>>> bf.add('two')
>>> bf.add('three')
>>> 'two' in bf
True
>>> 'four' in bf
False
"""
