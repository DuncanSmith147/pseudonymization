
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

import math

"""
A collection of functions for operating on bitstrings.
"""

def digits(n, pad=None):
    """
    Returns Python string representing I{n}.
    """
    vals = [str(val) for val in iterbits(n)]
    if pad is not None:
        vals.extend(['0']*(pad-len(vals)))
    vals.reverse()
    return ''.join(vals)

def getbit(n, i):
    """
    Returns 0 or 1, the bit-value of bit I{i} of I{n}.
    """
    return (n >> i) & 1

def hamdist(x, y):
    """
    Returns the Hamming distance (number of bit-positions
    where the bits differ) between I{x} and I{y}.
    """
    return popcount(x ^ y)

def lowbits(n, i):
    """
    Returns the I{i} lowest bits of I{n}.
    """
    if not i > 0:
        raise ValueError, 'number of bits must be > 0'
    res = 0
    for ind, bit in enumerate(iterbits(n)):
        if ind == i:
            break
        res += bit * 2**ind
    return res

def numdigits(n):
    """
    Returns length of string representing I{n}.
    """
    if n == 0:
        return 1
    else:
        return int(math.log(n, 2)) + 1

def popcount(n):
    """
    Returns the number of 1-bits set in I{n}.
    """
    cnt = 0
    while n:
        n &= n - 1
        cnt += 1
    return cnt

def scan0(n, i=0):
    """
    Returns the bit-index of the first 0-bit of I{n} (that
    is at least I{i}).
    """
    n = n >> i
    return int(math.log((n ^ n+1) + 1, 2) - 1) + i

def scan1(n, i=0):
    """
    Returns the bit-index of the first 1-bit of I{n} (that
    is at least I{i}).
    """
    n = n >> i
    if n:
        return int(math.log(n & -n, 2)) + i
    else:
        return None

def setbit(n, i, val=1):
    """
    Returns a copy of the value of I{n}, with bit I{i} set
    to value I{val}.
    """
    if val:
        return n | (1 << i)
    else:
        return n & ~ (1 << i)

def flipbit(n, i):
    """
    Returns a copy of the value of I{n}, with bit I{i} flipped.
    """
    return n ^ (1 << i)

def iterbits(n):
    """
    Returns a generator of the bits of I{n}, up to the most
    significant bit.
    """
    while n:
        yield n & 1
        n = n >> 1


if __name__ == '__main__':
    try:
        import gmpy
    except ImportError:
        pass
    else:
        ints = range(400)
        for n in ints:
            try:
                assert gmpy.digits(n, 2) == digits(n)
            except AssertionError:
                print 'digits fail %d' % n
                raise
            try:
                assert gmpy.numdigits(n, 2) == numdigits(n)
            except AssertionError:
                print 'numdigits fail %d' % n
                raise
            try:
                assert gmpy.popcount(n) == popcount(n)
            except AssertionError:
                print 'popcount fail %d' % n
                raise
        for n in list(ints):
            for i in ints:
                try:
                    assert gmpy.getbit(n, i) == getbit(n, i)
                except AssertionError:
                    print 'getbit fail %d, %d' % (n, i)
                    raise
                try:
                    assert gmpy.setbit(n, i) == setbit(n, i)
                except AssertionError:
                    print 'setbit fail %d', n
                    raise
                try:
                    assert gmpy.lowbits(n, i + 1) == lowbits(n, i + 1)
                except AssertionError:
                    print 'lowbits fail %d', n
                    raise
                try:
                    assert gmpy.hamdist(n, i) == hamdist(n, i)
                except AssertionError:
                    print 'hamdist fail %d', n
                    raise
                try:
                    assert abs(flipbit(n, i) - n) == 2**i
                except AssertionError:
                    print 'flipbit fail %d', n
                    raise
                try:
                    assert gmpy.scan0(n, i) == scan0(n, i)
                except AssertionError:
                    print 'scan0 fail %d', n
                    raise
                try:
                    assert gmpy.scan1(n, i + 1) == scan1(n, i + 1)
                except AssertionError:
                    print 'scan0 fail %d', n
                    raise
        assert scan1(0) is None




