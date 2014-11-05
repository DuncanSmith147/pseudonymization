
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

from hashlib import md5

import numpy as np


class SimpleTabulation(object):
    # simple tabulation hash for variable
    # length strings using md5 intermediate
    # representation
    # admissible hash lengths are 8, 16, 32 and 64
    sizes = [8, 16, 32, 64]
    int_types = {8:np.uint8, 16:np.uint16, 32:np.uint32, 64:np.uint64}
    def __init__(self, q=64, seed=None):
        # if seed is not None, then its value will
        # be used to seed the PRNG
        # q denotes bit-length of hash returned by the hash method
        if not q in self.sizes:
            raise ValueError ('Invalid parameter value for q')
        self.q = q
        self.int_type = self.int_types[q]
        if not seed is None:
            np.random.seed(seed)
        if self.q == 64:
            self.tables = (np.random.random_integers(0, 2**32-1, size=(8, 256)).astype(np.uint64) +
                           (np.random.random_integers(0, 2**32-1, size=(8, 256)).astype(np.uint64) << 32))
        else:
            self.tables = np.random.random_integers(0, 2**q-1, size=(self.q//8, 256)).astype(self.int_type)
        self._cache = {}

    def hash(self, s):
        if not s in self._cache:
            x = self.int_type(int(md5(s).hexdigest(), 16) >> (128-self.q))
            h = self.int_type(0)
            for i in range(self.tables.shape[0]):
                c = np.uint8(x)
                h ^= self.tables[i,c]
                #x = x >> 8
                x = self.int_type(x//256) # bit shifting exposes bug in numpy for 64 bit unsigned int
            self._cache[s] = int(h)
        return self._cache[s]
