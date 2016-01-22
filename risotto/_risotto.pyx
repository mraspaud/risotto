# cython: profile=True
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2015 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Rice (USES) algorithm (also known as SZIP) decoder

TODO:
* implement restricted mode
* implement low entropy mode

"""

import numpy as np
cimport numpy as np
cimport cython

ctypedef np.uint8_t DTYPE_t


cdef postprocess(dec, size_t n, int ref_sample=-1):
    cdef int current_sample
    cdef np.ndarray[np.int_t, ndim = 1] samples = np.empty(len(dec), dtype=np.int)
    cdef size_t idx = 0
    if ref_sample == -1:
        current_sample = dec[0]
        dec = dec[1:]
        samples[idx] = current_sample
        idx += 1
    else:
        current_sample = ref_sample

    cdef int delta_i
    cdef int theta_i
    cdef int ti
    for i, delta_i in enumerate(dec):
        theta_i = min(2 ** n - 1 - current_sample, current_sample - 0)

        if delta_i <= (2 * theta_i):
            if delta_i % 2:
                ti = -(delta_i + 1) / 2
            else:
                ti = delta_i / 2
        else:
            if theta_i == (current_sample - 0):
                ti = delta_i - theta_i
            else:
                ti = theta_i - delta_i
        current_sample += ti
        samples[idx] = current_sample
        idx += 1
    return samples


cdef read_block(size_t nbits, size_t block_size):
    cdef np.ndarray[np.int_t, ndim = 1] res = np.empty(block_size, dtype=int)
    for i in range(block_size):
        res[i] = read(nbits)
    return res


cdef unsigned char[:] masks


def init_arrays():
    global masks
    masks = np.zeros(9, dtype=np.uint8)
    masks[0] = 0
    for i in range(1, len(masks)):
        masks[i] = <unsigned char > ((1 << i) - 1)

init_arrays()


cdef unsigned char[:] data
cdef size_t offset = 0


@cython.boundscheck(False)
cdef inline read(size_t nbits):
    global data  # unsigned char[:]
    global offset  # size_t
    cdef unsigned int res = 0
    cdef size_t bits_left
    while nbits > 0:
        if len(data) == 0:
            raise ValueError("no more data to read")
        bits_left = 8 - offset
        if bits_left >= nbits:
            if res != 0:
                res <<= nbits
            res += ((1 << nbits) - 1) & (data[0] >> (bits_left - nbits))
            offset = (8 - (bits_left - nbits)) % 8
            nbits = 0
        else:
            if res != 0:
                res <<= bits_left
            res += ((1 << bits_left) - 1) & data[0]
            nbits -= bits_left
            offset = 0
        if offset == 0:
            data = data[1:]
    return res


@cython.boundscheck(False)
cdef inline decode_split(size_t n, size_t k, size_t block_size, bint is_ref):
    cdef int m
    cdef int sample
    cdef np.ndarray[np.int_t, ndim = 1] dec
    cdef int bit
    cdef int count
    cdef size_t start_offset = 0
    cdef np.ndarray[np.int_t, ndim = 1] fs
    m = 1 << k
    fs = np.zeros(block_size, dtype=np.int)
    # dec = []
    # ref_sample = None
    dec = np.empty(block_size, dtype=np.int)

    if is_ref:
        start_offset = 1
        dec[0] = read(n)

    for i in range(start_offset, block_size):
        bit = read(1)
        while bit == 0:
            fs[i] += 1
            bit = read(1)
    for i in range(start_offset, block_size):
        dec[i] = fs[i] * m + read(k)

    return dec


cpdef decode(unsigned char[:] data_array, size_t n, size_t j, size_t reference=0):

    global data
    data = data_array
    chunks = []
    cdef size_t option_len
    cdef size_t bcount
    cdef size_t block_size
    cdef size_t tot_len = 0
    cdef np.ndarray[np.int_t, ndim = 1] dec
    cdef np.ndarray[np.int_t, ndim = 1] decsamples
    cdef size_t kidx = 0
    cdef int gamma
    cdef int beta
    cdef int ms
    if n <= 8:
        option_len = 3
    elif n <= 16:
        option_len = 4
    else:
        option_len = 5

    sample = -1

    bcount = 0
    try:
        while len(data):
            option = read(option_len)
            block_size = j
            if option == masks[option_len]:
                # no compression
                # print "no compression"
                #if reference and bcount == 0:
                #    block_size -= 1
                #    # ref_sample = read(n)
                #    sample = read(n)
                # dec = [read(n) for i in range(block_size)]
                dec = read_block(n, block_size)
                decsamples = postprocess(dec, n, sample)
                sample = decsamples[-1]
                bcount += len(decsamples)
                if reference:
                    bcount %= reference
                    # print "bcount", bcount

            elif option > 0:
                # # split-k compression
                # print "split-k", option - 1
                if reference and bcount == 0:
                    sample = -1
                dec = decode_split(n, option - 1, block_size, reference and bcount == 0)
                decsamples = postprocess(dec, n, sample)
                # if ref_sample:
                #    decsamples = [ref_sample] + decsamples
                bcount += len(decsamples)
                if reference:
                    bcount %= reference
                    # print "bcount", bcount
                sample = decsamples[-1]
            elif option == 0:
                bit = read(1)
                if bit == 1:
                    # second extension
                    # print "second"
                    # ref_sample = None
                    if reference and bcount == 0:
                        sample = read(n)
                    # dec = []
                    dec = np.empty(j, dtype=np.int)
                    kidx = 0
                    for i in range(j / 2):
                        bit = read(1)
                        gamma = 0
                        while bit == 0:
                            gamma += 1
                            bit = read(1)
                        if gamma == 0:
                            beta = 0
                            ms = 0
                        elif gamma < 3:
                            beta = 1
                            ms = 1
                        elif gamma < 6:
                            beta = 2
                            ms = 3
                        elif gamma < 10:
                            beta = 3
                            ms = 6
                        elif gamma < 15:
                            beta = 4
                            ms = 10
                        elif gamma < 21:
                            beta = 5
                            ms = 15
                        elif gamma < 28:
                            beta = 6
                            ms = 21
                        elif gamma < 36:
                            beta = 7
                            ms = 28
                        elif gamma < 45:
                            beta = 8
                            ms = 36
                        else:
                            raise ValueError("wtf?")

                        delta_i1 = gamma - ms
                        delta_i = beta - delta_i1
                        # dec.append(delta_i)
                        # dec.append(delta_i1)

                        dec[kidx] = delta_i
                        kidx += 1
                        dec[kidx] = delta_i1
                        kidx += 1
                    decsamples = postprocess(dec, n, sample)
                    # if ref_sample:
                    #    decsamples = [ref_sample] + decsamples
                    bcount += len(decsamples)
                    if reference:
                        bcount %= reference
                        # print "bcount", bcount
                    sample = decsamples[-1]
                else:
                    # zero blocks
                    # print "zeros"
                    # ref_sample = None
                    if reference and bcount == 0:
                        sample = read(n)
                        block_size -= 1
                    bit = read(1)
                    count = 0
                    while bit == 0:
                        count += 1
                        bit = read(1)
                    if count < 4:
                        zeros = (count + 1) * j
                    if count == 4:
                        zeros = reference - bcount
                    if count > 4:
                        zeros = count * j
                    dec = np.zeros(zeros, dtype=np.int)
                    decsamples = postprocess(dec, n, sample)
                    # if ref_sample:
                    #    decsamples = [ref_sample] + decsamples
                    bcount += len(decsamples)
                    if reference:
                        bcount %= reference
                        # print "bcount", bcount
                    sample = decsamples[-1]
            chunks.append(decsamples)
            tot_len += len(decsamples)
    except (IndexError, ValueError) as err:
        # print "finished ?"
        if len(chunks) != j:
            # print len(chunks), tot_len
            # raise
            pass
    return chunks


def main():

    for N in range(5, 32):
        nstr = "{0:02d}".format(N)

        if N <= 8:
            dtype = np.uint8
            basename = "/home/a001673/usr/src/rice/121B2TestData/AllOptions/test_p256n"
            ref = 256
        elif N <= 16:
            dtype = np.uint16
            basename = "/home/a001673/usr/src/rice/121B2TestData/AllOptions/test_p256n"
            ref = 256
        elif N <= 32:
            dtype = np.uint32
            basename = "/home/a001673/usr/src/rice/121B2TestData/AllOptions/test_p512n"
            ref = 512
        else:
            raise ValueError("Something wrong happened")
        ref_data = np.fromfile(basename + nstr + ".dat", dtype=dtype)

        #comp = np.fromfile(basename + nstr + ".rz", dtype=np.uint8)
        with open(basename + nstr + ".rz") as fd:
            comp = bytearray(fd.read())
            # print "filesize:", len(comp)

        J = 16
        print basename + nstr + ".dat"
        for bcount, chunk in enumerate(decode(comp, N, J, ref)):
            refsamples = ref_data[bcount * J:(bcount + 1) * J]
            assert(all(chunk[:len(refsamples)] == refsamples[:len(chunk)]))
            # print "passed", bcount, "with", len(chunk), "samples"
