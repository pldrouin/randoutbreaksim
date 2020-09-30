"""@package ctout_read
Read a ctout output file from randoutbreaksim.

Load timeline information from a ctout output file.

"""

import numpy as np

def read(filename):
    entrytype=np.dtype('<i4, <i4, <u4, <u4, <u4')

    f = open(filename, 'rb')

    entries = np.fromfile(f, dtype=entrytype)

    f.close()

    return entries
