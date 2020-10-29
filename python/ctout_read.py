"""@package ctout_read
Read a ctout output file from randoutbreaksim.

Load timeline information from a ctout output file.
Entry format is time of positive test time (in minutes), pre-symptomatic period duration(in minutes), child's ID, parent ID (negative if child's infectious period is not interrupted through CT), and the number of traced contacts.

"""

import numpy as np

def read(filename):
    entrytype=np.dtype('<i4, <i4, <i4, <i4, <u4')

    f = open(filename, 'rb')

    entries = np.fromfile(f, dtype=entrytype)

    f.close()

    return entries
