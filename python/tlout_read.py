"""@package tlout_read
Read a tlout output file from randoutbreaksim.

Load timeline information from a tlout output file, and call a callback function for each path.

Callback function: First argument contains the timelines for the current path, the second argument the index of the bin corresponding to t=0, and the third argument indicates if the path reached extinction or not.
The first argument can contain two or three timelines, depending on the output file. The first timeline is the number of currently infected individuals, the second timeline is the number of new infections, and the third timeline is the number of new positive tests.
"""

import numpy as np

def read(filename, callback):
    headertype=np.dtype('<u4, b')

    f = open(filename, 'rb')

    npers, hbf = np.fromfile(f, dtype=headertype, count=1)[0]

    hasreltime = not((hbf&7)==1)
    postestresults = ((hbf&8)==8)

    tltype=np.dtype('<u4')

    if(hasreltime):
        rechdrtype = np.dtype('<u4, <u4, u1')
        rechdr=np.fromfile(f, dtype=rechdrtype, count=1)

        if(postestresults):

            while(rechdr.size):
                nbins, t0index, extinction = rechdr[0]
                ret = np.array_split(np.fromfile(f, dtype=tltype, count=3*nbins),3)
                callback(ret, t0index, extinction)
                rechdr = np.fromfile(f, dtype=rechdrtype, count=1)

        else:

            while(rechdr.size):
                nbins, t0index, extinction = rechdr[0]
                ret = np.array_split(np.fromfile(f, dtype=tltype, count=2*nbins),2)
                callback(ret, t0index, extinction)
                rechdr = np.fromfile(f, dtype=rechdrtype, count=1)

    else:
        rechdrtype = np.dtype('<u4, u1')
        rechdr=np.fromfile(f, dtype=rechdrtype, count=1)

        if(postestresults):

            while(rechdr.size):
                nbins, extinction = rechdr[0]
                ret = np.array_split(np.fromfile(f, dtype=tltype, count=3*nbins),3)
                callback(ret, 0, extinction)
                rechdr = np.fromfile(f, dtype=rechdrtype, count=1)

        else:

            while(rechdr.size):
                nbins, extinction = rechdr[0]
                ret = np.array_split(np.fromfile(f, dtype=tltype, count=2*nbins),2)
                callback(ret, 0, extinction)
                rechdr = np.fromfile(f, dtype=rechdrtype, count=1)

    f.close()