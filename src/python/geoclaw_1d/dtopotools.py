r"""
GeoClaw dtopotools Module  `$CLAW/geoclaw_1d/src/python/geoclaw_1d/dtopotools.py`

First pass at 1d version of dtopotools

:Classes:

  - DTopography

"""

import os,sys
import numpy


# ==============================================================================
#  DTopography Base Class for 1-dimensional dtopo
# ==============================================================================

class DTopography1d(object):
    r"""Basic object representing moving topography in 1d

    """


    def __init__(self, path=None, dtopo_type=None):
        r"""DTopography initialization routine.

        See :class:`DTopography` for more info.

        """

        self.dZ = None
        self.times = []
        self.x = None
        self.delta = None
        self.path = path
        if path:
            self.read(path, dtopo_type)


    def read(self, path=None, dtopo_type=None, verbose=False):
        r"""
        Read in a dtopo file and use to set attributes of this object.

        :input:

         - *path* (path) - Path to existing dtopo file to read in.
         - *dtopo_type* (int) - Type of topography file to read.  Default is 3
            if not specified or apparent from file extension.
        """

        if path is not None:
            self.path = path
        else:
            if self.path is None:
                raise ValueError("Need to specify a path to a file.")
            else:
                path = self.path

        if dtopo_type is None:
            dtopo_type = topotools.determine_topo_type(path, default=3)


        if dtopo_type == 2 or dtopo_type == 3:
            fid = open(path)
            mx = int(fid.readline().split()[0])
            mt = int(fid.readline().split()[0])
            xlower = float(fid.readline().split()[0])
            t0 = float(fid.readline().split()[0])
            dx = float(fid.readline().split()[0])
            dt = float(fid.readline().split()[0])
            fid.close()

            xupper = xlower + (mx-1)*dx
            x=numpy.linspace(xlower,xupper,mx)
            times = numpy.linspace(t0, t0+(mt-1)*dt, mt)

            dZvals = numpy.array(numpy.loadtxt(path, skiprows=6), ndmin=2)
            if dtopo_type==2:
                dZ = reshape(dZvals,(mt,mx))
            elif dtopo_type==3:
                dZ = dZvals

            self.x = x
            self.times = times
            self.dZ = dZ

        else:
            raise ValueError("Only dtopo types 2, and 3 are supported,",
                             " given %s." % dtopo_type)


    def write(self, path=None, dtopo_type=None):
        r"""Write out subfault resulting dtopo to file at *path*.

        :input:

         - *path* (path) - Path to the output file to written to.
         - *dtopo_type* (int) - Type of topography file to write out.  Default
           is 3.

        """

        if path is not None:
            self.path = path
        if self.path is None:
            raise IOError("*** need to specify path to file for writing")
        path = self.path

        if dtopo_type is None:
            dtopo_type = topotools.determine_topo_type(path, default=3)

        x = self.x
        dx = x[1] - x[0]

        # Construct each interpolating function and evaluate at new grid
        ## Shouldn't need to interpolate in time.
        with open(path, 'w') as data_file:

            if dtopo_type == 2 or dtopo_type == 3:
                if len(self.times) == 1:
                    dt = 0.
                else:
                    dt = float(self.times[1] - self.times[0])
                # Write out header
                data_file.write("%7i       mx \n" % x.shape[0])
                data_file.write("%7i       mt \n" % len(self.times))
                data_file.write("%20.14e   xlower\n" % x[0])
                data_file.write("%20.14e   t0\n" % self.times[0])
                data_file.write("%20.14e   dx\n" % dx)
                data_file.write("%20.14e   dt\n" % dt)

                if dtopo_type == 2:
                    for (n, time) in enumerate(self.times):
                        for j in range(len(x)):
                            data_file.write('%012.6e\n' % self.dZ[n,j])

                elif dtopo_type == 3:
                    for (n, time) in enumerate(self.times):
                        data_file.write(self.x.shape[0] * '%012.6e  '
                                              % tuple(self.dZ[n,:]))
                        data_file.write("\n")

            else:
                raise ValueError("Only dtopo_type 2 and 3 are ",
                                 "supported in 1d, not dtopo_type=%s." \
                                 % dtopo_type)


    def dZ_at_t(self, t):
        """
        Interpolate dZ to specified time t and return deformation.
        """
        if t <= self.times[0]:
            return self.dZ[0,:]
        elif t >= self.times[-1]:
            return self.dZ[-1,:]
        else:
            n = max(numpy.where(self.times <= t)[0])
            t1 = self.times[n]
            t2 = self.times[n+1]
            dz = (t2-t)/(t2-t1) * self.dZ[n,:] + \
                 (t-t1)/(t2-t1) * self.dZ[n+1,:]
            return dz