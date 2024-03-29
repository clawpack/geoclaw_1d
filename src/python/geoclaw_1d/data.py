#!/usr/bin/env python

"""

Classes representing parameters for 1D GeoClaw runs

:Classes:

 - GeoClawData1D -- has been removed, now using usual GeoClawData
 - GaugeData1D
 - GridData1D

:Constants:

 - Rearth - Radius of earth in meters
 - DEG2RAD factor to convert degrees to radians
 - RAD2DEG factor to convert radians to degrees

"""

import os
import numpy
import clawpack.clawutil.data

# Radius of earth in meters.
# For consistency, should always use this value when needed, e.g.
# in setrun.py or topotools:
Rearth = 6367.5e3  # average of polar and equatorial radii

DEG2RAD = numpy.pi / 180.0
RAD2DEG = 180.0 / numpy.pi
LAT2METER = Rearth * DEG2RAD


#  Gauge data object

#  Gauge data object
class GaugeData1D(clawpack.clawutil.data.ClawData):
    r"""
     Gauge data object for 1d.
     input specs for gauges are in 1d in setrun.py...output is like that of
      2d amr (with level=1 and y=0) so that same reading/plotting tools can be used.
    """

    @property
    def gauge_numbers(self):
        if len(self.gauges) == 1:
            return [self.gauges[0][0]]
        else:
            return [gauge[0] for gauge in self.gauges]

    def __init__(self, num_dim=2):
        super(GaugeData1D,self).__init__()

        self.add_attribute('num_dim',num_dim)
        self.add_attribute('gauges',[])

    def __str__(self):
        output = "Gauges: %s\n" % len(self.gauges)
        for gauge in self.gauges:
            output = "\t".join((output,"%4i:" % gauge[0]))
            output = " ".join((output,"%19.10e" % gauge[1]))
            output = " ".join((output,"%17.10e" % gauge[2]))
            output = " ".join((output,"%13.6e\n" % gauge[3]))
        return output

    def write(self,out_file='gauges.data',data_source='setrun.py'):
        r"""Write out gague information data file."""


        # Check to make sure we have only unique gauge numebrs
        if len(self.gauges) > 0:
            if len(self.gauge_numbers) != len(set(self.gauge_numbers)):
                raise Exception("Non unique gauge numbers specified.")

        # Write out gauge data file
        self.open_data_file(out_file,data_source)
        self.data_write(name='ngauges',value=len(self.gauges))
        for gauge in self.gauges:
            self._out_file.write("%4i %19.10e %19.10e %13.6e  %13.6e\n" % tuple(gauge))
        self.close_data_file()

    def read(self,data_path="./",file_name='gauges.data'):
        r"""Read gauge data file"""
        path = os.path.join(data_path, file_name)
        gauge_file = open(path,'r')

        # Read past comments and blank lines
        header_lines = 0
        ignore_lines = True
        while ignore_lines:
            line = gauge_file.readline()
            if line[0] == "#" or len(line.strip()) == 0:
                header_lines += 1
            else:
                break

        # Read number of gauges, should be line that was last read in
        num_gauges = int(line.split()[0])

        # Read in each gauge line
        for n in xrange(num_gauges):
            line = gauge_file.readline().split()
            self.gauges.append([int(line[0]),float(line[1]),float(line[2]),
                                                              float(line[3]),float(line[4])])

        gauge_file.close()


class GridData1D(clawpack.clawutil.data.ClawData):
    r"""
    1D data object for grid info

    """
    def __init__(self):
        super(GridData1D,self).__init__()

        self.add_attribute('grid_type',0)
        self.add_attribute('fname_celledges',None)
        self.add_attribute('monitor_fgmax',False)
        self.add_attribute('monitor_runup',False)
        self.add_attribute('monitor_total_zeta',False)

    def write(self,out_file='grid.data',data_source='setrun.py'):

        self.open_data_file(out_file,data_source)

        self.data_write('grid_type')
        if self.grid_type == 2:
            if self.fname_celledges is None:
                self.fname_celledges = 'celledges.txt'
                print('*** grid_type ==2 and fname_celledges not specified,')
                print('*** using celledges.txt')
            # if path is relative in setrun, assume it's relative to the
            # same directory that out_file comes from
            fname = os.path.abspath(os.path.join(os.path.dirname(out_file),
                                    self.fname_celledges))
            self._out_file.write("\n'%s'   =: fname_celledges\n " % fname)

        self._out_file.write("\n%s   =: monitor_fgmax" \
                             % str(self.monitor_fgmax)[0])
        self._out_file.write("\n%s   =: monitor_runup" \
                             % str(self.monitor_runup)[0])
        self._out_file.write("\n%s   =: monitor_total_zeta" \
                             % str(self.monitor_total_zeta)[0])
        self.close_data_file()

    def read(self, path, force=False):
        with open(os.path.abspath(path), 'r') as data_file:
            for line in data_file:
                if "=:" in line:
                    value, tail = line.split("=:")
                    varname = tail.split()[0]
                    if varname == 'grid_type':
                        self.grid_type = int(value)
                    elif varname == 'fname_celledges':
                        self.fname_celledges = value.strip()


class BoussData1D(clawpack.clawutil.data.ClawData):
    r"""
    1D data object for Boussinesq info

    """
    def __init__(self):
        super(BoussData1D,self).__init__()

        self.add_attribute('boussEquations',2)
        self.add_attribute('boussMinDepth',20.)

    def write(self,out_file='bouss.data',data_source='setrun.py'):

        self.open_data_file(out_file,data_source)

        self.data_write('boussEquations')
        self.data_write('boussMinDepth')

        self.close_data_file()

