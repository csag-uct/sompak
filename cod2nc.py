import netCDF4
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser('Write out node centroid data from SOMPAK .cod file as NetCDF4 file')
parser.add_argument('codfile', type=str, help='cod filename')
parser.add_argument('--lat0', type=float, required=True, help='starting latitude')
parser.add_argument('--lon0', type=float, required=True, help='starting longitude')
parser.add_argument('--ysize', type=int, required=True, help='latitude grid points')
parser.add_argument('--xsize', type=int, required=True, help='longitude grid points')
parser.add_argument('--dlat', type=float, required=True, help='latitude spacing')
parser.add_argument('--dlon', type=float, required=True, help='longitude spacing')
parser.add_argument('--variables', nargs='*', type=str, default=None, help='list of variable names')
parser.add_argument('--scale', nargs='*', type=float, default=None, help='scaling values for each variable: output = variable * scale + offset')
parser.add_argument('--offset', nargs='*', type=float, default=None, help='offset values for each variable: output = variable * scale + offset')
parser.add_argument('--timevar', action='store_true', help='write nodes as time steps (useful for GrADS visualisation)')
parser.add_argument('-o', '--output', type=str, required=True, help='output NetCDF filename')
args = parser.parse_args()

# Open the COD file
infile = open(args.codfile)

# Get the sompack parameters
size, shape, somx, somy, method = infile.readline().split()
size = int(size)
print size, shape, somx, somy, method

# Construct the grid shape and size and calculate number of variables
gridshape = (args.ysize, args.xsize)
gridsize = gridshape[0] * gridshape[1]
varcount = size/gridsize
print gridshape, gridsize, varcount

# Process scale and offset arguments or set defaults
scale = np.ones((varcount,))
if args.scale:
	scale[:] = np.array([float(s) for s in args.scale])[:varcount]

offset = np.zeros((varcount,))
if args.offset:
	offset[:] = np.array([float(s) for s in args.offset])[:varcount]

# If we don't have variable names then construct them as varN
if args.variables:
	varnames = args.variables
else:
	varnames = ['var{}'.format(n) for n in range(varcount)]

# Calculate number of nodes and data shapes
nodes = int(somx) * int(somy)
datashape = (varcount, args.ysize, args.xsize)
varshape = (nodes, args.ysize, args.xsize)

# Setup variables dict
variables = {}
vi = 0
for varname in varnames:
	variables[varname] = np.zeros(varshape, dtype=np.float32)
	print varname, scale[vi], offset[vi]

# Now process the lines of the file
ni = 1
for line in infile.readlines():
	print("#node={}".format(ni))

	# Read the full line into a temporary array and reshape to datashape
	tmp = np.array([float(v) for v in line.split()], dtype=np.float32).reshape(datashape)

	# Now extract each variable section
	vi = 0
	for varname in varnames:
		variables[varname][ni-1,:] = (tmp[vi,:] * scale[vi]) + offset[vi]
		v = variables[varname][ni-1,:]
		print("\t{} (min, mean, max) = {:2f} {:2f} {:2f}").format(varname, v.min(), v.mean(), v.max())
		vi += 1

	ni += 1


# Create the output NetCDF4 file
ds = netCDF4.Dataset(args.output, 'w', format='NETCDF4_CLASSIC')

# Create time or node dimension depending on argument
if args.timevar:
	ds.createDimension('time', nodes)
else:
	ds.createDimension('node', nodes)

# Other dimensions	
ds.createDimension('lat', args.ysize)
ds.createDimension('lon', args.xsize)

# Latitude and longitude variables
lat = ds.createVariable('lat', np.float32, ['lat'])
lat.units = 'degrees_north'
lat[:] = np.arange(args.lat0, args.lat0 + (args.ysize * args.dlat), args.dlat)

lon = ds.createVariable('lon', np.float32, ['lon'])
lon.units = 'degrees_east'
lon[:] = np.arange(args.lon0, args.lon0 + (args.xsize * args.dlon), args.dlon)

# Create time values if needed
if args.timevar:
	time = ds.createVariable('time', np.int32, ['time'])
	time.units = 'days since 1900-01-01'
	time[:] = np.arange(nodes, dtype=np.int32)

# Create each output variable
for varname in variables:
	if args.timevar:
		var = ds.createVariable(varname, np.float32, ['time', 'lat', 'lon'])
	else:
		var = ds.createVariable(varname, np.float32, ['node', 'lat', 'lon'])

	var.short_name = varname
	var[:] = variables[varname][:]

# Finished!
ds.close()
	
