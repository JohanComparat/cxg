import time
t0 = time.time()
import numpy as np
import numpy as n
import healpy
from astropy.table import Table, Column, vstack, hstack
import sys, os, glob

import astropy.units as u
import astropy.constants as cc
import astropy.io.fits as fits

from astropy.wcs import WCS
from sklearn.neighbors import BallTree

topdir    = sys.argv[1]
p2_data   = sys.argv[2]
C_topdir    = sys.argv[3]
C_p2_data   = sys.argv[4]
C_p2_random = sys.argv[5]
p_2_OUT  = os.path.join(topdir, sys.argv[6])
#if os.path.isfile(p_2_OUT)==False:
basename = sys.argv[4][:-4]
mag_1 = sys.argv[7]
mag_2 = sys.argv[8]
z_min = float(sys.argv[9])
z_max = float(sys.argv[10])

D2 = Table.read(os.path.join(C_topdir, C_p2_data ), format='fits')
R2 = Table.read(os.path.join(C_topdir, C_p2_random ), format='fits')
D2=D2[(D2['Z_LAMBDA']>=z_min)&(D2['Z_LAMBDA']<=z_max)]
R2=R2[(R2['redshift']>=z_min)&(R2['redshift']<=z_max)]

#import astropy.units as u
#import astropy.constants as cc
#from astropy.cosmology import FlatLambdaCDM
#cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
#h = 0.6774
#L_box = 1000.0 / h
#cosmo = cosmoUNIT
z_bar = np.mean(D2['Z_LAMBDA'])

print('reading', os.path.join(topdir, p2_data), time.time()-t0, 's')
GAL_i = Table.read(os.path.join(topdir, p2_data))
color = GAL_i[mag_1]-GAL_i[mag_2]
color_shift = 0.05
color_step = 0.01
color_grid = np.arange(np.round(color.min(),2), np.round(color.max(),2)-color_step, color_step)

deg_to_rad = np.pi/180.
tree_array = {}
for jj, c_min in enumerate(color_grid):
	selection = (color>=c_min)&(color<c_min+color_shift)
	coord_GAL = deg_to_rad * np.transpose([GAL_i['DEC'][selection], GAL_i['RA'][selection] ])
	Tree_GAL = BallTree(coord_GAL, metric='haversine')
	tree_array[jj] = Tree_GAL
print('trees constructed', time.time()-t0, 's')

coord_CLU = deg_to_rad * np.transpose([ D2['DEC'], D2['RA'] ])

amin = 10.
radius_max = amin / 60. * n.pi / 180. # 10 arcminutes
asec = 10.
radius_step = asec / 3600. * n.pi / 180. # 10 arcsecond
radii = np.arange(radius_step, radius_max+radius_step, radius_step)

CT_mat={}
for jj in np.arange(len(color_grid)):
	count_matrix = np.zeros((len(radii), len(coord_CLU)))
	for ii, rmax in enumerate(radii):
		counts = tree_array[jj].query_radius(coord_CLU, rmax, count_only=True)
		count_matrix[ii]=counts
		#print(ii, rmax)
	CT_mat[jj] = count_matrix.sum(axis=1)


coord_RR = deg_to_rad * np.transpose([ R2['DEC'], R2['RA'] ])

amin = 10.
radius_max = amin / 60. * n.pi / 180. # 10 arcminutes
asec = 10.
radius_step = asec / 3600. * n.pi / 180. # 10 arcsecond
radii = np.arange(radius_step, radius_max+radius_step, radius_step)

RR_mat={}
for jj in np.arange(len(color_grid)):
	count_matrix = np.zeros((len(radii), len(coord_RR)))
	for ii, rmax in enumerate(radii):
		counts = tree_array[jj].query_radius(coord_RR, rmax, count_only=True)
		count_matrix[ii]=counts
		#print(ii, rmax)
	RR_mat[jj] = count_matrix.sum(axis=1)

RES = {}
RES['N_D'] = len(D2)
RES['N_R'] = len(R2)
RES['z_bar'] = z_bar
RES['color_grid'] = color_grid
RES['color_step'] = color_step
RES['color_shift'] = color_shift
RES['theta_grid'] = radii
RES['theta_step'] = radius_step
RES['DATA'] = CT_mat
RES['RAND'] = RR_mat
##p_2_OUT  = os.path.join(topdir, "Counts_DATA_S0_RAND_gr_010_z_012.npy")
np.save(p_2_OUT, RES)
print(p_2_OUT, 'written')

#RES = np.load(p_2_OUT, allow_pickle='TRUE').item()
