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
deg_to_rad = np.pi/180.

topdir    = sys.argv[1]
p2_data   = sys.argv[2]
C_topdir    = sys.argv[3]
C_p2_data   = sys.argv[4]
C_p2_random = sys.argv[5]
p_2_OUT  = os.path.join(topdir, sys.argv[6])
if os.path.isfile(p_2_OUT)==False:
	basename = sys.argv[4][:-4]
	color_name = sys.argv[7]
	z_min = float(sys.argv[8])
	z_max = float(sys.argv[9])

	D2 = Table.read(os.path.join(C_topdir, C_p2_data ), format='fits')
	R2 = Table.read(os.path.join(C_topdir, C_p2_random ), format='fits')
	print(len(D2), len(R2), 'full')
	D2=D2[(D2['Z_LAMBDA']>=z_min)&(D2['Z_LAMBDA']<=z_max)]
	R2=R2[(R2['redshift']>=z_min)&(R2['redshift']<=z_max)]
	print(len(D2), len(R2), 'zcut')

	coord_D2 = deg_to_rad * np.transpose([D2['DEC'], D2['RA'] ])
	Tree_D2 = BallTree(coord_D2, metric='haversine')
	coord_R2 = deg_to_rad * np.transpose([R2['DEC'], R2['RA'] ])
	Tree_R2 = BallTree(coord_R2, metric='haversine')

	print('reading', os.path.join(topdir, p2_data), time.time()-t0, 's')
	GAL_i = Table.read(os.path.join(topdir, p2_data))
	#keep = (t["VIS_DET"]==1)&(t["DET_QUALITY_FLAG"]==0) & (t["FLUX_VIS_PSF"]/t["FLUXERR_VIS_PSF"]>=5) & (t["FLAG_VIS"]==0) & (t["FLAG_Y"]==0)&(t["FLAG_G_EXT_DECAM"]==0)&(t["FLAG_R_EXT_DECAM"]==0) & ( -2.5*np.log10(t["FLUX_G_EXT_DECAM_2FWHM_APER"]/t["FLUX_R_EXT_DECAM_2FWHM_APER"]) > -0.5 ) & ( -2.5*np.log10(t["FLUX_VIS_2FWHM_APER"]/t["FLUX_Y_2FWHM_APER"]) > -0.5 )  & ( -2.5*np.log10(t["FLUX_G_EXT_DECAM_2FWHM_APER"]/t["FLUX_R_EXT_DECAM_2FWHM_APER"]) < 3 ) & ( -2.5*np.log10(t["FLUX_VIS_2FWHM_APER"]/t["FLUX_Y_2FWHM_APER"]) < 3 )
	# & (t["FLUX_Y_2FWHM_APER"]/t["FLUXERR_Y_2FWHM_APER"]>=5) & (t["FLUX_G_EXT_DECAM_2FWHM_APER"]/t["FLUXERR_G_EXT_DECAM_2FWHM_APER"]>=5) & (t["FLUX_R_EXT_DECAM_2FWHM_APER"]/t["FLUXERR_R_EXT_DECAM_2FWHM_APER"]>=5) & ( -2.5*np.log10(t["FLUXERR_G_EXT_DECAM_2FWHM_APER"]/t["FLUXERR_R_EXT_DECAM_2FWHM_APER"]) > -0.5 ) & ( -2.5*np.log10(t["FLUX_VIS_2FWHM_APER"]/t["FLUX_Y_2FWHM_APER"]) > -0.5 )  & ( -2.5*np.log10(t["FLUXERR_G_EXT_DECAM_2FWHM_APER"]/t["FLUXERR_R_EXT_DECAM_2FWHM_APER"]) < 3 ) & ( -2.5*np.log10(t["FLUX_VIS_2FWHM_APER"]/t["FLUX_Y_2FWHM_APER"]) < 3 ) #
	#GAL_i = t[keep]
	#GAL_i['gr']=-2.5*np.log10(GAL_i["FLUX_G_EXT_DECAM_2FWHM_APER"]/GAL_i["FLUX_R_EXT_DECAM_2FWHM_APER"])
	#GAL_i['vy']=-2.5*np.log10(GAL_i["FLUX_VIS_2FWHM_APER"]/GAL_i["FLUX_Y_2FWHM_APER"])
	#GAL_i.keep_columns(['RIGHT_ASCENSION', 'DECLINATION', 'gr', 'vy'])
	#GAL_i.write(os.path.join(topdir, 'concat_short.fits'))
	print(len(GAL_i))
	color = GAL_i[color_name]
	coord_GAL = deg_to_rad * np.transpose([GAL_i['DECLINATION'], GAL_i['RIGHT_ASCENSION'] ])

	dd, ii = Tree_R2.query(coord_GAL, return_distance=True)
	has_gal_1amin = np.unique( np.hstack((ii)) [( np.hstack((dd))/deg_to_rad < 1/60. )] )
	R2 = R2[has_gal_1amin]

	dd, ii = Tree_D2.query(coord_GAL, return_distance=True)
	has_gal_1amin = np.unique( np.hstack((ii)) [( np.hstack((dd))/deg_to_rad < 1/60. )] )
	D2 = D2[has_gal_1amin]
	print(len(D2), len(R2), 'euclid cut')
	z_bar = np.mean(D2['Z_LAMBDA'])
	if len(D2)>1 and len(R2)>1:

		color_step = 0.05
		color_grid = np.arange(np.round(color.min(),2), np.round(color.max(),2)-color_step, color_step)
		color_bins = np.arange(np.round(color.min(),2), np.round(color.max(),2), color_step)
		print('color_grid', color_grid)
		N_gal_in_bin = np.histogram(color, bins=color_bins)[0]
		print('N_gal_in_bin', N_gal_in_bin)
		color_grid = color_grid[N_gal_in_bin>10]
		print('color_grid', color_grid)
		tree_array = {}
		for jj, c_min in enumerate(color_grid):
			selection = (color>=c_min)&(color<c_min+color_step)
			coord_GAL = deg_to_rad * np.transpose([GAL_i['DECLINATION'][selection], GAL_i['RIGHT_ASCENSION'][selection] ])
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
		RES['theta_grid'] = radii
		RES['theta_step'] = radius_step
		RES['DATA'] = CT_mat
		RES['RAND'] = RR_mat
		##p_2_OUT  = os.path.join(topdir, "Counts_DATA_S0_RAND_gr_010_z_012.npy")
		np.save(p_2_OUT, RES)
		print(p_2_OUT, 'written')

		#RES = np.load(p_2_OUT, allow_pickle='TRUE').item()
