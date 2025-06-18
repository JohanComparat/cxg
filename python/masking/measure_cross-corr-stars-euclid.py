import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
import numpy as n
import astropy.io.fits as fits
import os, sys, glob
from os.path import join
#import pymangle as mangle
import numpy as np
import matplotlib.pyplot as p
from scipy.interpolate import interp1d

import astropy.units as u
import astropy.cosmology as cc
cosmo = cc.Planck13
from astropy.coordinates import SkyCoord

from sklearn.neighbors import BallTree
from astropy.table import Table,unique,Column
from math import radians, cos, sin, asin, sqrt, pi

deg_to_rad = np.pi/180.
cat_dir = "/Users/Johan Comparat/Documents/Shared/data/Euclid/data/test_tile"
gaia_dir = "/Users/Johan Comparat/Documents/Shared/data/gaia_cat"
fig_dir  ='../../figures/mask'

mask_flag = int(sys.argv[2])
str_mask_flag = '_mask'+str(mask_flag)

d0 = Table.read( os.path.join( cat_dir, "EUC_LE3_CL-LE2-CATALOG_5E47CD_20250616T132858.719404Z_00.00.fits"  )             )
d_mask = Table.read( os.path.join( cat_dir, "EUC_LE3_CL-LE2-CATALOG_AC72D1_20250616T132855.875770Z_00.00.fits"  )             )
d_all = Table.read( os.path.join( cat_dir, "EUC_MER_FINAL-CAT_TILE102015162-8741AA_20250307T152950.857147Z_00.00.fits.gz" )  )
#os.path.join( cat_dir, "OutputCL.xml"     )
#os.path.join( cat_dir, "OutputCombCL.xml" )
if mask_flag == 1 :
	d1 = d_mask
if mask_flag == 0 :
	d1 = d_all
# galaxy list
E_sel = (d1['FLUX_VIS_2FWHM_APER']>0)#(d1['PHZ_MEDIAN']>=0)&
d1 = d1[E_sel]
ra_gal, dec_gal = d1['RIGHT_ASCENSION'], d1['DECLINATION']
gal_coordinates = deg_to_rad * np.transpose([dec_gal, ra_gal])
print('measures distances')
# Tree_obj_GAL = BallTree(gal_coordinates, metric='haversine')
ra_max = np.max(ra_gal)
ra_min = np.min(ra_gal)
dec_max = np.max(dec_gal)
dec_min = np.min(dec_gal)

option = sys.argv[1]
#option = "12_g_13"

# g0 = Table.read( os.path.join( gaia_dir, "table_5_g_6.fits"   ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_6_g_7.fits"   ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_7_g_8.fits"   ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_8_g_9.fits"   ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_9_g_10.fits"  ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_-10_g_5.fits" ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_10_g_11.fits" ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_11_g_12.fits" ) )
g0 = Table.read( os.path.join( gaia_dir, "table_" + option + ".fits" ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_13_g_14.fits" ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_14_g_15.fits" ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_15_g_16.fits" ) )
# g0 = Table.read( os.path.join( gaia_dir, "table_16_g_17.fits" ) )
#
s_field = (g0['ra']>=ra_min)&(g0['ra']<=ra_max)&(g0['dec']>=dec_min)&(g0['dec']<=dec_max)
ra_data  = g0['ra'][s_field]
dec_data = g0['dec'][s_field]
gaia_coordinates = deg_to_rad * np.transpose([dec_data, ra_data])
Tree_obj_STA = BallTree(gaia_coordinates , metric='haversine')

# random catalogue
RR_ra_all = []
RR_dec_all = []
N=0
while N<len(ra_data)*5:
	size = int(1e6)
	uu = np.random.uniform(size=size)
	dec_fs = np.arccos(1 - 2 * uu) * 180 / np.pi - 90.
	ra_fs = np.random.uniform(size=size) * 2 * np.pi * 180 / np.pi
	sR_field = (ra_fs>=ra_min)&(ra_fs<=ra_max)&(dec_fs>=dec_min)&(dec_fs<=dec_max)
	RR_ra_all.append(ra_fs[sR_field])
	RR_dec_all.append(dec_fs[sR_field])
	N+=len(ra_fs[sR_field])

RR = Table()
RR['RA'] = np.hstack(( RR_ra_all ))
RR['DEC'] = np.hstack(( RR_dec_all ))
coord_cat_RR  = deg_to_rad * n.transpose([RR['DEC'], RR['RA'] ])
Tree_obj_RR = BallTree(coord_cat_RR, metric='haversine')

N_gal = len(gal_coordinates)
N_gaia = len(gaia_coordinates)
N_RR = len(coord_cat_RR)

# event counts as a function of angular distance
arcsec = 1. / 3600.
rs = np.logspace(0, 3, 50)
r_rad = deg_to_rad * rs * arcsec

# CT bins
# n_bins = 10
n_bins = 10

CT_bins = np.logspace(-4, 4, n_bins)
CT_bins_min = CT_bins[:-1]
CT_bins_max = CT_bins[1:]

for c0, c1 in zip( CT_bins_min[::-1], CT_bins_max[::-1]):
	s1 = ( d1['FLUX_VIS_2FWHM_APER'] >= c0 ) & ( d1['FLUX_VIS_2FWHM_APER'] < c1 )
	N_data = len(gal_coordinates[s1])
	print(c0, c1, N_data)
	test_c = np.array([ Tree_obj_STA.query_radius(gal_coordinates[s1], r = rr, count_only=True) for rr in r_rad ])
	N_pairs_total = test_c.sum(axis=1)
	Delta_N_pairs = N_pairs_total[1:]-N_pairs_total[:-1]
	area = 4.*np.pi*(rs[1:]**2 - rs[:-1]**2) # arcsec**2
	pair_density = Delta_N_pairs/(area*(N_data*N_gaia)**0.5)
	prefix = str(np.round(c0, 5)) + '_FLUX_VIS_2FWHM_APER_' + str(np.round(c1, 5)) + str_mask_flag
	out_data = os.path.join(fig_dir , 'xcorr_gal_'+option+'_sourceCat_' + prefix + '.data')
	np.savetxt(out_data, np.transpose([rs[1:], pair_density]), header='radius_arcsec density' )
	#
	test_c = np.array([ Tree_obj_RR.query_radius(gal_coordinates[s1], r = rr, count_only=True) for rr in r_rad ])
	N_pairs_total = test_c.sum(axis=1)
	Delta_N_pairs = N_pairs_total[1:]-N_pairs_total[:-1]
	area = 4.*np.pi*(rs[1:]**2 - rs[:-1]**2) # arcsec**2
	pair_density = Delta_N_pairs/(area*(N_data*N_RR)**0.5)
	prefix = str(np.round(c0, 5)) + '_FLUX_VIS_2FWHM_APER_' + str(np.round(c1, 5)) + str_mask_flag
	out_data = os.path.join(fig_dir , 'xcorr_RR_'+option+'_sourceCat_' + prefix + '.data')
	np.savetxt(out_data, np.transpose([rs[1:], pair_density]), header='radius_arcsec density' )

