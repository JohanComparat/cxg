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
from sklearn.neighbors import BallTree

speed_light = cc.c.to(u.km/u.s).value
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
import Corrfunc
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from Corrfunc.io import read_catalog
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.theory.DD import DD

from scipy.interpolate import interp1d

pimax = int(sys.argv[1])
topdir    = sys.argv[2]
p2_data   = sys.argv[3]
p2_random = sys.argv[4]
basename = sys.argv[3][:-10]#.split('-')[0]
name_corr = sys.argv[5]#.split('-')[0]
Ms_min = float(sys.argv[6])

print('reading', os.path.join(os.environ['LSDR10'], 'sweep/MergeALL_BGSlike_LPH.fits'), time.time()-t0, 's')
GAL_i = Table.read(os.path.join(os.environ['LSDR10'], 'sweep/MergeALL_BGSlike_LPH.fits'))

m0 = float(basename.split('_')[3])
m1 = float(basename.split('_')[5])
z0 = float(basename.split('_')[6])
z1 = float(basename.split('_')[8])
GAL = GAL_i[(GAL_i['LPH_MASS_BEST']>Ms_min)&(GAL_i['LPH_MASS_BEST']<m1)&(GAL_i['Z_PHOT_MEAN']>z0)&(GAL_i['Z_PHOT_MEAN']<z1)]
print('opened', time.time()-t0, 's')
#ls10['BEST_Z'] = ls10['Z_PHOT_MEAN']
RS_model = Table.read( os.path.join( os.environ['GIT_STMOD_DATA'], 'data', 'models', 'model_GAL', 'legacy_dr10_south_v0.3_grz_z_cal_zspec_redgals_model.fit') )
z_RS = np.hstack(( 0., RS_model['nodes'][0] ))
gz_RS = np.hstack(( 1.3, RS_model['meancol'][0].sum(axis=1) ))
RS_color_gz = interp1d(z_RS, gz_RS)
#RS_color_gz = lambda redshift : redshift * 3 + 1.3
gz_med_RS = RS_color_gz(GAL['Z_PHOT_MEAN'])
#gz_min = gz_med_RS - 2 * scat
GAL['gz_med_RS'] = gz_med_RS
GAL['is_RS'] = ( GAL['g_mag']-GAL['z_mag']> GAL['gz_med_RS'] - 0.15 )
GAL['is_BC'] = ( GAL['g_mag']-GAL['z_mag']< GAL['gz_med_RS'] - 0.23 )
GAL['is_GV'] = (~GAL['is_RS'])&(~GAL['is_BC'])
GAL['UID'] = (n.round(GAL['RA'],6) * 10_000_000).astype('int64') * 1_000_000_000 + (n.round(GAL['DEC'],6) * 1_000_000).astype('int64')
print('RS, BC split', time.time()-t0, 's')


def tabulate_wprp_clustering_noW(RA, DEC, Z, rand_RA , rand_DEC, rand_Z, out_file='test.fits', CV_frac=0.01, pimax = 100.0, N_JK=100 ):
	"""
	wprp direct estimate
	path_2_data : path to the catalogue to correlate
	path_2_random
	"""
	#
	CZ = Z * speed_light
	rand_CZ = rand_Z * speed_light
	N = len(RA)
	rand_N = len(rand_RA)
	#print(N, rand_N, out_file, time.time()-t0)
	bins = 10**np.arange(-2.0, 1.81, 0.1)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 16
	# Auto pairs counts in DD
	autocorr=1
	DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
								np.array(list(RA)).astype('float'),
								np.array(list(DEC)).astype('float'),
								np.array(list(CZ)).astype('float') )#, is_comoving_dist=True)
	#print('DD',DD_counts['npairs'], DD_counts['npairs'].shape)
	# Auto pairs counts in DR
	autocorr=0
	DR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							np.array(list(RA)).astype('float'),
							np.array(list(DEC)).astype('float'),
							np.array(list(CZ)).astype('float'),
							RA2=rand_RA.astype('float'),
							DEC2=rand_DEC.astype('float'),
							CZ2=rand_CZ.astype('float'))
	#print('DR',DR_counts['npairs'])
	# Auto pairs counts in RR
	autocorr=1
	RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							rand_RA.astype('float'),
							rand_DEC.astype('float'),
							rand_CZ.astype('float'))
	#print('RR',RR_counts['npairs'])
	# All the pair counts are done, get the angular correlation function
	wp = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts, nbins, pimax)
	t = Table()
	t.add_column(Column(data = bins[:-1], name='rp_min', unit='Mpc'  ) )
	t.add_column(Column(data = bins[1:], name='rp_max', unit='Mpc'  ) )
	t.add_column(Column(data = x, name='rp_mid', unit='Mpc'  ) )
	t.add_column(Column(data = wp, name='wprp', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * N, name='N_data', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * rand_N, name='N_random', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * pimax, name='pimax', unit=''  ) )
	#t.add_column(Column(data = np.ones_like(x) * CV_frac, name='CV_frac', unit=''  ) )
	# repeat when removing random 10%
	wprp_JK = np.zeros((N_JK, len(wp)))
	for jj in np.arange(N_JK):
		s1=(np.random.random(len(RA))<0.9)
		ra1, dec1, cz1 = np.array(list(RA)).astype('float')[s1], np.array(list(DEC)).astype('float')[s1], np.array(list(CZ)).astype('float')[s1]
		# Auto pairs counts in DD
		autocorr=1
		DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins, ra1, dec1, cz1)
		# Auto pairs counts in DR
		autocorr=0
		DR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins, ra1, dec1, cz1, RA2=rand_RA.astype('float'), DEC2=rand_DEC.astype('float'), CZ2=rand_CZ.astype('float'))
		# Auto pairs counts in RR
		autocorr=1
		RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins, rand_RA.astype('float'), rand_DEC.astype('float'), rand_CZ.astype('float'))
		# All the pair counts are done, get the angular correlation function
		wprp_JK[jj] = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts, nbins, pimax)
	#print ( "wprp_JK.mean(axis=0)", wprp_JK.mean(axis=0).shape, wprp_JK.mean(axis=0),  wprp_JK.mean(axis=0) / wp  )
	#print ( "wprp_JK.std(axis=0) ", wprp_JK.std(axis=0) .shape, wprp_JK.std(axis=0),  wprp_JK.std(axis=0) / wp  )
	#print ( "wp                ", wp                .shape, wp                  )
	#t['wprp_JK'] = wprp_JK
	t.add_column(Column(data = wprp_JK.mean(axis=0), name='wprp_JK_mean', unit=''  ) )
	t.add_column(Column(data = wprp_JK.std(axis=0), name='wprp_JK_std', unit=''  ) )
	#print(out_file)
	t.write(out_file, overwrite=True, format='fits')
	print(out_file, time.time()-t0, 's')

print('opening data, randoms', time.time()-t0, 's')
data1 = Table.read(os.path.join(topdir, p2_data ), format='fits')
rand = Table.read(os.path.join(topdir, p2_random ), format='fits')

print('filtering pixels', time.time()-t0, 's')
NSIDES = [4, 8, 16]
for NSIDE in NSIDES:
	NSIDE_str = str(NSIDE).zfill(2)
	data1['HPX_'+NSIDE_str] = healpy.ang2pix(NSIDE, np.pi/2. - data1['DEC']*np.pi/180. , data1['RA']*np.pi/180. , nest=True)
	rand['HPX_'+NSIDE_str]  = healpy.ang2pix(NSIDE, np.pi/2. - rand['DEC']*np.pi/180. , rand['RA']*np.pi/180. , nest=True)
	mask_dir = os.path.join(os.environ['GIT_DR10W_DATA'], 'figures/BGS_wprp/HPX_'+NSIDE_str)
	REJ_file = os.path.join( mask_dir, basename+'-wprp-pimax100-bin0p05-JK100.fits-wprp-all-bin0p05_REJoutlierClipped3Sigma.fits')
	if os.path.isfile(REJ_file) :
		t = Table.read( REJ_file )
		data1['keep_HPX_'+NSIDE_str] = (np.isin(data1['HPX_'+NSIDE_str], t['hpx_idx'], invert = True))
		rand['keep_HPX_'+NSIDE_str]  = (np.isin(rand['HPX_'+NSIDE_str], t['hpx_idx'], invert = True))
	else:
		data1['keep_HPX_'+NSIDE_str] = True
		rand['keep_HPX_'+NSIDE_str]  = True

data_keep = ( data1['keep_HPX_04'] ) & ( data1['keep_HPX_08'] ) & ( data1['keep_HPX_16'] ) & ( data1['LPH_MASS_BEST'] > Ms_min)
rand_keep = ( rand['keep_HPX_04'] ) & ( rand['keep_HPX_08'] ) & ( rand['keep_HPX_16'] )

t0 = time.time()
print('before masking ND, NR = ', len(data1), len(rand), time.time()-t0, 's')
DDD = data1[ (data_keep) ]
RRR = rand[  (rand_keep)  ]
print('full set ready', time.time()-t0, 's')


deg_to_rad = np.pi/180.
coord_GAL_BC = deg_to_rad * np.transpose([GAL['DEC'][GAL['is_BC']], GAL['RA'][GAL['is_BC']] ])
Tree_GAL_BC = BallTree(coord_GAL_BC, metric='haversine')
coord_GAL_RS = deg_to_rad * np.transpose([GAL['DEC'][GAL['is_RS']], GAL['RA'][GAL['is_RS']] ])
Tree_GAL_RS = BallTree(coord_GAL_RS, metric='haversine')
print('trees constructed', time.time()-t0, 's')

coord_GAL = deg_to_rad * np.transpose([ DDD['DEC'], DDD['RA'] ])

asec = 0.01
radius = asec / 3600. * n.pi / 180. # 0.1 arcsecond
dist_BC, ind_BC = Tree_GAL_BC.query(coord_GAL, k=1, return_distance=True)
print('Tree_GAL_BC.query', time.time()-t0, 's')
dist_RS, ind_RS = Tree_GAL_RS.query(coord_GAL, k=1, return_distance=True)
print('Tree_GAL_RS.query', time.time()-t0, 's')

is_BC_4_2pcf = np.arange(len(coord_GAL))[(np.hstack((dist_BC))<radius)]
is_RS_4_2pcf = np.arange(len(coord_GAL))[(np.hstack((dist_RS))<radius)]

N_gal = len(DDD[is_BC_4_2pcf])
N_RD = len(RRR)
RRR['Z'] = n.tile(DDD['BEST_Z'][is_BC_4_2pcf], int(N_RD*1./N_gal)+1)[:N_RD]
R_tF = n.random.random(N_RD)
N_R_F=5
sR = ( R_tF < N_gal * N_R_F / N_RD )
RRR = RRR[sR]

p_2_2PCF  = os.path.join(topdir, 'BC_' + name_corr)
p_2_DATA_OUT  = os.path.join(topdir, 'BC_' + p2_data)
p_2_RAND_OUT  = os.path.join(topdir, 'BC_' + p2_random)
print('BC after masking & Mmin selection ND, NR = ', len(DDD[is_BC_4_2pcf]), len(RRR), time.time()-t0, 's')
if os.path.isfile(p_2_2PCF)==False :
	try:
		tabulate_wprp_clustering_noW(DDD['RA'][is_BC_4_2pcf], DDD['DEC'][is_BC_4_2pcf], DDD['BEST_Z'][is_BC_4_2pcf], RRR['RA'] , RRR['DEC'], RRR['Z'],  out_file=p_2_2PCF, CV_frac=0.01, pimax = 100.0, N_JK = 2 )
	except(RuntimeError):
		print('RuntimeError')

print('before masking ND, NR = ', len(data1), len(rand), time.time()-t0, 's')
DDD = data1[ (data_keep) ]
RRR = rand[  (rand_keep)  ]
print('full set ready', time.time()-t0, 's')

DDD[is_BC_4_2pcf].write(p_2_DATA_OUT, overwrite = True)
print(p_2_DATA_OUT, 'written', time.time()-t0, 's')
RRR.write(p_2_RAND_OUT, overwrite = True)
print(p_2_RAND_OUT, 'written', time.time()-t0, 's')
RRR = rand[  (rand_keep)  ]

N_gal = len(DDD[is_RS_4_2pcf])
N_RD = len(RRR)
RRR['Z'] = n.tile(DDD['BEST_Z'][is_RS_4_2pcf], int(N_RD*1./N_gal)+1)[:N_RD]
R_tF = n.random.random(N_RD)
N_R_F=5
sR = ( R_tF < N_gal * N_R_F / N_RD )
RRR = RRR[sR]

p_2_2PCF  = os.path.join(topdir, 'RS_' + name_corr)
p_2_DATA_OUT  = os.path.join(topdir, 'RS_' + p2_data)
p_2_RAND_OUT  = os.path.join(topdir, 'RS_' + p2_random)
print('RS after masking & Mmin selection ND, NR = ', len(DDD[is_RS_4_2pcf]), len(RRR), time.time()-t0, 's')
if os.path.isfile(p_2_2PCF)==False :
	try:
		tabulate_wprp_clustering_noW(DDD['RA'][is_RS_4_2pcf], DDD['DEC'][is_RS_4_2pcf], DDD['BEST_Z'][is_RS_4_2pcf], RRR['RA'] , RRR['DEC'], RRR['Z'],  out_file=p_2_2PCF, CV_frac=0.01, pimax = 100.0, N_JK = 2 )
	except(RuntimeError):
		print('RuntimeError')

DDD[is_RS_4_2pcf].write(p_2_DATA_OUT, overwrite = True)
print(p_2_DATA_OUT, 'written', time.time()-t0, 's')
RRR.write(p_2_RAND_OUT, overwrite = True)
print(p_2_RAND_OUT, 'written', time.time()-t0, 's')
