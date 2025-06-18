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

speed_light = cc.c.to(u.km/u.s).value
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
import Corrfunc
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from Corrfunc.io import read_catalog
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.theory.DD import DD
from numpy.random import choice

pimax = int(sys.argv[1])
topdir    = sys.argv[2]
p2_data   = sys.argv[3]
p2_random = sys.argv[4]
C_topdir    = sys.argv[5]
C_p2_data   = sys.argv[6]
C_p2_random = sys.argv[7]
PCF_basename  = sys.argv[8]
basename = sys.argv[8].split('-')[0]
Ms_min = float(sys.argv[9])
z_min = float(sys.argv[10])
z_max = float(sys.argv[11])
NSIDE_val = int(sys.argv[12])
out_subfolder = sys.argv[13]

JK_dir = os.path.join(topdir, out_subfolder)
os.system('mkdir -p '+JK_dir)
print(JK_dir)
#p_2_2PCF  = os.path.join(JK_dir, PCF_basename)


def tabulate_wprp_clustering_noW(RA, DEC, Z, rand_RA , rand_DEC, rand_Z, RA2, DEC2, Z2, rand_RA2 , rand_DEC2, rand_Z2, out_file='test.fits', pimax = 100.0, N_JK=100 ):
	CZ = Z * speed_light
	rand_CZ = rand_Z * speed_light
	CZ2 = Z2 * speed_light
	rand_CZ2 = rand_Z2 * speed_light
	N = len(RA)
	rand_N = len(rand_RA)
	N2 = len(RA2)
	rand_N2 = len(rand_RA2)
	#print(N, rand_N, out_file, time.time()-t0)
	bins = 10**np.arange(-2.0, 1.81, 0.1)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 32
	# Auto pairs counts in DD
	autocorr=0
	DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
								np.array(list(RA)).astype('float'),
								np.array(list(DEC)).astype('float'),
								np.array(list(CZ)).astype('float'),
								RA2= np.array(list(RA2)).astype('float'),
								DEC2=np.array(list(DEC2)).astype('float'),
								CZ2= np.array(list(CZ2)).astype('float')
								)#, is_comoving_dist=True)
	# Auto pairs counts in DR
	D1R_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							np.array(list(RA)).astype('float'),
							np.array(list(DEC)).astype('float'),
							np.array(list(CZ)).astype('float'),
							RA2 =rand_RA2.astype('float'),
							DEC2=rand_DEC2.astype('float'),
							CZ2 =rand_CZ2.astype('float'))
	D2R_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							np.array(list(RA2)).astype('float'),
							np.array(list(DEC2)).astype('float'),
							np.array(list(CZ2)).astype('float'),
							RA2 =rand_RA.astype('float'),
							DEC2=rand_DEC.astype('float'),
							CZ2 =rand_CZ.astype('float'))
	# Auto pairs counts in RR
	RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							rand_RA.astype('float'),
							rand_DEC.astype('float'),
							rand_CZ.astype('float'),
							RA2 =rand_RA2.astype('float'),
							DEC2=rand_DEC2.astype('float'),
							CZ2 =rand_CZ2.astype('float'))
	#print('RR',RR_counts['npairs'])
	# All the pair counts are done, get the angular correlation function
	wp = convert_rp_pi_counts_to_wp(N, N2, rand_N, rand_N2, DD_counts, D1R_counts, D2R_counts, RR_counts, nbins, pimax)
	t = Table()
	t.add_column(Column(data = bins[:-1], name='rp_min', unit='Mpc'  ) )
	t.add_column(Column(data = bins[1:], name='rp_max', unit='Mpc'  ) )
	t.add_column(Column(data = x, name='rp_mid', unit='Mpc'  ) )
	t.add_column(Column(data = wp, name='wprp', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * N, name='N_data', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * rand_N, name='N_random', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * pimax, name='pimax', unit=''  ) )
	# wprp_JK = np.zeros((N_JK, len(wp)))
	# for jj in np.arange(N_JK):
	# 		s1=(np.random.random(len(RA))<0.9)
	# 		s2=(np.random.random(len(RA2))<0.9)
	# 		autocorr=0
	# 		DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
	# 									np.array(list(RA[s1])).astype('float'),
	# 									np.array(list(DEC[s1])).astype('float'),
	# 									np.array(list(CZ[s1])).astype('float'),
	# 									RA2= np.array(list(RA2[s2])).astype('float'),
	# 									DEC2=np.array(list(DEC2[s2])).astype('float'),
	# 									CZ2= np.array(list(CZ2[s2])).astype('float')
	# 									)#, is_comoving_dist=True)
	# 		# Auto pairs counts in DR
	# 		D1R_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
	# 								np.array(list(RA[s1])).astype('float'),
	# 								np.array(list(DEC[s1])).astype('float'),
	# 								np.array(list(CZ[s1])).astype('float'),
	# 								RA2 =rand_RA2.astype('float'),
	# 								DEC2=rand_DEC2.astype('float'),
	# 								CZ2 =rand_CZ2.astype('float'))
	# 		D2R_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
	# 								np.array(list(RA2[s2])).astype('float'),
	# 								np.array(list(DEC2[s2])).astype('float'),
	# 								np.array(list(CZ2[s2])).astype('float'),
	# 								RA2 =rand_RA.astype('float'),
	# 								DEC2=rand_DEC.astype('float'),
	# 								CZ2 =rand_CZ.astype('float'))
	# 		# Auto pairs counts in RR
	# 		RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
	# 								rand_RA.astype('float'),
	# 								rand_DEC.astype('float'),
	# 								rand_CZ.astype('float'),
	# 								RA2 =rand_RA2.astype('float'),
	# 								DEC2=rand_DEC2.astype('float'),
	# 								CZ2 =rand_CZ2.astype('float'))
	# 		#print('RR',RR_counts['npairs'])
	# 		# All the pair counts are done, get the angular correlation function
	# 		wprp_JK[jj] = convert_rp_pi_counts_to_wp(N, N2, rand_N, rand_N2, DD_counts, D1R_counts, D2R_counts, RR_counts, nbins, pimax)
	# #t['wprp_JK'] = wprp_JK
	# t.add_column(Column(data = wprp_JK.mean(axis=0), name='wprp_JK_mean', unit=''  ) )
	# t.add_column(Column(data = wprp_JK.std(axis=0), name='wprp_JK_std', unit=''  ) )
	# print(out_file)
	t.write(out_file, overwrite=True, format='fits')
	print(out_file, time.time()-t0, 's')

data1 = Table.read(os.path.join(topdir, p2_data ), format='fits')
rand = Table.read(os.path.join(topdir, p2_random ), format='fits')

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

data_keep = ( data1['keep_HPX_04'] ) & ( data1['keep_HPX_08'] ) & ( data1['keep_HPX_16'] ) & ( data1['LPH_MASS_BEST'] > Ms_min) & (data1['BEST_Z']>=z_min) & (data1['BEST_Z']<=z_max)
rand_keep = ( rand['keep_HPX_04'] ) & ( rand['keep_HPX_08'] ) & ( rand['keep_HPX_16'] )

t0 = time.time()
print('before masking ND, NR = ', len(data1), len(rand))
DDD = data1[ (data_keep) ]
RRR = rand[  (rand_keep)  ]

N_gal = len(DDD)
N_RD = len(RRR)
RRR['Z'] = n.tile(DDD['BEST_Z'], int(N_RD*1./N_gal)+1)[:N_RD]

R_tF = n.random.random(N_RD)
N_R_F=5
sR = ( R_tF < N_gal * N_R_F / N_RD )
RRR = RRR[sR]

D2 = Table.read(os.path.join(C_topdir, C_p2_data ), format='fits')
R2 = Table.read(os.path.join(C_topdir, C_p2_random ), format='fits')
#print(D2.info())
#print(R2.info())
D2=D2[(D2['Z_LAMBDA']>=z_min)&(D2['Z_LAMBDA']<=z_max)]
R2=R2[(R2['redshift']>=z_min)&(R2['redshift']<=z_max)]

NSIDES = [2, 4, 8, 16, 32]
for NSIDE in NSIDES:
	NSIDE_str = str(NSIDE).zfill(2)
	DDD['HPX_'+NSIDE_str] = healpy.ang2pix(NSIDE, np.pi/2. - DDD['DEC']*np.pi/180. , DDD['RA']*np.pi/180. , nest=True)
	RRR['HPX_'+NSIDE_str]  = healpy.ang2pix(NSIDE, np.pi/2. - RRR['DEC']*np.pi/180. , RRR['RA']*np.pi/180. , nest=True)
	D2['HPX_'+NSIDE_str] = healpy.ang2pix(NSIDE, np.pi/2. - D2['DEC']*np.pi/180. , D2['RA']*np.pi/180. , nest=True)
	R2['HPX_'+NSIDE_str]  = healpy.ang2pix(NSIDE, np.pi/2. - R2['DEC']*np.pi/180. , R2['RA']*np.pi/180. , nest=True)

NSIDE = NSIDE_val
print('='*100)
print(NSIDE)
NSIDE_str = str(NSIDE).zfill(2)
for jk_i in np.arange(100):
	p_2_2PCF = os.path.join(JK_dir, PCF_basename + '_NSIDE_' + NSIDE_str + '_J_'+str(jk_i).zfill(4) + '.fits')
	print(jk_i, p_2_2PCF, os.path.isfile(p_2_2PCF))
	if os.path.isfile(p_2_2PCF) == False:
		U_NSIDE_list = np.unique(DDD['HPX_' + NSIDE_str])
		N_select = int(len(U_NSIDE_list)*0.9)
		pixels_keep = choice(U_NSIDE_list, size=N_select, replace=False)#, p=None)
		D_in = np.isin(DDD['HPX_' + NSIDE_str], pixels_keep)
		R_in = np.isin(RRR['HPX_' + NSIDE_str], pixels_keep)
		D2_in = np.isin(D2['HPX_' + NSIDE_str], pixels_keep)
		R2_in = np.isin(R2['HPX_' + NSIDE_str], pixels_keep)
		print(len(U_NSIDE_list), N_select, len(pixels_keep), pixels_keep[:15],
			  	len(DDD['RA'][D_in]), len(DDD['RA'][D_in])/len(DDD['RA'][D_in]), len(RRR['RA'][R_in])/len(RRR['RA'][R_in]),
				len(D2['RA'][D2_in]), len(D2['RA'][D2_in]) / len(D2['RA'][D2_in]), len(R2['RA'][R2_in]) / len(R2['RA'][R2_in]) )
		tabulate_wprp_clustering_noW(
			DDD['RA'][D_in], DDD['DEC'][D_in], DDD['BEST_Z'][D_in],
			RRR['RA'][R_in] , RRR['DEC'][R_in], RRR['Z'][R_in],
			D2['RA'][D2_in], D2['DEC'][D2_in], D2['Z_LAMBDA'][D2_in],
			R2['RA'][R2_in] , R2['DEC'][R2_in], R2['redshift'][R2_in],
			out_file=p_2_2PCF, pimax = 100.0)

