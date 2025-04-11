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
import numpy as np
from Corrfunc.theory.wp import wp
from Corrfunc.theory.xi import xi
from Corrfunc.io import read_catalog

# Simple python script to read L-Galaxies output (Ayromlou+ 2021b)
# In[]
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5py
h_const = 0.6774

def hdf5_to_dict(hdf5_group, print_header=False):
    result = {}
    for key in hdf5_group.keys():
        item = hdf5_group[key]
        if isinstance(item, h5py.Dataset):
            result[key] = item[()]
        elif isinstance(item, h5py.Group):
            result[key] = hdf5_to_dict(item)
        if print_header and key == 'Header':
            print("Header:")
            for attr_name, attr_value in item.attrs.items():
                print(f"{attr_name}: {attr_value}")
    return result

def open_hdf5_dict(file_path, print_header=False):
    with h5py.File(file_path, 'r') as hdf:
        return hdf5_to_dict(hdf, print_header)
filename = '/media/sf_Shared/data/LGAL/LGal2021_galaxies_z0p2.hdf5'
filename = '/data36s/simulation_1/LGal2021_galaxies_z0p2.hdf5'
filename_cluster = '/data36s/simulation_1/LGal2021_BCG_z0p2.hdf5'
gal = open_hdf5_dict(filename)
#cluster = open_hdf5_dict(filename_cluster)

#Millennium (Planck1) 0.315 0.049 0.685 67.3 0.96 0.826 21603 1.43 × 109 714 64 56
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.3 * u.km / u.s / u.Mpc, Om0=0.315)
h = 0.673
L_box = 714.0
z_array = n.arange(0, 7.5, 0.001)
dc_to_z = interp1d(cosmo.comoving_distance(z_array), z_array)


def get_radecz(x,y,z,vxn, vyn, vzn):
	rr = (x**2 + y**2 + z**2)**0.5
	theta = n.arccos(z / rr) * 180 / n.pi
	phi = n.arctan2(y, x) * 180 / n.pi
	ra = phi + 180.
	dec = theta - 90.
	redshift_R = dc_to_z(rr)
	vPara = (vxn * x + vyn * y + vzn * z) / rr
	rr_s = rr + vPara / cosmo.H(redshift_R).value
	rr_s[rr_s <= 0] = rr[rr_s <= 0]
	redshift_S = dc_to_z(rr_s)
	return ra, dec, redshift_S

def get_radecz_R(x,y,z):
	rr = (x**2 + y**2 + z**2)**0.5
	theta = n.arccos(z / rr) * 180 / n.pi
	phi = n.arctan2(y, x) * 180 / n.pi
	ra = phi + 180.
	dec = theta - 90.
	redshift_R = dc_to_z(rr)
	return ra, dec, redshift_R



Xg=gal['Pos'].T[0]/1000.
Yg=gal['Pos'].T[1]/1000.
Zg=gal['Pos'].T[2]/1000.
VXg=gal['Vel'].T[0]
VYg=gal['Vel'].T[1]
VZg=gal['Vel'].T[2]
ra_g, dec_g, redshift_g = get_radecz(Xg, Yg, Zg, VXg, VYg, VZg)

def tabulate_wprp_clustering_noW(RA, DEC, Z, rand_RA , rand_DEC, rand_Z, RA2, DEC2, Z2, rand_RA2 , rand_DEC2, rand_Z2, out_file='test.fits', pimax = 100.0, N_JK=100 ):
	CZ = Z * speed_light
	rand_CZ = rand_Z * speed_light
	CZ2 = Z2 * speed_light
	rand_CZ2 = rand_Z2 * speed_light
	N = len(RA)
	rand_N = len(rand_RA)
	N2 = len(RA2)
	rand_N2 = len(rand_RA2)
	print(N, rand_N, N2, rand_N2, out_file, time.time()-t0)
	bins = 10**np.arange(-2.0, 1.6, 0.25)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 16
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
	t.add_column(Column(data = np.ones_like(x) * N2, name='N_data2', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * rand_N, name='N_random', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * rand_N2, name='N_random2', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * pimax, name='pimax', unit=''  ) )
	wprp_JK = np.zeros((N_JK, len(wp)))
	for jj in np.arange(N_JK):
			s1=(np.random.random(len(RA))<0.9)
			s2=(np.random.random(len(RA2))<0.9)
			autocorr=0
			DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
										np.array(list(RA[s1])).astype('float'),
										np.array(list(DEC[s1])).astype('float'),
										np.array(list(CZ[s1])).astype('float'),
										RA2= np.array(list(RA2[s2])).astype('float'),
										DEC2=np.array(list(DEC2[s2])).astype('float'),
										CZ2= np.array(list(CZ2[s2])).astype('float')
										)#, is_comoving_dist=True)
			# Auto pairs counts in DR
			D1R_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
									np.array(list(RA[s1])).astype('float'),
									np.array(list(DEC[s1])).astype('float'),
									np.array(list(CZ[s1])).astype('float'),
									RA2 =rand_RA2.astype('float'),
									DEC2=rand_DEC2.astype('float'),
									CZ2 =rand_CZ2.astype('float'))
			D2R_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
									np.array(list(RA2[s2])).astype('float'),
									np.array(list(DEC2[s2])).astype('float'),
									np.array(list(CZ2[s2])).astype('float'),
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
			wprp_JK[jj] = convert_rp_pi_counts_to_wp(N, N2, rand_N, rand_N2, DD_counts, D1R_counts, D2R_counts, RR_counts, nbins, pimax)
	#t['wprp_JK'] = wprp_JK
	t.add_column(Column(data = wprp_JK.mean(axis=0), name='wprp_JK_mean', unit=''  ) )
	t.add_column(Column(data = wprp_JK.std(axis=0), name='wprp_JK_std', unit=''  ) )
	#print(out_file)
	t.write(out_file, overwrite=True, format='fits')
	print(out_file, time.time()-t0, 's')

density_S1=2993/13116 # per square degree
density_S1_gpc3=12.85*10**(-7) # per square degree
density_G1075=164.8 # per square degree
density_G1075_gpc3=79.7*10**(-5)
area_sim = 129600./(8*np.pi)
volume_sim = L_box**3

for M_min in np.arange(13.5, 15.6, 0.1):
	s_g_t = (gal['M_Crit200']>10**M_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	s_g_RF = (gal['M_Crit200']>10**M_min)&(gal['Type']==0)
	print(np.round(M_min,2), len(ra_g[s_g_t]),
	   #', sky density', np.round(len(ra_g[s_g_t])/area_sim, 3), np.round(density_S1,3 ), 'deg-2,',
	   ', rho', np.round(np.log10(len(ra_g[s_g_RF])/volume_sim),2), np.round(np.log10(density_S1_gpc3),2), 'Mpc-3',
	   ', R_Crit200=', np.round(gal['R_Crit200'][s_g_RF].mean(),1),
	   ', StellarMass=', np.round(np.log10(gal['StellarMass'][s_g_RF].mean()),1) )

for MS_min in np.arange(10, 11.5, 0.01):
	s_g_t = (gal['StellarMass']+gal['HaloStellarMass']>10**MS_min)&(redshift_g>=0.1)&(redshift_g<=0.3)
	s_g_RF = (gal['StellarMass']+gal['HaloStellarMass']>10**MS_min)
	print(np.round(MS_min,2), len(ra_g[s_g_t]),
	   #', sky density', np.round(len(ra_g[s_g_t])/area_sim, 3), np.round(density_G1075,3 ), 'deg-2,',
	   ', rho', np.round(np.log10(len(ra_g[s_g_RF])/volume_sim),3), np.round(np.log10(density_G1075_gpc3),3), 'Mpc-3' )

pcf_dir ='/home/comparat/software/cxg/data/lgal/'

Ms_val_GAM = 10.67

s_g = (gal['StellarMass']>10**Ms_val_GAM)&(redshift_g>=0.1)&(redshift_g<=0.3)
RR_x, RR_y, RR_z = np.random.uniform(0, L_box, size=(3,5*len(ra_g[s_g])))
RR_ra_g, RR_dec_g, RR_redshift_g = get_radecz_R(RR_x, RR_y, RR_z)
RR_g_sel = (RR_redshift_g>=0.1)&(RR_redshift_g<=0.3)
for M_min in np.arange(14., 15., 0.1):
	p_2_2PCF = os.path.join(pcf_dir, 'LGAL_01z03_M200c_'+str(np.round(M_min,2))+'_Ms_1067_JK_100.fits')
	s_c = (gal['M_Crit200']>10**M_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	RRc_x, RRc_y, RRc_z = np.random.uniform(0, L_box, size=(3,20*len(ra_g[s_c])))
	RR_ra_c, RR_dec_c, RR_redshift_c = get_radecz_R(RRc_x, RRc_y, RRc_z)
	RR_c_sel = (RR_redshift_c>=0.1)&(RR_redshift_c<=0.3)
	tabulate_wprp_clustering_noW(
		ra_g[s_g], dec_g[s_g], redshift_g[s_g], RR_ra_g[RR_g_sel], RR_dec_g[RR_g_sel], RR_redshift_g[RR_g_sel],
		ra_g[s_c], dec_g[s_c], redshift_g[s_c], RR_ra_c[RR_c_sel], RR_dec_c[RR_c_sel], RR_redshift_c[RR_c_sel],
		out_file=p_2_2PCF, pimax = 100.0, N_JK=100)

s_g = (gal['StellarMass']>10**Ms_val_GAM)&(redshift_g>=0.1)&(redshift_g<=0.3)&(np.log10(gal['StarFormationRate']/gal['StellarMass'])<-11)
RR_x, RR_y, RR_z = np.random.uniform(0, L_box, size=(3,5*len(ra_g[s_g])))
RR_ra_g, RR_dec_g, RR_redshift_g = get_radecz_R(RR_x, RR_y, RR_z)
RR_g_sel = (RR_redshift_g>=0.1)&(RR_redshift_g<=0.3)
for M_min in np.arange(14., 15., 0.1):
	p_2_2PCF = os.path.join(pcf_dir, 'LGAL_01z03_M200c_'+str(np.round(M_min,2))+'_RS_1067_JK_100.fits')
	s_c = (gal['M_Crit200']>10**M_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	RRc_x, RRc_y, RRc_z = np.random.uniform(0, L_box, size=(3,20*len(ra_g[s_c])))
	RR_ra_c, RR_dec_c, RR_redshift_c = get_radecz_R(RRc_x, RRc_y, RRc_z)
	RR_c_sel = (RR_redshift_c>=0.1)&(RR_redshift_c<=0.3)
	tabulate_wprp_clustering_noW(
		ra_g[s_g], dec_g[s_g], redshift_g[s_g], RR_ra_g[RR_g_sel], RR_dec_g[RR_g_sel], RR_redshift_g[RR_g_sel],
		ra_g[s_c], dec_g[s_c], redshift_g[s_c], RR_ra_c[RR_c_sel], RR_dec_c[RR_c_sel], RR_redshift_c[RR_c_sel],
		out_file=p_2_2PCF, pimax = 100.0, N_JK=100)


s_g = (gal['StellarMass']>10**Ms_val_GAM)&(redshift_g>=0.1)&(redshift_g<=0.3)&(np.log10(gal['StarFormationRate']/gal['StellarMass'])>=-11)
RR_x, RR_y, RR_z = np.random.uniform(0, L_box, size=(3,5*len(ra_g[s_g])))
RR_ra_g, RR_dec_g, RR_redshift_g = get_radecz_R(RR_x, RR_y, RR_z)
RR_g_sel = (RR_redshift_g>=0.1)&(RR_redshift_g<=0.3)
for M_min in np.arange(14., 15., 0.1):
	p_2_2PCF = os.path.join(pcf_dir, 'LGAL_01z03_M200c_'+str(np.round(M_min,2))+'_BC_1067_JK_100.fits')
	s_c = (gal['M_Crit200']>10**M_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	RRc_x, RRc_y, RRc_z = np.random.uniform(0, L_box, size=(3,20*len(ra_g[s_c])))
	RR_ra_c, RR_dec_c, RR_redshift_c = get_radecz_R(RRc_x, RRc_y, RRc_z)
	RR_c_sel = (RR_redshift_c>=0.1)&(RR_redshift_c<=0.3)
	tabulate_wprp_clustering_noW(
		ra_g[s_g], dec_g[s_g], redshift_g[s_g], RR_ra_g[RR_g_sel], RR_dec_g[RR_g_sel], RR_redshift_g[RR_g_sel],
		ra_g[s_c], dec_g[s_c], redshift_g[s_c], RR_ra_c[RR_c_sel], RR_dec_c[RR_c_sel], RR_redshift_c[RR_c_sel],
		out_file=p_2_2PCF, pimax = 100.0, N_JK=100)


for LX_min in np.arange(42.7, 45.6, 0.05):
	s_c_t = (gal['XrayLum']>LX_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	s_c_RF = (gal['XrayLum']>LX_min)&(gal['Type']==0)
	print(np.round(LX_min,2), len(ra_g[s_c_t]),
	   #', sky density', np.round(len(ra_g[s_c_t])/area_sim, 3), np.round(density_S1,3 ), 'deg-2,',
	   ', rho', np.round(np.log10(len(ra_g[s_c_RF])/volume_sim),2), np.round(np.log10(density_S1_gpc3),2), 'Mpc-3',
	   ', R_Crit200=', np.round(gal['R_Crit200'][s_c_RF].mean(),1),
	   ', M_Crit200=', np.round(np.log10(gal['M_Crit200'][s_c_RF].mean()),1),
	   ', StellarMass=', np.round(np.log10(gal['StellarMass'][s_c_RF].mean()),1)
	   )


from colossus.cosmology import cosmology
cosmology.setCosmology('planck18')
from colossus.halo import mass_so
from colossus.halo import mass_defs
from colossus.halo import concentration

c200c = concentration.concentration(gal['M_Crit200']*h, '200c', 0.2, model = 'ishiyama21')
M500coh, R500coh, c500coh = mass_defs.changeMassDefinition(gal['M_Crit200']*h, c200c, 0.2, '200c', '500c')
M_Crit500 = M500coh / h

fun_M_L =  lambda log10M500c : 44.7 + 1.61 * (log10M500c-15)
sigma_LX = 0.3
fun_M_T =  lambda log10M500c : 0.6 * log10M500c - 8.
sigma_kT = 0.2
# covariance
rho = 0.95
cov_kT_LX = np.array([[sigma_kT**2, rho*sigma_kT*sigma_LX ],[rho*sigma_kT*sigma_LX, sigma_LX**2]])
# generates means
EZ = cosmo.efunc(redshift_g)
kT_Mean_oEzm23 = ( fun_M_T(np.log10(M_Crit500)) ) #+ 2./3. * np.log10(EZ)
LX_Mean_eZm2 = fun_M_L( np.log10(M_Crit500) ) #+ 2*np.log10(EZ)
# generate values with correlated scatter
corr_scat = np.random.multivariate_normal([0,0], cov_kT_LX, size=len(M_Crit500))
corr_scat_kT = corr_scat.T[0]
corr_scat_LX = corr_scat.T[1]
CLUSTER_kT = 10**( kT_Mean_oEzm23 + corr_scat_kT + 2./3. * np.log10(EZ) )
CLUSTER_LX_soft_RF_R500c = LX_Mean_eZm2 + corr_scat_LX + 2*np.log10(EZ)


for LX_min in np.arange(42.7, 45.6, 0.05):
	s_c_t = (CLUSTER_LX_soft_RF_R500c>LX_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	s_c_RF = (CLUSTER_LX_soft_RF_R500c>LX_min)&(gal['Type']==0)
	print(np.round(LX_min,2), len(ra_g[s_c_t]),
	   #', sky density', np.round(len(ra_g[s_c_t])/area_sim, 3), np.round(density_S1,3 ), 'deg-2,',
	   ', rho', np.round(np.log10(len(ra_g[s_c_RF])/volume_sim),2), np.round(np.log10(density_S1_gpc3),2), 'Mpc-3',
	   ', R_Crit200=', np.round(gal['R_Crit200'][s_c_RF].mean(),1),
	   ', M_Crit200=', np.round(np.log10(gal['M_Crit200'][s_c_RF].mean()),1),
	   ', StellarMass=', np.round(np.log10(gal['StellarMass'][s_c_RF].mean()),1)
	   )

for MS_min in np.arange(10, 11.5, 0.01):
	s_g_t = (gal['StellarMass']>10**MS_min)&(redshift_g>=0.1)&(redshift_g<=0.3)
	s_g_RF = (gal['StellarMass']>10**MS_min)
	print(np.round(MS_min,2), len(ra_g[s_g_t]),
	   #', sky density', np.round(len(ra_g[s_g_t])/area_sim, 3), np.round(density_G1075,3 ), 'deg-2,',
	   ', rho', np.round(np.log10(len(ra_g[s_g_RF])/volume_sim),2), np.round(np.log10(density_G1075_gpc3),2), 'Mpc-3' )

pcf_dir ='/home/comparat/software/cxg/data/lgal/'

Ms_val_GAM = 10.67

s_g = (gal['StellarMass']>10**Ms_val_GAM)&(redshift_g>=0.1)&(redshift_g<=0.3)
RR_x, RR_y, RR_z = np.random.uniform(0, L_box, size=(3,5*len(ra_g[s_g])))
RR_ra_g, RR_dec_g, RR_redshift_g = get_radecz_R(RR_x, RR_y, RR_z)
RR_g_sel = (RR_redshift_g>=0.1)&(RR_redshift_g<=0.3)
for LX_min in np.arange(42.7, 43.6, 0.1):
	p_2_2PCF = os.path.join(pcf_dir, 'LGAL_01z03_L0520_'+str(np.round(LX_min,2))+'_Ms_1067_JK_100.fits')
	s_c = (CLUSTER_LX_soft_RF_R500c>LX_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	RRc_x, RRc_y, RRc_z = np.random.uniform(0, L_box, size=(3,20*len(ra_g[s_c])))
	RR_ra_c, RR_dec_c, RR_redshift_c = get_radecz_R(RRc_x, RRc_y, RRc_z)
	RR_c_sel = (RR_redshift_c>=0.1)&(RR_redshift_c<=0.3)
	tabulate_wprp_clustering_noW(
		ra_g[s_g], dec_g[s_g], redshift_g[s_g], RR_ra_g[RR_g_sel], RR_dec_g[RR_g_sel], RR_redshift_g[RR_g_sel],
		ra_g[s_c], dec_g[s_c], redshift_g[s_c], RR_ra_c[RR_c_sel], RR_dec_c[RR_c_sel], RR_redshift_c[RR_c_sel],
		out_file=p_2_2PCF, pimax = 100.0, N_JK=100)

s_g = (gal['StellarMass']>10**Ms_val_GAM)&(redshift_g>=0.1)&(redshift_g<=0.3)&(np.log10(gal['StarFormationRate']/gal['StellarMass'])<-11)
RR_x, RR_y, RR_z = np.random.uniform(0, L_box, size=(3,5*len(ra_g[s_g])))
RR_ra_g, RR_dec_g, RR_redshift_g = get_radecz_R(RR_x, RR_y, RR_z)
RR_g_sel = (RR_redshift_g>=0.1)&(RR_redshift_g<=0.3)
for LX_min in np.arange(42.7, 43.6, 0.1):
	p_2_2PCF = os.path.join(pcf_dir, 'LGAL_01z03_L0520_'+str(np.round(LX_min,2))+'_RS_1067_JK_100.fits')
	s_c = (CLUSTER_LX_soft_RF_R500c>LX_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	RRc_x, RRc_y, RRc_z = np.random.uniform(0, L_box, size=(3,20*len(ra_g[s_c])))
	RR_ra_c, RR_dec_c, RR_redshift_c = get_radecz_R(RRc_x, RRc_y, RRc_z)
	RR_c_sel = (RR_redshift_c>=0.1)&(RR_redshift_c<=0.3)
	tabulate_wprp_clustering_noW(
		ra_g[s_g], dec_g[s_g], redshift_g[s_g], RR_ra_g[RR_g_sel], RR_dec_g[RR_g_sel], RR_redshift_g[RR_g_sel],
		ra_g[s_c], dec_g[s_c], redshift_g[s_c], RR_ra_c[RR_c_sel], RR_dec_c[RR_c_sel], RR_redshift_c[RR_c_sel],
		out_file=p_2_2PCF, pimax = 100.0, N_JK=100)

s_g = (gal['StellarMass']>10**Ms_val_GAM)&(redshift_g>=0.1)&(redshift_g<=0.3)&(np.log10(gal['StarFormationRate']/gal['StellarMass'])>=-11)
RR_x, RR_y, RR_z = np.random.uniform(0, L_box, size=(3,5*len(ra_g[s_g])))
RR_ra_g, RR_dec_g, RR_redshift_g = get_radecz_R(RR_x, RR_y, RR_z)
RR_g_sel = (RR_redshift_g>=0.1)&(RR_redshift_g<=0.3)
for LX_min in np.arange(42.7, 43.6, 0.1):
	p_2_2PCF = os.path.join(pcf_dir, 'LGAL_01z03_L0520_'+str(np.round(LX_min,2))+'_BC_1067_JK_100.fits')
	s_c = (CLUSTER_LX_soft_RF_R500c>LX_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	RRc_x, RRc_y, RRc_z = np.random.uniform(0, L_box, size=(3,20*len(ra_g[s_c])))
	RR_ra_c, RR_dec_c, RR_redshift_c = get_radecz_R(RRc_x, RRc_y, RRc_z)
	RR_c_sel = (RR_redshift_c>=0.1)&(RR_redshift_c<=0.3)
	tabulate_wprp_clustering_noW(
		ra_g[s_g], dec_g[s_g], redshift_g[s_g], RR_ra_g[RR_g_sel], RR_dec_g[RR_g_sel], RR_redshift_g[RR_g_sel],
		ra_g[s_c], dec_g[s_c], redshift_g[s_c], RR_ra_c[RR_c_sel], RR_dec_c[RR_c_sel], RR_redshift_c[RR_c_sel],
		out_file=p_2_2PCF, pimax = 100.0, N_JK=100)


Ms_val_GAM = 10.23


s_g = (gal['StellarMass']>10**Ms_val_GAM)&(redshift_g>=0.1)&(redshift_g<=0.3)
RR_x, RR_y, RR_z = np.random.uniform(0, L_box, size=(3,5*len(ra_g[s_g])))
RR_ra_g, RR_dec_g, RR_redshift_g = get_radecz_R(RR_x, RR_y, RR_z)
RR_g_sel = (RR_redshift_g>=0.1)&(RR_redshift_g<=0.3)
for LX_min in np.arange(42.5, 44.5, 0.1):
	p_2_2PCF = os.path.join(pcf_dir, 'LGAL_01z03_L0520_'+str(np.round(LX_min,2))+'_Ms_1023_JK_100.fits')
	s_c = (CLUSTER_LX_soft_RF_R500c>LX_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	RRc_x, RRc_y, RRc_z = np.random.uniform(0, L_box, size=(3,20*len(ra_g[s_c])))
	RR_ra_c, RR_dec_c, RR_redshift_c = get_radecz_R(RRc_x, RRc_y, RRc_z)
	RR_c_sel = (RR_redshift_c>=0.1)&(RR_redshift_c<=0.3)
	tabulate_wprp_clustering_noW(
		ra_g[s_g], dec_g[s_g], redshift_g[s_g], RR_ra_g[RR_g_sel], RR_dec_g[RR_g_sel], RR_redshift_g[RR_g_sel],
		ra_g[s_c], dec_g[s_c], redshift_g[s_c], RR_ra_c[RR_c_sel], RR_dec_c[RR_c_sel], RR_redshift_c[RR_c_sel],
		out_file=p_2_2PCF, pimax = 100.0, N_JK=100)

s_g = (gal['StellarMass']>10**Ms_val_GAM)&(redshift_g>=0.1)&(redshift_g<=0.3)&(np.log10(gal['StarFormationRate']/gal['StellarMass'])<-11)
RR_x, RR_y, RR_z = np.random.uniform(0, L_box, size=(3,5*len(ra_g[s_g])))
RR_ra_g, RR_dec_g, RR_redshift_g = get_radecz_R(RR_x, RR_y, RR_z)
RR_g_sel = (RR_redshift_g>=0.1)&(RR_redshift_g<=0.3)
for LX_min in np.arange(42.5, 44.5, 0.1):
	p_2_2PCF = os.path.join(pcf_dir, 'LGAL_01z03_L0520_'+str(np.round(LX_min,2))+'_RS_1023_JK_100.fits')
	s_c = (CLUSTER_LX_soft_RF_R500c>LX_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	RRc_x, RRc_y, RRc_z = np.random.uniform(0, L_box, size=(3,20*len(ra_g[s_c])))
	RR_ra_c, RR_dec_c, RR_redshift_c = get_radecz_R(RRc_x, RRc_y, RRc_z)
	RR_c_sel = (RR_redshift_c>=0.1)&(RR_redshift_c<=0.3)
	tabulate_wprp_clustering_noW(
		ra_g[s_g], dec_g[s_g], redshift_g[s_g], RR_ra_g[RR_g_sel], RR_dec_g[RR_g_sel], RR_redshift_g[RR_g_sel],
		ra_g[s_c], dec_g[s_c], redshift_g[s_c], RR_ra_c[RR_c_sel], RR_dec_c[RR_c_sel], RR_redshift_c[RR_c_sel],
		out_file=p_2_2PCF, pimax = 100.0, N_JK=100)

s_g = (gal['StellarMass']>10**Ms_val_GAM)&(redshift_g>=0.1)&(redshift_g<=0.3)&(np.log10(gal['StarFormationRate']/gal['StellarMass'])>=-11)
RR_x, RR_y, RR_z = np.random.uniform(0, L_box, size=(3,5*len(ra_g[s_g])))
RR_ra_g, RR_dec_g, RR_redshift_g = get_radecz_R(RR_x, RR_y, RR_z)
RR_g_sel = (RR_redshift_g>=0.1)&(RR_redshift_g<=0.3)
for LX_min in np.arange(42.5, 44.5, 0.1):
	p_2_2PCF = os.path.join(pcf_dir, 'LGAL_01z03_L0520_'+str(np.round(LX_min,2))+'_BC_1023_JK_100.fits')
	s_c = (CLUSTER_LX_soft_RF_R500c>LX_min)&(redshift_g>=0.1)&(redshift_g<=0.3)&(gal['Type']==0)
	RRc_x, RRc_y, RRc_z = np.random.uniform(0, L_box, size=(3,20*len(ra_g[s_c])))
	RR_ra_c, RR_dec_c, RR_redshift_c = get_radecz_R(RRc_x, RRc_y, RRc_z)
	RR_c_sel = (RR_redshift_c>=0.1)&(RR_redshift_c<=0.3)
	tabulate_wprp_clustering_noW(
		ra_g[s_g], dec_g[s_g], redshift_g[s_g], RR_ra_g[RR_g_sel], RR_dec_g[RR_g_sel], RR_redshift_g[RR_g_sel],
		ra_g[s_c], dec_g[s_c], redshift_g[s_c], RR_ra_c[RR_c_sel], RR_dec_c[RR_c_sel], RR_redshift_c[RR_c_sel],
		out_file=p_2_2PCF, pimax = 100.0, N_JK=100)

