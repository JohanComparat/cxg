from astropy.table import Table, Column, hstack, vstack
import time
t0 = time.time()

import numpy as n
import numpy as np
import os, sys, glob
import astropy.io.fits as fits

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

import astropy.units as u
import astropy.constants as cc
from astropy.cosmology import FlatLambdaCDM
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT

from scipy.interpolate import interp1d
z_array = n.arange(0.01, 0.4, 0.0001)
dl2_4pi_itp = interp1d (z_array, (4*n.pi*cosmo.luminosity_distance(z_array).to(u.cm)**2).value )
root_directory = '/data56s/comparat/erosim/data_s4_c030/'
root_directory = os.environ['DATA_S4'] #'/data56s/comparat/erosim/data_s4_c030/'

prefix = 'LS10VLIM'
#cat_name = "ALL-ANY-LS_DR10_8R20_10.5_m_12_0.05_zPHOT_0.17-Nmax-R0-DATA_WITH_PROFILES"
cat_name = sys.argv[1]
E_min = 500
E_max = 2000
#E_min = int(sys.argv[2])
#E_max = int(sys.argv[3])
mergedOUT_dir = os.path.join(root_directory, 'mergedCubes_'+str(E_min)+'_'+str(E_max), prefix )
#mergedCube_dir = os.path.join( root_directory, 'mergedCubes', prefix )
#p_2_incat = os.path.join(mergedCube_dir, cat_name + '_wEVT.fits')
#gal_info = Table.read(p_2_incat)

GAL_i = Table.read(os.path.join(os.environ['LSDR10'], 'sweep/MergeALL_BGSlike_LPH.fits'))

m0 = float(cat_name.split('_')[3])
m1 = float(cat_name.split('_')[5])
z0 = float(cat_name.split('_')[6])
z1 = float(cat_name.split('_')[8])
GAL = GAL_i[(GAL_i['LPH_MASS_BEST']>m0)&(GAL_i['LPH_MASS_BEST']<m1)&(GAL_i['Z_PHOT_MEAN']>z0)&(GAL_i['Z_PHOT_MEAN']<z1)]

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
print('is_RS', len(GAL['UID'][GAL['is_RS']]))
print('is_BC', len(GAL['UID'][GAL['is_BC']]))
print('is_GV', len(GAL['UID'][GAL['is_GV']]))
#dir_2_cat = os.path.join(os.environ['LSDR10'], 'sweep/BGS_VLIM_Mstar')
#p2_input_Cat = os.path.join(dir_2_cat, cat_name + '_DATA.fits')
#GAL = Table.read(p2_input_Cat)

p2_cat = os.path.join(mergedOUT_dir, cat_name +'.fits')
print(p2_cat)

KC = ['UID', 'BEST_Z', 'LPH_MASS_BEST', 'DD_Nw', 'R0_Nw', 'RR_Nw', 'PS_Nw', 'DD_N', 'R0_N', 'RR_N', 'PS_N']

t1 = []
if os.path.isfile(p2_cat)==False:
    p2_cats = np.array(glob.glob(os.path.join(mergedOUT_dir, cat_name+'*-SPLIT_N_*.fits')))
    if len(p2_cats)>=1:
        print(p2_cats[0])
        t0 = Table()
        hdu1 = fits.open(p2_cats[0])[1]
        for el in KC:
            #print(el)
            t0[el] = hdu1.data[el]
        t1.append( t0 )
        print(len(t0))
        for p2_cats_i in p2_cats[1:]:
            print(p2_cats_i)
            t0 = Table()
            hdu1 = fits.open(p2_cats_i)[1]
            for el in KC:
                #print(el)
                t0[el] = hdu1.data[el]
            print(len(t0))
            t1.append( t0 )
        print('LOADED')
        full_cat_any = vstack(( t1 ))
        t0=0
        t1=0


full_cat_any['rand'] = n.random.random(len(full_cat_any))
t0 = time.time()

def stack_me(SF_flag = 'BC'):
    full_cat = full_cat_any[np.isin(full_cat_any['UID'], GAL['UID'][GAL['is_'+SF_flag]])]
    N_gal = len(full_cat)
    print(N_gal)
    p2_out = os.path.join(mergedOUT_dir, cat_name.replace('ANY', SF_flag) + '_STACKEDprofiles.fits' )
    print(p2_out)

    t_out = Table()

    redshift_array = full_cat['BEST_Z']
    #M_halo = n.log10(full_cat['M_halo'])
    MS = full_cat['LPH_MASS_BEST']
    #is_QU  = full_cat['color_flag']==1
    #is_SF  = full_cat['color_flag']==0
    if SF_flag=='ANY' :
        sSFR = n.ones_like(MS)* (-11)
    elif SF_flag=='rmRS' or SF_flag=='cmRS' or SF_flag=='RS' :
        sSFR = n.ones_like(MS)* (-12.5)
    elif SF_flag=='rmBC' or SF_flag=='cmBC' or SF_flag=='BC':
        sSFR = n.ones_like(MS)* (-9.5)
    SFR = sSFR + MS

    l_bins = n.hstack(( 0., 10**n.arange(1.2, n.log10(3000.), 0.1), 3000. ))
    x_l_bins = ( l_bins[1:] + l_bins[:-1] ) / 2.
    bins = l_bins # n.hstack(( 0., 10**n.arange(1.2, n.log10(3000.), 0.2), 3000. ))
    t_out['x_lo'] = bins[:-1]
    t_out['x_up'] = bins[1:]
    area_shell = ( bins[1:]**2 - bins[:-1]**2 ) * n.pi
    sys_err  = 0.02
    norm_sys = 0.02

    print('JK 90 percent of sample estimation')
    DD_jk = []
    R0_jk = []
    RR_jk = []
    PS_jk = []
    for jj in n.arange(1000):
        s_jk = ( n.random.random(N_gal) < 0.9 )
        DD_jk .append( full_cat['DD_Nw'][s_jk].sum(axis=0) / (0.9*N_gal) / area_shell )
        R0_jk .append( full_cat['R0_Nw'][s_jk].sum(axis=0) / (0.9*N_gal) / area_shell )
        RR_jk .append( full_cat['RR_Nw'][s_jk].sum(axis=0) / (0.9*N_gal) / area_shell )
        PS_jk .append( full_cat['PS_Nw'][s_jk].sum(axis=0) / (0.9*N_gal) / area_shell )
        print(jj, 'dt-',time.time()-t0)
    DD_jk_err = n.std(n.array(DD_jk), axis=0) / n.mean(n.array(DD_jk), axis=0)
    R0_jk_err = n.std(n.array(R0_jk), axis=0) / n.mean(n.array(R0_jk), axis=0)
    RR_jk_err = n.std(n.array(RR_jk), axis=0) / n.mean(n.array(RR_jk), axis=0)
    PS_jk_err = n.std(n.array(PS_jk), axis=0) / n.mean(n.array(PS_jk), axis=0)

    t_out['dd_profile'    ] = full_cat['DD_Nw'].sum(axis=0) / N_gal / area_shell
    t_out['dd_profile_err'] = DD_jk_err * t_out['dd_profile'    ]
    t_out['dd_profile_N' ] = full_cat['DD_N'].sum(axis=0)
    t_out['dd_xb'         ] = x_l_bins
    t_out['bg'         ]    = t_out['dd_profile'    ][-1]
    t_out['dd_profile_1k_JK'] = np.transpose(DD_jk)

    t_out['r0_profile'    ] = full_cat['R0_Nw'].sum(axis=0) / N_gal / area_shell
    t_out['r0_profile_err'] = R0_jk_err * t_out['r0_profile'    ]
    t_out['r0_profile_N' ] = full_cat['R0_N'].sum(axis=0)
    t_out['r0_xb'         ] = x_l_bins
    t_out['r0_profile_1k_JK'] = np.transpose(R0_jk)

    t_out['rr_profile'    ] = full_cat['RR_Nw'].sum(axis=0) / N_gal / area_shell
    t_out['rr_profile_err'] = RR_jk_err * t_out['rr_profile'    ]
    t_out['rr_profile_N' ] = full_cat['RR_N'].sum(axis=0)
    t_out['rr_xb'         ] = x_l_bins
    t_out['rr_profile_1k_JK'] = np.transpose(RR_jk)

    t_out['ps_profile'    ] = full_cat['PS_Nw'].sum(axis=0) / N_gal / area_shell
    t_out['ps_profile_err'] = PS_jk_err * t_out['ps_profile'    ]
    t_out['ps_profile_N' ] = full_cat['PS_N'].sum(axis=0)
    t_out['ps_xb'         ] = x_l_bins
    t_out['ps_profile_1k_JK'] = np.transpose(PS_jk)

    # slope 2.0,
    frac_op20_slope2 = n.log(1.5/0.7) / n.log(2/0.5)
    # slope 1.8, Basu-Zych et al. 2020 paper
    frac_op20 = ( (1.5**0.2 - 0.7**0.2)/0.2  ) / ( (2**0.2 - 0.5**0.2)/0.2 )

    # slope 2.0,
    frac_2_10 = n.log(1.5/0.7) / n.log(2/0.5)
    # slope 1.8, Basu-Zych et al. 2020 paper
    frac_2_10 = ( (2.0**0.2 - 0.5**0.2)/0.2  ) / ( (10**0.2 - 2**0.2)/0.2 ) # 0.63
    # adding a nH 0f 4e20 (mean of the galaxies), it decreases to :
    frac_2_10 = 0.56

    # Hard X-ray emission, after Aird et al. 2017
    def galaxy_lx_2_10(redshift, log10_mass, log10_sfr):
        log10_mass_dep = 28.81 + 3.9*n.log10(1 + redshift) + log10_mass
        log10_sfr_dep = 39.5 + 0.67*n.log10( 1 + redshift) + 0.86 * log10_sfr
        return log10_mass_dep , log10_sfr_dep

    def galaxy_lx_2_10_up(redshift, log10_mass, log10_sfr):
        log10_mass_dep = 28.81 + 0.08 + ( 3.9 + 0.36 ) * n.log10(1 + redshift) + log10_mass
        log10_sfr_dep = 39.5 + 0.06 + ( 0.67 + 0.31 ) *n.log10( 1 + redshift) + ( 0.86 + 0.05 ) * log10_sfr
        return log10_mass_dep , log10_sfr_dep

    def galaxy_lx_2_10_lo(redshift, log10_mass, log10_sfr):
        log10_mass_dep = 28.81 - 0.08 + ( 3.9 - 0.36 ) * n.log10(1 + redshift) + log10_mass
        log10_sfr_dep = 39.5 - 0.06 + ( 0.67 - 0.31 ) *n.log10( 1 + redshift) + ( 0.86 - 0.05 ) * log10_sfr
        return log10_mass_dep , log10_sfr_dep

    # Galaxy hard X-ray luminosity
    log10_LX_2_10_gal_mass   , log10_LX_2_10_gal_sfr    = galaxy_lx_2_10   (redshift_array, MS, SFR )
    log10_LX_2_10_gal_mass_up, log10_LX_2_10_gal_sfr_up = galaxy_lx_2_10_up(redshift_array, MS, SFR )
    log10_LX_2_10_gal_mass_lo, log10_LX_2_10_gal_sfr_lo = galaxy_lx_2_10_lo(redshift_array, MS, SFR )
    LX_tot    = n.log10(10**(log10_LX_2_10_gal_mass-35) + 10**(log10_LX_2_10_gal_sfr-35)) + 35 + n.log10(frac_2_10) # log10 erg/s
    LX_tot_up = n.log10(10**(log10_LX_2_10_gal_mass_up-35) + 10**(log10_LX_2_10_gal_sfr_up-35)) + 35 + n.log10(frac_2_10) #
    LX_tot_lo = n.log10(10**(log10_LX_2_10_gal_mass_lo-35) + 10**(log10_LX_2_10_gal_sfr_lo-35)) + 35 + n.log10(frac_2_10) #

    med_LX_all      = n.mean(LX_tot)
    med_LX_all_up   = n.mean(LX_tot_up)
    med_LX_all_lo   = n.mean(LX_tot_lo)

    t_out['XRB_Ai17_med_LX_all'] = med_LX_all
    t_out['XRB_Ai17_med_LX_all_up'] = med_LX_all_up
    t_out['XRB_Ai17_med_LX_all_lo'] = med_LX_all_lo

    def meanSM(Mh, z): return n.log10(Mh * 2. * (0.0351 - 0.0247 * z / (1. + z)) / (
                (Mh / (10 ** (11.79 + 1.5 * z / (1. + z)))) ** (- 0.9 + 0.5 * z / (1. + z)) + (
                    Mh / (10 ** (11.79 + 1.5 * z / (1. + z)))) ** (0.67 + 0.2 * z / (1. + z))))

    all_mh = n.arange(10, 18, 0.1)
    med_Z = n.median(redshift_array) * n.ones_like(all_mh)
    SM_val = meanSM(10**all_mh, med_Z)
    itp = interp1d(SM_val, all_mh)
    mean_halo_mass = itp(n.median(MS))

    t_out['mean_halo_mass'] = mean_halo_mass
    t_out['N_gal'] = N_gal
    t_out['MS_mean'] = n.mean(MS)
    t_out['MS_std'] = n.std(MS)
    t_out['SFR_mean'] = n.mean(SFR)
    t_out['SFR_std'] = n.std(SFR)
    t_out['redshift_mean'] = n.mean(redshift_array)
    t_out['redshift_std']  = n.std(redshift_array)

    t_out.write(p2_out, overwrite = True)
    print(p2_out, 'written')

    mat = np.transpose(DD_jk)
    cv = np.cov(mat)
    cc = np.corrcoef(mat)
    p2_out_cv = os.path.join(mergedOUT_dir, cat_name.replace('ANY', SF_flag)  + '_STACKEDprofiles_CovarianceMatrix_DD.ascii' )
    np.savetxt(p2_out_cv, cv )
    p2_out_cc = os.path.join(mergedOUT_dir, cat_name.replace('ANY', SF_flag)  + '_STACKEDprofiles_CorrelationCoefficient_DD.ascii' )
    np.savetxt(p2_out_cc, cc )

    mat = np.transpose(R0_jk)
    cv = np.cov(mat)
    cc = np.corrcoef(mat)
    p2_out_cv = os.path.join(mergedOUT_dir, cat_name.replace('ANY', SF_flag)  + '_STACKEDprofiles_CovarianceMatrix_R0.ascii' )
    np.savetxt(p2_out_cv, cv )
    p2_out_cc = os.path.join(mergedOUT_dir, cat_name.replace('ANY', SF_flag)  + '_STACKEDprofiles_CorrelationCoefficient_R0.ascii' )
    np.savetxt(p2_out_cc, cc )

    mat = np.transpose(RR_jk)
    cv = np.cov(mat)
    cc = np.corrcoef(mat)
    p2_out_cv = os.path.join(mergedOUT_dir, cat_name.replace('ANY', SF_flag)  + '_STACKEDprofiles_CovarianceMatrix_RR.ascii' )
    np.savetxt(p2_out_cv, cv )
    p2_out_cc = os.path.join(mergedOUT_dir, cat_name.replace('ANY', SF_flag)  + '_STACKEDprofiles_CorrelationCoefficient_RR.ascii' )
    np.savetxt(p2_out_cc, cc )

    mat = np.transpose(PS_jk)
    cv = np.cov(mat)
    cc = np.corrcoef(mat)
    p2_out_cv = os.path.join(mergedOUT_dir, cat_name.replace('ANY', SF_flag)  + '_STACKEDprofiles_CovarianceMatrix_PS.ascii' )
    np.savetxt(p2_out_cv, cv )
    p2_out_cc = os.path.join(mergedOUT_dir, cat_name.replace('ANY', SF_flag)  + '_STACKEDprofiles_CorrelationCoefficient_PS.ascii' )
    np.savetxt(p2_out_cc, cc )

stack_me(SF_flag = 'BC')
stack_me(SF_flag = 'RS')
