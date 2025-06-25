print('+'*100)
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import os, sys
import glob
import numpy as np
from astropy.table import Table, Column, vstack, hstack
from scipy.optimize import curve_fit
#import healpy
from scipy.interpolate import interp1d
colors = ["#67E568","#FFF000","#FFB62B","#E56124",
		  "#E53E30",
		  "#7F2353","#F911FF","#9F8CA6","#257F27","#08420D"]

fig_dir  ='../figures/'
# "C:\Users\Johan Comparat\Documents\Shared\software\st_mod_data\data"
# export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
os.environ['GIT_STMOD_DATA'] = os.path.join(os.environ['USERPROFILE'], "Documents\Shared\software\st_mod_data") # visible in this process + all children
ZuMa_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/benchmark/zu-mandelbaum-1505.02781v1')

os.environ['LSDR10'] = os.path.join(os.environ['USERPROFILE'], "Documents\Shared\data\legacysurvey\dr10")
LS10_bgs_vlim_dir = os.path.join(os.environ['LSDR10'], 'sweep\BGS_VLIM_Mstar')

d1025 = '../data/_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100'
d1075 = '../data/_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100'
d1025_BC = '../DATA/_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100'
d1075_BC = '../DATA/_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100'
d1025_RS = '../DATA/_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100'
d1075_RS = '../DATA/_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100'
dC0xG1025 = '../data/C0xG1025'
dC1xG1075 = '../data/C1xG1075'
dC1xC1 = '../data/C1xC1'
dC0xC0 = '../data/C0xC0'

ZuMa = {}
ZuMa["esd_10.2_M_10.6"]  = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_esd_10.2_M_10.6_measurements.txt" ), unpack = True)
ZuMa["esd_10.6_M_11.0"]  = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_esd_10.6_M_11.0_measurements.txt" ), unpack = True)
ZuMa["esd_11.0_M_11.2"]  = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_esd_11.0_M_11.2_measurements.txt" ), unpack = True)
ZuMa["esd_11.2_M_11.4"]  = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_esd_11.2_M_11.4_measurements.txt" ), unpack = True)
ZuMa["esd_11.4_M_12.0"]  = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_esd_11.4_M_12.0_measurements.txt" ), unpack = True)
ZuMa["esd_9.4_M_9.8"]    = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_esd_9.4_M_9.8_measurements.txt"   ), unpack = True)
ZuMa["esd_9.8_M_10.2"]   = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_esd_9.8_M_10.2_measurements.txt"  ), unpack = True)

ZuMa["wprp_10.2_M_10.6"] = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_wprp_10.2_M_10.6_measurements.txt"), unpack = True)
ZuMa["wprp_10.6_M_11.0"] = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_wprp_10.6_M_11.0_measurements.txt"), unpack = True)
ZuMa["wprp_11.0_M_11.2"] = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_wprp_11.0_M_11.2_measurements.txt"), unpack = True)
ZuMa["wprp_11.2_M_11.4"] = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_wprp_11.2_M_11.4_measurements.txt"), unpack = True)
ZuMa["wprp_11.4_M_12.0"] = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_wprp_11.4_M_12.0_measurements.txt"), unpack = True)
ZuMa["wprp_9.4_M_9.8"]   = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_wprp_9.4_M_9.8_measurements.txt"  ), unpack = True)
ZuMa["wprp_9.8_M_10.2"]  = np.loadtxt( os.path.join(ZuMa_dir, "Fig6_wprp_9.8_M_10.2_measurements.txt" ), unpack = True)

BGS = {}
BGS["ANY_10.25"] = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask.fits") )
BGS["ANY_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask.fits" ) )

BGS["RS_10.25"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.25-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
BGS["RS_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_10.75-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

BGS["BC_10.25"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.25-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
BGS["BC_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_10.75-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

CLU={}
CLU["S0_0.05_z_0.22"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.22-wprp-pimax100-bin0p1.fits"  ) )
CLU["S0_0.05_z_0.35"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.35-wprp-pimax100-bin0p1.fits"  ) )

CLU["S1_0.05_z_0.31"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.31-wprp-pimax100-bin0p1.fits"  ) )
CLU["S1_0.05_z_0.35"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.35-wprp-pimax100-bin0p1.fits"  ) )

CxG={}
CxG["S0_ANY_10.25"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S0_BC_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
CxG["S0_RS_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

CxG["S1_ANY_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S1_BC_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
CxG["S1_RS_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

z13_r, z13_wp, z13_wp_up = np.loadtxt('../data/zu-weinberg-fig13-red-datapoints-M1140-M1190.csv', unpack=True, delimiter=',')
z13_wp_ferr = z13_wp_up/z13_wp-1


#
#
# Galaxies alone
#
#

# plot wprp 1025
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS.png')
plt.figure(11, (6,5))

t_ref = BGS["ANY_10.25"]
itp_ref = interp1d(t_ref['rp_mid'], t_ref['wprp'])
t_wp = Table.read( os.path.join( d1025, 'LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02.fits' ) )
s_5_20 = (t_wp['rp_mid']>=5)&(t_wp['rp_mid']<=20)
# print('1025 b corr',t_wp['wprp'][s_5_20]/itp_ref(t_wp['rp_mid'][s_5_20]), np.mean(t_wp['wprp'][s_5_20]/itp_ref(t_wp['rp_mid'][s_5_20])))
b_corr_1025 = np.mean(t_wp['wprp'][s_5_20]/itp_ref(t_wp['rp_mid'][s_5_20])) **0.5
print('1025 b corr',b_corr_1025)
pcf_list = np.array( glob.glob( os.path.join( d1025, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
# t_wp = BGS["ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='solid', label='All, N='+str(int(t_wp['N_data'][0])), color='black')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']* t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='black', alpha=0.3)

# t_wp = BGS["RS_10.25"]
t_wp = Table.read( os.path.join( d1025_RS, 'LS10_VLIM_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02.fits' ) )
pcf_list = np.array( glob.glob( os.path.join( d1025_RS, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='darkred', alpha=0.3)

# t_wp = BGS["BC_10.25"]
t_wp = Table.read( os.path.join( d1025_BC, 'LS10_VLIM_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask-01z02.fits' ) )
pcf_list = np.array( glob.glob( os.path.join( d1025_BC, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='darkblue', alpha=0.3)

# plt.ylim((0.5, 1e4))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$ [Mpc$^2$/h$^2$]")
plt.legend(loc=3, fontsize=12, title='0.05<z<0.22,\n' +r'10.25$<\log_{10}(M*/M_\odot)<$12')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


# plot wprp 1075
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS.png')
plt.figure(11, (6,5))

itp_ref = interp1d(t_ref['rp_mid'], t_ref['wprp'])
t_wp = Table.read( os.path.join( d1075, 'LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03.fits' ) )
s_5_20 = (t_wp['rp_mid']>=5)&(t_wp['rp_mid']<=20)
# print('1075 b corr', t_wp['wprp'][s_5_20]/itp_ref(t_wp['rp_mid'][s_5_20]), np.mean(t_wp['wprp'][s_5_20]/itp_ref(t_wp['rp_mid'][s_5_20])))
b_corr_1075 = np.mean(t_wp['wprp'][s_5_20]/itp_ref(t_wp['rp_mid'][s_5_20])) **0.5
print('1075 b corr',b_corr_1075)
t_wp = Table.read( os.path.join( d1075, 'LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03.fits' ) )
pcf_list = np.array( glob.glob( os.path.join( d1075, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
# t_wp = BGS["ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='solid', label='All, N='+str(int(t_wp['N_data'][0])), color='black')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']* t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='black', alpha=0.3)

# t_wp = BGS["RS_10.25"]
t_wp = Table.read( os.path.join( d1075_RS, 'LS10_VLIM_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03.fits' ) )
pcf_list = np.array( glob.glob( os.path.join( d1075_RS, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='darkred', alpha=0.3)

# t_wp = BGS["BC_10.25"]
t_wp = Table.read( os.path.join( d1075_BC, 'LS10_VLIM_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03.fits' ) )
pcf_list = np.array( glob.glob( os.path.join( d1075_BC, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='darkblue', alpha=0.3)


plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$ [Mpc$^2$/h$^2$]")
plt.legend(loc=3, fontsize=12, title='0.05<z<0.31,\n'+r'10.75$<\log_{10}(M*/M_\odot)<$12')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


#
#
# S1 10.75
#
#


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS-S1-clusters.png')
plt.figure(3, (5., 5.))

t_wp = CLU["S1_0.05_z_0.31"]
# plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='solid', label=r'Clusters, $\log_{10}(L_X\; [\mathrm{erg/s}])>43.1$', color='orange')
pcf_list = np.array( glob.glob( os.path.join( dC1xC1, '*NSIDE_16*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
s1 = (t_wp['rp_mid']>=1)
f_err = all_wprp.std(axis=0)[s1]/all_wprp.mean(axis=0)[s1]
y_val = all_wprp.mean(axis=0)[s1]
plt.plot(t_wp['rp_mid'][s1], t_wp['rp_mid'][s1] * y_val , lw=1,  ls='dashed', label='C1', color='orange')
plt.fill_between(t_wp['rp_mid'][s1], t_wp['rp_mid'][s1] * y_val*(1-f_err), t_wp['rp_mid'][s1]*y_val*(1+f_err), color='orange', alpha=0.4)


t_wp = Table.read( os.path.join( d1075, 'LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03.fits' ) )
pcf_list = np.array( glob.glob( os.path.join( d1075, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
t_ref = BGS["ANY_10.75"]
f_err_resamp = t_ref['wprp_JK_std']/t_ref['wprp_JK_mean']
print(f_err,f_err_resamp)
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='dashed', label='G1075', color='black')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']* t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='black', alpha=0.4)


# t_wp = BGS["RS_10.25"]
t_wp = Table.read( os.path.join( d1075_RS, 'LS10_VLIM_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03.fits' ) )
pcf_list = np.array( glob.glob( os.path.join( d1075_RS, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
t_ref = BGS["RS_10.75"]
f_err_resamp = t_ref['wprp_JK_std']/t_ref['wprp_JK_mean']
print(f_err,f_err_resamp)
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='dashed', color='darkred')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='darkred', alpha=0.4)

# t_wp = BGS["BC_10.25"]
t_wp = Table.read( os.path.join( d1075_BC, 'LS10_VLIM_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03.fits' ) )
pcf_list = np.array( glob.glob( os.path.join( d1075_BC, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
t_ref = BGS["BC_10.75"]
f_err_resamp = t_ref['wprp_JK_std']/t_ref['wprp_JK_mean']
print(f_err,f_err_resamp)
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='dashed', color='darkblue')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='darkblue', alpha=0.4)

pcf_list = np.array( glob.glob( os.path.join( d1075_BC, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
t_ref = BGS["BC_10.75"]
f_err_resamp = t_ref['wprp_JK_std']/t_ref['wprp_JK_mean']
print(f_err,f_err_resamp)
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=1,  ls='dashed', color='darkblue')
plt.fill_between(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp']*(1-f_err), t_wp['rp_mid']*t_wp['wprp']*(1+f_err), color='darkblue', alpha=0.4)


# t_ref = CxG["S1_ANY_10.75"]
pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_ANY_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
# f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'],all_wprp.mean(axis=0)*t_wp['rp_mid'], lw=1,  ls='solid', label='C1 x G1075', color='black')
plt.fill_between(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid']*(1-f_err), all_wprp.mean(axis=0)*t_wp['rp_mid']*(1+f_err), color='black', alpha=0.4)

# t_wp = CxG["S1_BC_10.75"]
pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_BC_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
# f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid'], lw=1,  ls='solid', color='darkblue')
plt.fill_between(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid']*(1-f_err), all_wprp.mean(axis=0)*t_wp['rp_mid']*(1+f_err), color='darkblue', alpha=0.4)#, label='blue-cloud')

t_wp = CxG["S1_RS_10.75"]
pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_RS_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp = Table.read(el)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
# f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid'], lw=1,  ls='solid', color='darkred')
plt.fill_between(t_wp['rp_mid'], all_wprp.mean(axis=0)*t_wp['rp_mid']*(1-f_err), all_wprp.mean(axis=0)*t_wp['rp_mid']*(1+f_err), color='darkred', alpha=0.4)#, label='red-sequence')

plt.plot(z13_r, z13_r*z13_wp, label='Zu 13 (M$^*\sim11.7$))', ls='dotted', lw=1, color='green', zorder=0)
plt.fill_between(z13_r, z13_r*z13_wp*(1-z13_wp_ferr), z13_r*z13_wp*(1+z13_wp_ferr), color='green', alpha=0.4)

plt.ylim((9, 2e3))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$")#")
plt.legend(loc=4, fontsize=12,ncol=1)#, title=r'$0.1<z<0.3$, ' '$\log_{10}(L_X\; [\mathrm{erg/s}])>43.1$ x'+ "\n" +r'$10.75<\log_{10}(M*[M_\odot])<12$')#, title='LS10, r<19.5, 0.05<z<0.22')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS-S1-clusters-err.png')
plt.figure(13, (5., 5.))

t_wp = CxG["S1_ANY_10.75"]
pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_ANY_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='C0 x G1025', color='black')

t_wp = CxG["S1_BC_10.75"]
pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_BC_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='blue-cloud', color='darkblue')

t_wp = CxG["S1_RS_10.75"]
pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_RS_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='red-sequence', color='darkred')

plt.plot(z13_r, z13_wp_ferr, lw=2, label='Zu 13 (M$^*\sim11.7$)', ls='dotted', color='green', zorder=0)

plt.ylim((0.005, 0.3))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$\sigma(w_p)/w_p$")#")
#plt.legend(loc=4, fontsize=12,ncol=2)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS-S1-clusters-ratioCROSSAUTO.png')
plt.figure(13, (5., 5.))


# ANY G 1075 x ANY G 1075 REF
t_ref = Table.read( os.path.join( d1075, 'LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask-01z03.fits' ) )
pcf_list = np.array( glob.glob( os.path.join( d1075, '*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref['wprp_JK_mean'] = all_wprp.mean(axis=0)
t_ref['wprp_JK_std'] = all_wprp.std(axis=0)
f_err1 = all_wprp.std(axis=0)/all_wprp.mean(axis=0)

# ANY G 10.75 x C1
t_wp = CxG["S1_ANY_10.75"]
pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_ANY_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_wp['wprp_JK_mean'] = all_wprp.mean(axis=0)
t_wp['wprp_JK_std'] = all_wprp.std(axis=0)
f_err2 = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
f_err = (f_err1**2+f_err2**2)**0.5

plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']/t_ref['wprp_JK_mean']*(1-f_err), t_wp['wprp_JK_mean']/t_ref['wprp_JK_mean']*(1+f_err), color='black', alpha=0.4)
plt.plot(t_wp['rp_mid'], t_wp['wprp_JK_mean']/t_ref['wprp_JK_mean'], lw=2,  ls='solid', color='black')#, label='all galaxies')

t_ref = BGS["RS_10.75"]
t_wp = CxG["S1_RS_10.75"]
f_err1 = t_ref['wprp_JK_std']/t_ref['wprp_JK_mean']
pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_RS_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_err2 = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
# f_err2 = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
f_err = (f_err1**2+f_err2**2)**0.5
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp']*(1-f_err), t_wp['wprp']/t_ref['wprp']*(1+f_err), color='darkred', alpha=0.4)
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='solid', color='darkred', label='red-sequence')

t_ref = BGS["BC_10.75"]
t_wp = CxG["S1_BC_10.75"]
f_err1 = t_ref['wprp_JK_std']/t_ref['wprp_JK_mean']
pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_BC_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
f_ef_err2rr = all_wprp.std(axis=0)/all_wprp.mean(axis=0)
# f_err2 = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
f_err = (f_err1**2+f_err2**2)**0.5
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp']*(1-f_err), t_wp['wprp']/t_ref['wprp']*(1+f_err), color='darkblue', alpha=0.4)
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='solid', color='darkblue', label='blue-cloud')

b_gal = 1.43#*b_corr_1075
b_clu = 3.35
b_clu_e = 0.23 + 0.02
xx = t_wp['rp_mid']
plt.fill_between(xx, y1=np.ones_like(xx)*( ((b_clu-b_clu_e)/b_gal)), y2=np.ones_like(xx)*( ((b_clu+b_clu_e)/b_gal)), color='orange', alpha=0.2)
plt.axhline((b_clu/b_gal), color='orange', ls='dashed', label='$b_{C}/b_{G}$')

plt.ylim((0.9, 20))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w^{C1xG1075}_p/w^{G1075}_p$")
plt.legend(loc=1, fontsize=12,ncol=1)#, title='C1xG1075')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

sys.exit()
#
#
# S0 10.25
#
#

z13_r, z13_wp, z13_wp_up = np.loadtxt('../data/zu-weinberg-fig13-green-datapoints-M1120-M1140.csv', unpack=True, delimiter=',')
z13_wp_ferr = z13_wp_up/z13_wp-1

# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS-S0-clusters.png')
plt.figure(13, (5., 5.))
t_wp = BGS["ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed',  label=r'G1025', color='black')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='black', alpha=0.4)

t_wp = BGS["RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed', color='darkred')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkred', alpha=0.4)

t_wp = BGS["BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed', color='darkblue')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkblue', alpha=0.4)

#t_wp = CLU["S0_0.05_z_0.22"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='solid', label='C1', color='orange')
t_wp = CxG["S0_ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='C1xG1025', color='black')
t_wp = CxG["S0_BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', color='darkblue')#, label='blue-cloud'
t_wp = CxG["S0_RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', color='darkred')# , label='red-sequence'

plt.plot(z13_r, z13_r*z13_wp, label=r'Zu 13 (M$^*\sim11.3$)', ls='dotted', lw=2, color='green', zorder=0)

plt.ylim((1, 2e3))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$")#")
plt.legend(loc=4, fontsize=12,ncol=1, title=r'C1xG1025')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS-S0-clusters.png')
plt.figure(13, (5., 5.))
#plt.plot(z13M1120_r, z13M1120_r*z13M1120_wp, label='Zu 13, M*BCG=11.3', ls='dotted')
t_wp = BGS["ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed',  label=r'G1025 x G1025', color='black')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='black', alpha=0.4)
t_wp = BGS["RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed', color='darkred')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkred', alpha=0.4)
t_wp = BGS["BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed', color='darkblue')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkblue', alpha=0.4)

#t_wp = CLU["S0_0.05_z_0.31"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='solid', label=r'Clusters, $\log_{10}(L_X\; [\mathrm{erg/s}])>43.1$', color='orange')
t_wp = CxG["S0_ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='solid', label='C0 x G1025', color='black')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='black', alpha=0.4)
t_wp = CxG["S0_BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='solid', color='darkblue')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkblue', alpha=0.4)#, label='blue-cloud')
t_wp = CxG["S0_RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='solid', color='darkred')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkred', alpha=0.4)#, label='red-sequence')

plt.plot(z13_r, z13_r*z13_wp, label=r'Zu 13 (M$^*\sim11.3$)', ls='dotted', color='green', zorder=0)
plt.fill_between(z13_r, z13_r*z13_wp*(1-z13_wp_ferr), z13_r*z13_wp*(1+z13_wp_ferr), color='green', alpha=0.4)

plt.ylim((9, 2e3))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$")#")
plt.legend(loc=4, fontsize=12,ncol=1)#, title=r'$0.1<z<0.2$, ' '$\log_{10}(L_X\; [\mathrm{erg/s}])>42.7$ x'+ "\n" +r'$10.25<\log_{10}(M*[M_\odot])<12$')#, title='LS10, r<19.5, 0.05<z<0.22')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS-S0-clusters-err.png')
plt.figure(13, (5., 5.))
#plt.plot(z13M1120_r, z13M1120_r*z13M1120_wp, label='Zu 13, M*BCG=11.3', ls='dotted')
# t_wp = BGS["ANY_10.25"]
# f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
# plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='dashed',  label=r'G1025 x G1025', color='black')
# t_wp = BGS["RS_10.25"]
# f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
# plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='dashed', color='darkred')
# t_wp = BGS["BC_10.25"]
# f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
# plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='dashed', color='darkblue')

#t_wp = CLU["S0_0.05_z_0.31"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='solid', label=r'Clusters, $\log_{10}(L_X\; [\mathrm{erg/s}])>43.1$', color='orange')
t_wp = CxG["S0_ANY_10.25"]
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='C0 x G1025', color='black')
t_wp = CxG["S0_BC_10.25"]
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='blue-cloud', color='darkblue')
t_wp = CxG["S0_RS_10.25"]
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='red-sequence', color='darkred')

plt.plot(z13_r, z13_wp_ferr, label='Zu 13 (SDSS)', ls='dotted', color='green', zorder=0)

plt.ylim((0.01, 0.3))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$\sigma(w_p)/w_p$")#")
#plt.legend(loc=4, fontsize=12,ncol=2)#, title=r'$0.1<z<0.2$, ' '$\log_{10}(L_X\; [\mathrm{erg/s}])>42.7$ x'+ "\n" +r'$10.25<\log_{10}(M*[M_\odot])<12$')#, title='LS10, r<19.5, 0.05<z<0.22')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS-S0-clusters-ratio.png')
plt.figure(13, (5., 5.))
t_ref = BGS["ANY_10.25"]
t_wp = BGS["RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='RS, Gal x Gal', color='darkred')
t_wp = BGS["BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='BC, Gal x Gal ', color='darkblue')

t_ref = CxG["S0_ANY_10.25"]
t_wp = CxG["S0_BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S0_RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')

plt.axhline(1,color='k', ls='dotted')
plt.ylim((0, 2))
plt.xlim((0.03, 60))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)/w^{REF}_p$")
plt.legend(loc=1, fontsize=10, title='10.25<M*<12, S0 LX>42.7')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS-S0-clusters-ratioCROSSAUTO.png')
plt.figure(13, (5., 5.))
t_ref = BGS["ANY_10.25"]
t_wp = CxG["S0_ANY_10.25"]
f_err1 = t_ref['wprp_JK_std']/t_ref['wprp_JK_mean']
f_err2 = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
f_err = (f_err1[::2]**2+f_err2**2)**0.5
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2]*(1-f_err), t_wp['wprp']/t_ref['wprp'][::2]*(1+f_err), color='black', alpha=0.4)
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=2,  ls='solid', color='black')#, label='all galaxies')

t_ref = BGS["RS_10.25"]
t_wp = CxG["S0_RS_10.25"]
f_err1 = t_ref['wprp_JK_std']/t_ref['wprp_JK_mean']
f_err2 = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
f_err = (f_err1**2+f_err2**2)**0.5
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp']*(1-f_err), t_wp['wprp']/t_ref['wprp']*(1+f_err), color='darkred', alpha=0.4)
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='solid', color='darkred', label='red-sequence')

t_ref = BGS["BC_10.25"]
t_wp = CxG["S0_BC_10.25"]
f_err1 = t_ref['wprp_JK_std']/t_ref['wprp_JK_mean']
f_err2 = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
f_err = (f_err1**2+f_err2**2)**0.5
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp']*(1-f_err), t_wp['wprp']/t_ref['wprp']*(1+f_err), color='darkblue', alpha=0.4)
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='solid', color='darkblue', label='blue-cloud')

b_gal = 1.29
b_clu = 2.95
b_clu_e = 0.21
xx = t_wp['rp_mid']
plt.fill_between(xx, y1=np.ones_like(xx)*( ((b_clu-b_clu_e)/b_gal)), y2=np.ones_like(xx)*( ((b_clu+b_clu_e)/b_gal)), color='orange', alpha=0.2)
plt.axhline((b_clu/b_gal), color='orange', ls='dashed', label='$b_{C}/b_{G}$')

plt.ylim((0.9, 20))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w^{C0xG1025}_p/w^{G1025}_p$")
plt.legend(loc=1, fontsize=12,ncol=1)#, title=r'C0xG1025')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



sys.exit()

#
#
# S2 11.0
#
#

# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M110-ANY-BC-RS-S2-clusters.png')
plt.figure(13, (5., 5.))
t_wp = BGS["ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='11.0<M*<12, N='+str(int(t_wp['N_data'][0])), color='black')
t_wp = BGS["RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')

t_wp = CLU["S2_0.05_z_0.35"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='Clusters, S2 LX>43.3, N='+str(int(t_wp['N_data'][0])), color='orange')
t_wp = CxG["S2_ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='Galaxies x Clusters ', color='black')
t_wp = CxG["S2_BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S2_RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')
plt.ylim((0.1, 3e4))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)$ [Mpc/h]")
plt.legend(loc=3, fontsize=10, title='LS10, r<19.5, 0.05<z<0.35')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'wprp-obs-M110-ANY-BC-RS-S2-clusters-ratio.png')
plt.figure(13, (5., 5.))
t_ref = BGS["ANY_11.0"]
#plt.plot(t_ref['rp_mid'], t_ref['wprp'], lw=4,  ls='solid', label='11.0<M*<12, N='+str(int(t_wp['N_data'][0])), color='black')
t_wp = BGS["RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', label='RS, Gal x Gal', color='darkred')
t_wp = BGS["BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', label='BC, Gal x Gal', color='darkblue')

#t_wp = CLU["S2_0.05_z_0.35"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='Clusters, S2 LX>43.3, N='+str(int(t_wp['N_data'][0])), color='orange')
t_ref = CxG["S2_ANY_11.0"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='Galaxies x Clusters ', color='black')
t_wp = CxG["S2_BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S2_RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')
plt.axhline(1,color='k', ls='dotted')
plt.ylim((0, 2))
plt.xlim((0.03, 60))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)/w^{REF}_p$")
plt.legend(loc=1, fontsize=10, title='11.0<M*<12, S2 LX>43.3')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



p2_fig = os.path.join( fig_dir, 'wprp-obs-M110-ANY-BC-RS-S2-clusters-ratioCROSSAUTO.png')
plt.figure(13, (5., 5.))
t_ref = BGS["ANY_11.0"]
t_wp = CxG["S2_ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='Any', color='black')

t_ref = BGS["RS_11.0"]
t_wp = CxG["S2_RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='RS', color='darkred')

t_ref = BGS["BC_11.0"]
t_wp = CxG["S2_BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='BC ', color='darkblue')

plt.ylim((0.1, 100))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w^{CROSS}_p(r_p)/w^{AUTO}_p$")
plt.legend(loc=1, fontsize=10, title='11.0<M*<12, S2 LX>43.3')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

#
#
# S1 10.5
#
#

# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M105-ANY-BC-RS-S1-clusters.png')
plt.figure(13, (5., 5.))
t_wp = BGS["ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='10.5<M*<12, N='+str(int(t_wp['N_data'][0])), color='black')
t_wp = BGS["RS_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')

t_wp = CLU["S1_0.05_z_0.35"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='Clusters, S1 LX>43.1, N='+str(int(t_wp['N_data'][0])), color='orange')
t_wp = CxG["S1_ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='Galaxies x Clusters ', color='black')
t_wp = CxG["S1_BC_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S1_RS_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')

plt.ylim((0.1, 3e4))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)$ [Mpc/h]")
plt.legend(loc=3, fontsize=10, title='LS10, r<19.5, 0.05<z<0.26')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'wprp-obs-M105-ANY-BC-RS-S1-clusters-ratio.png')
plt.figure(13, (5., 5.))
t_ref = BGS["ANY_10.5"]
t_wp = BGS["RS_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', label='RS, Gal x Gal', color='darkred')
t_wp = BGS["BC_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', label='BC, Gal x Gal ', color='darkblue')

t_ref = CxG["S1_ANY_10.5"]
t_wp = CxG["S1_BC_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S1_RS_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')

plt.axhline(1,color='k', ls='dotted')
plt.ylim((0, 2))
plt.xlim((0.03, 60))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)/w^{REF}_p$")
plt.legend(loc=1, fontsize=10, title='10.5<M*<12, S1 LX>43.1')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'wprp-obs-M105-ANY-BC-RS-S1-clusters-ratioCROSSAUTO.png')
plt.figure(13, (5., 5.))
t_ref = BGS["ANY_10.5"]
t_wp = CxG["S1_ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='Any', color='black')

t_ref = BGS["RS_10.5"]
t_wp = CxG["S1_RS_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='RS', color='darkred')

t_ref = BGS["BC_10.5"]
t_wp = CxG["S1_BC_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='BC ', color='darkblue')

plt.ylim((0.1, 100))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w^{CROSS}_p(r_p)/w^{AUTO}_p$")
plt.legend(loc=1, fontsize=10, title='10.5<M*<12, S1 LX>43.1')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



#
#
# S0 10.0
#
#



# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M100-ANY-BC-RS-S0-clusters.png')
plt.figure(13, (5., 5.))
t_wp = BGS["ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='10.0<M*<12, N='+str(int(t_wp['N_data'][0])), color='black')
t_wp = BGS["RS_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')

t_wp = CLU["S0_0.05_z_0.35"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='Clusters, S0 LX>42.7, N='+str(int(t_wp['N_data'][0])), color='orange')
t_wp = CxG["S0_ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='Galaxies x Clusters ', color='black')
t_wp = CxG["S0_BC_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S0_RS_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')

plt.ylim((0.1, 3e4))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)$ [Mpc/h]")
plt.legend(loc=3, fontsize=10, title='LS10, r<19.5, 0.05<z<0.18')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'wprp-obs-M100-ANY-BC-RS-S0-clusters-ratio.png')
plt.figure(13, (5., 5.))
t_ref = BGS["ANY_10.0"]
t_wp = BGS["RS_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', label='RS, Gal x Gal', color='darkred')
t_wp = BGS["BC_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', label='BC, Gal x Gal ', color='darkblue')

t_ref = CxG["S0_ANY_10.0"]
t_wp = CxG["S0_BC_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S0_RS_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')

plt.axhline(1,color='k', ls='dotted')
plt.ylim((0, 2))
plt.xlim((0.03, 60))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)/w^{REF}_p$")
plt.legend(loc=1, fontsize=10, title='10.0<M*<12, S0 LX>42.7')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'wprp-obs-M100-ANY-BC-RS-S0-clusters-ratioCROSSAUTO.png')
plt.figure(13, (5., 5.))
t_ref = BGS["ANY_10.0"]
t_wp = CxG["S0_ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='Any', color='black')

t_ref = BGS["RS_10.0"]
t_wp = CxG["S0_RS_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='RS', color='darkred')

t_ref = BGS["BC_10.0"]
t_wp = CxG["S0_BC_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='BC ', color='darkblue')

plt.ylim((0.1, 100))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w^{CROSS}_p(r_p)/w^{AUTO}_p$")
plt.legend(loc=1, fontsize=10, title='10.0<M*<12, S0 LX>42.7')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


sys.exit()

p2_fig = os.path.join( fig_dir, 'all_wp_LS10_pimax100.png')
plt.figure(1, (6,5))
for kk, cc in zip(list(BGS.keys()), colors):
    print(kk)
    t_wp = BGS[kk]
    plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label=kk.split('_')[1], color=cc)

plt.text(0.015, 340, r"LS10 galaxies" )
plt.ylim((-100,400))
plt.xlim((1e-2,40))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc='lower right', ncol=3)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

p2_fig = os.path.join( fig_dir, 'all_wp_UCHUU_pimax100.png')
plt.figure(1, (6,5))
for kk, cc in zip(list(UCHUU.keys()), colors):
    print(kk)
    t_wp = UCHUU[kk]
    plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label=kk.split('_')[1], color=cc)

plt.text(0.02, 340, r"Uchuu mock galaxies" )
plt.ylim((-100,400))
plt.xlim((1e-2,40))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc='lower right', ncol=3)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

sys.exit()

# make a table containing all correlation functions
ALL_WP = {}
ALL_WP['rp_min'] = np.round(np.log10 (t_wp['rp_min']),2)
ALL_WP['rp_max'] = np.round(np.log10 (t_wp['rp_max']),2)
ALL_WP['rp_mid'] = np.round(t_wp['rp_mid'],3)
for kk, cc in zip(list(BGS.keys()), colors):
    print(kk)
    ALL_WP[kk] = np.round(BGS[kk]['wprp'],1)

t_all = Table(ALL_WP)
t_all.write(os.path.join( os.environ['GIT_DR10W'], 'latex', 'all-wprp.latex'), format='latex', overwrite = True)
print(os.path.join( os.environ['GIT_DR10W'], 'latex', 'all-wprp.latex'))



# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-lit-M115.png')
plt.figure(12, (5,5))
t_wp = BGS_noSys["ANY_11.5"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='dashed', label='LS10, 11.5<M*<12')
t_wp = BGS["ANY_11.5"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label='& pixel mask')
t_wp = BGS_ebv["ANY_11.5"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='dashed', label='& EBV<0.07')

t_wp = UCHUU['ALL_11.5']
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2, ls='solid', label='UCHUU mock')

plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc=0, fontsize=12)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-lit-M1125.png')
plt.figure(12, (5,5))

t_wp = BGS_noSys["ANY_11.25"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='dashed', label='LS10, 11.25<M*<12')
t_wp = BGS["ANY_11.25"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label='& pixel mask')
t_wp = BGS_ebv["ANY_11.25"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='dashed', label='& EBV<0.07')

t_wp = UCHUU['ALL_11.25']
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2, ls='solid', label='UCHUU mock')

plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc=0, fontsize=12)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-lit-M110.png')
plt.figure(12, (5,5))
t_wp = BGS_noSys["ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='dashed', label='LS10, 11.0<M*<12')
t_wp = BGS["ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label='& pixel mask')
t_wp = BGS_ebv["ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='dashed', label='& EBV<0.07')

t_wp = UCHUU['ALL_11.0']
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2, ls='solid', label='UCHUU mock')

plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc=0, fontsize=12)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



#plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-lit-M100.png')
plt.figure(12, (5,5))
t_wp = BGS_noSys["ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='dashed', label='LS10, 10.0<M*<12')
t_wp = BGS["ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label='& pixel mask')
t_wp = BGS_ebv["ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='dashed', label='& EBV<0.07')
lit = ZuMa["wprp_9.8_M_10.2"]
if len(lit)==3:
    plt.fill_between(lit[0], lit[0]*lit[1], lit[0]*lit[2], label='Zu 16, 9.8<M*<10.2', alpha=0.5)
if len(lit)==2:
    plt.plot(lit[0], lit[0]*lit[1], label='Zu 16, 9.8<M*<10.2')
t_wp = UCHUU['ALL_10.0']
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2, ls='solid', label='UCHUU mock')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc=0, fontsize=12)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-lit-M1075.png')
plt.figure(12, (5,5))
t_wp = BGS_noSys["ANY_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='dashed', label='LS10, 10.75<M*<12')
t_wp = BGS["ANY_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label='& pixel mask')
t_wp = BGS_ebv["ANY_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='dashed', label='& EBV<0.07')
lit = ZuMa["wprp_10.6_M_11.0"]
if len(lit)==3:
    plt.fill_between(lit[0], lit[0]*lit[1], lit[0]*lit[2], label='Zu 16, 10.6<M*<11.0', alpha=0.5)
if len(lit)==2:
    plt.plot(lit[0], lit[0]*lit[1], label='Zu 16, 10.6<M*<11.0')
t_wp = UCHUU['ALL_10.75']
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2, ls='solid', label='UCHUU mock')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc=0, fontsize=12)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-lit-M105.png')
plt.figure(12, (5,5))
t_wp = BGS_noSys["ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='dashed', label='LS10, 10.5<M*<12')
t_wp = BGS["ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label='& pixel mask')
t_wp = BGS_ebv["ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='dashed', label='& EBV<0.07')
lit = ZuMa["wprp_10.2_M_10.6"]
if len(lit)==3:
    plt.fill_between(lit[0], lit[0]*lit[1], lit[0]*lit[2], label='Zu 16, 10.2<M*<10.6', alpha=0.5)
if len(lit)==2:
    plt.plot(lit[0], lit[0]*lit[1], label='Zu 16, 10.2<M*<10.6')
t_wp = UCHUU['ALL_10.5']
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2, ls='solid', label='UCHUU mock')
#t_wp = UCHUU_zerr_002_Merr_010['ALL_10.5']
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc=0, fontsize=12)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



## plot wprp
#p2_fig = os.path.join( fig_dir, 'wprp-obs-lit-M1025.png')
#plt.figure(12, (5,5))
#t_wp = BGS_noSys["ANY_10.25"]
#plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='dashed', label='LS10, 10.25<M*<12')
#t_wp = BGS["ANY_10.25"]
#plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label='& pixel mask')
#t_wp = BGS_ebv["ANY_10.25"]
#plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='dashed', label='& EBV<0.07')
#lit = ZuMa["wprp_10.2_M_10.6"]
#if len(lit)==3:
    #plt.fill_between(lit[0], lit[0]*lit[1], lit[0]*lit[2], label='Zu 16, 10.2<M*<10.6', alpha=0.5)
#if len(lit)==2:
    #plt.plot(lit[0], lit[0]*lit[1], label='Zu 16, 10.2<M*<10.6')
#t_wp = UCHUU['ALL_10.25']
#plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2, ls='solid', label='UCHUU mock')
#plt.xscale('log')
##plt.yscale('log')
#plt.xlabel(r"$r_p$ [Mpc/h]")
#plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
#plt.legend(loc=0, fontsize=12)
#plt.tight_layout()
#plt.savefig(p2_fig)
#plt.clf()
#print(p2_fig)



#plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-lit-M100.png')
plt.figure(12, (5,5))
t_wp = BGS_noSys["ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='dashed', label='LS10, 10.0<M*<12')
t_wp = BGS["ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label='& pixel mask')
t_wp = BGS_ebv["ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='dashed', label='& EBV<0.07')
lit = ZuMa["wprp_9.8_M_10.2"]
if len(lit)==3:
    plt.fill_between(lit[0], lit[0]*lit[1], lit[0]*lit[2], label='Zu 16, 9.8<M*<10.2', alpha=0.5)
if len(lit)==2:
    plt.plot(lit[0], lit[0]*lit[1], label='Zu 16, 9.8<M*<10.2')
t_wp = UCHUU['ALL_10.0']
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2, ls='solid', label='UCHUU mock')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc=0, fontsize=12)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

sys.exit()

# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-lit-M95.png')
plt.figure(12, (5,5))
t_wp = BGS_noSys["ANY_9.5"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='dashed', label='LS10, 9.5<M*<12')
t_wp = BGS["ANY_9.5"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='solid', label='& pixel mask')
t_wp = BGS_ebv["ANY_9.5"]
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=3,  ls='dashed', label='& EBV<0.07')
lit = ZuMa["wprp_9.4_M_9.8"]
if len(lit)==3:
    plt.fill_between(lit[0], lit[0]*lit[1], lit[0]*lit[2], label='Zu 16, 9.4<M*<9.8', alpha=0.5)
if len(lit)==2:
    plt.plot(lit[0], lit[0]*lit[1], label='Zu 16, 9.4<M*<9.8')
t_wp = UCHUU['ALL_9.5']
plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2, ls='solid', label='UCHUU mock')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p \times w_p(r_p)$ $\pi_{max}=100$")
plt.legend(loc=0, fontsize=12)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



