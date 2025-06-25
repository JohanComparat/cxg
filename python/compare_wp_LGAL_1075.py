"""
python compare_wp.py 0 40
python compare_wp.py 1 40
python compare_wp.py 2 40
python compare_wp.py 3 40
python compare_wp.py 4 40
python compare_wp.py 5 40
python compare_wp.py 6 40
python compare_wp.py 7 40

python compare_wp.py 0 60
python compare_wp.py 1 60
python compare_wp.py 2 60
python compare_wp.py 3 60
python compare_wp.py 4 60
python compare_wp.py 5 60
python compare_wp.py 6 60
python compare_wp.py 7 60

python compare_wp.py 0 100
python compare_wp.py 1 100
python compare_wp.py 2 100
python compare_wp.py 3 100
python compare_wp.py 4 100
python compare_wp.py 5 100
python compare_wp.py 6 100
python compare_wp.py 7 100

"""
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
import healpy
from scipy.interpolate import interp1d
colors = ["#67E568","#FFF000","#FFB62B","#E56124",
		  "#E53E30",
		  "#7F2353","#F911FF","#9F8CA6","#257F27","#08420D"]

fig_dir  ='../figures/uchuu'

os.environ['GIT_STMOD_DATA'] = os.path.join(os.environ['USERPROFILE'], "Documents\Shared\software\st_mod_data") # visible in this process + all children
ZuMa_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/benchmark/zu-mandelbaum-1505.02781v1')

os.environ['LSDR10'] = os.path.join(os.environ['USERPROFILE'], "Documents\Shared\data\legacysurvey\dr10")
LS10_bgs_vlim_dir = os.path.join(os.environ['LSDR10'], 'sweep\BGS_VLIM_Mstar')

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
#BGS["ANY_9.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_9.0_Mstar_12.0_0.05_z_0.08_N_0523486-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
#BGS["ANY_9.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_9.5_Mstar_12.0_0.05_z_0.12_N_1432502-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
BGS["ANY_10.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["ANY_10.25"] = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p05-HpxMask.fits") )
BGS["ANY_10.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["ANY_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["ANY_11.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["ANY_11.25"] = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_11.25_Mstar_12.0_0.05_z_0.35_N_0541855-wprp-pimax100-bin0p05-HpxMask.fits") )
BGS["ANY_11.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p05-HpxMask.fits" ) )

BGS["RS_10.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_10.0-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["RS_10.25"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.25-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["RS_10.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-MstarMin_10.5-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["RS_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_10.75-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["RS_11.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["RS_11.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p05-HpxMask.fits" ) )

BGS["BC_10.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_10.0-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["BC_10.25"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.25-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["BC_10.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-MstarMin_10.5-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["BC_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_10.75-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["BC_11.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["BC_11.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p05-HpxMask.fits" ) )




uchuu_dir_z0p14 = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/WPRP', 'z0p14')
uchuu_dir_z0p19 = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/WPRP', 'z0p19')
uchuu_dir_z0p25 = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/WPRP', 'z0p25')
U0p14 = {}
U0p14["BC_10.25_LX42.7"] = Table.read(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_BC_HaloLX_gt_42.7_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
U0p14["RS_10.25_LX42.7"] = Table.read(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_RS_HaloLX_gt_42.7_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
U0p14["ANY_10.25_LX42.7"] = Table.read(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_HaloLX_gt_42.7_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"    ))
U0p14["ANY_10.25"] = Table.read(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                   ))
U0p14["BC_10.25"] = Table.read(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_BC_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))
U0p14["RS_10.25"] = Table.read(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_RS_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))
U0p19 = {}
U0p19["BC_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_BC_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
U0p19["RS_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_RS_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
U0p19["ANY_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"    ))
U0p19["ANY_10.75"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                   ))
U0p19["RS_10.75"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_RS_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))
U0p19["BC_10.75"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_BC_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))
U0p25 = {}
U0p25["BC_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_BC_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
U0p25["RS_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_RS_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
U0p25["ANY_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"    ))
U0p25["ANY_10.75"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                   ))
U0p25["RS_10.75"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_RS_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))
U0p25["BC_10.75"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_BC_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))


## plot wprp
#p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS-uchuu.png')
#plt.figure(11, (6,5))
#t_wp = BGS["ANY_10.75"]
#plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=4,  ls='solid', label='All galaxies', color='grey')

#t_wp = U0p19["ANY_10.75"]
#plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2,  ls='dashed', color='k', label='Uchuu+UM z=0.19')

##plt.ylim((0.5, 1e4))
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel(r"$r_p$ [Mpc/h]")
#plt.ylabel(r"$r_p\times w_p(r_p)$ [Mpc/h], $\pi_{max}=100$ Mpc$/h$")
#plt.legend(loc=4, fontsize=12, title='0.05<z<0.31,\n'+r'10.75$<\log_{10}(M*/M_\odot)<$12')
#plt.tight_layout()
#plt.savefig(p2_fig)
#plt.clf()
#print(p2_fig)

CLU={}
CLU["S0_0.05_z_0.18"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.18-wprp-pimax100-bin0p1.fits"  ) )
CLU["S0_0.05_z_0.22"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.22-wprp-pimax100-bin0p1.fits"  ) )
CLU["S0_0.05_z_0.26"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.26-wprp-pimax100-bin0p1.fits"  ) )
CLU["S0_0.05_z_0.35"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S0_0.05_z_0.35-wprp-pimax100-bin0p1.fits"  ) )
CLU["S1_0.05_z_0.18"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.18-wprp-pimax100-bin0p1.fits"  ) )
CLU["S1_0.05_z_0.26"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.26-wprp-pimax100-bin0p1.fits"  ) )
CLU["S1_0.05_z_0.31"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.31-wprp-pimax100-bin0p1.fits"  ) )
CLU["S1_0.05_z_0.35"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S1_0.05_z_0.35-wprp-pimax100-bin0p1.fits"  ) )
CLU["S2_0.05_z_0.18"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S2_0.05_z_0.18-wprp-pimax100-bin0p1.fits"  ) )
CLU["S2_0.05_z_0.26"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S2_0.05_z_0.26-wprp-pimax100-bin0p1.fits"  ) )
CLU["S2_0.05_z_0.35"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S2_0.05_z_0.35-wprp-pimax100-bin0p1.fits"  ) )
CLU["S3_0.05_z_0.18"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S3_0.05_z_0.18-wprp-pimax100-bin0p1.fits"  ) )
CLU["S3_0.05_z_0.26"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S3_0.05_z_0.26-wprp-pimax100-bin0p1.fits"  ) )
CLU["S3_0.05_z_0.35"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S3_0.05_z_0.35-wprp-pimax100-bin0p1.fits"  ) )
CLU["S4_0.05_z_0.18"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S4_0.05_z_0.18-wprp-pimax100-bin0p1.fits"  ) )
CLU["S4_0.05_z_0.26"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S4_0.05_z_0.26-wprp-pimax100-bin0p1.fits"  ) )
CLU["S4_0.05_z_0.35"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "AUTOCORR_eRASS1_VLIM_CLUSTERS_S4_0.05_z_0.35-wprp-pimax100-bin0p1.fits"  ) )
CxG={}
CxG["S0_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S0_ANY_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S0_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S0_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S0_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask.fits"  ) )

CxG["S1_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S1_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S1_ANY_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S1_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S1_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask.fits"  ) )

CxG["S2_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S2_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S2_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S2_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask.fits"  ) )

CxG["S3_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S3_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S3_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S3_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask.fits"  ) )

CxG["S4_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S4_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S4_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S4_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask.fits"  ) )

CxG["S2_BC_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_BC_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits" ) )
CxG["S2_RS_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_RS_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits" ) )

CxG["S1_BC_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_BC_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits" ) )
CxG["S1_RS_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_RS_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits" ) )

CxG["S0_BC_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_BC_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits" ) )
CxG["S0_RS_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_RS_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits" ) )


CxG["S1_BC_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
CxG["S1_RS_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

CxG["S0_BC_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
CxG["S0_RS_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )


#uchuu_dir_z0p14 = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/WPRP', 'z0p14')
#uchuu_dir_z0p19 = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/WPRP', 'z0p19')
#uchuu_dir_z0p25 = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/WPRP', 'z0p25')
def get_data(file_list):
    t = Table.read(file_list[0])
    all_wp = np.array([ Table.read(el)['wprp'] for el in file_list ])
    t['wprp_JK_mean'] = np.mean(all_wp, axis=0)
    t['wprp_JK_std'] = np.std(all_wp, axis=0)
    return t

U0p14["BC_10.25_LX42.8"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_BC_HaloLX_gt_42.8_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p14["RS_10.25_LX42.8"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_RS_HaloLX_gt_42.8_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p14["ANY_10.25_LX42.8"] = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_HaloLX_gt_42.8_replication_*-wprp-pimax100-2pcf-bin0p1.fits"    )))
U0p14["BC_10.25_LX42.9"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_BC_HaloLX_gt_42.9_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p14["RS_10.25_LX42.9"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_RS_HaloLX_gt_42.9_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p14["ANY_10.25_LX42.9"] = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_HaloLX_gt_42.9_replication_*-wprp-pimax100-2pcf-bin0p1.fits"    )))


U0p19[ "BC_10.75_LX43.1"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_BC_HaloLX_gt_43.1_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p19[ "RS_10.75_LX43.1"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_RS_HaloLX_gt_43.1_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p19["ANY_10.75_LX43.1"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_HaloLX_gt_43.1_replication_*-wprp-pimax100-2pcf-bin0p1.fits"    )))
U0p19[ "BC_10.75_LX43.2"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_BC_HaloLX_gt_43.2_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p19[ "RS_10.75_LX43.2"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_RS_HaloLX_gt_43.2_replication_*-wprp-pimax100-2pcf-bin0p1.fits" )))
U0p19["ANY_10.75_LX43.2"]  = get_data(glob.glob(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.75_HaloLX_gt_43.2_replication_*-wprp-pimax100-2pcf-bin0p1.fits"    )))

#U0p14["ANY_10.25"] = Table.read(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                   ))
#U0p14["BC_10.25"] = Table.read(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_BC_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))
#U0p14["RS_10.25"] = Table.read(os.path.join(uchuu_dir_z0p14, "Ms_gt_10.25_RS_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))

#U0p19["BC_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_BC_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
#U0p19["RS_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_RS_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
#U0p19["ANY_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"    ))
#U0p19["ANY_10.75"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                   ))
#U0p19["RS_10.75"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_RS_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))
#U0p19["BC_10.75"] = Table.read(os.path.join(uchuu_dir_z0p19, "Ms_gt_10.75_BC_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))

#U0p25["BC_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_BC_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
#U0p25["RS_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_RS_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits" ))
#U0p25["ANY_10.75_LX43.1"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_HaloLX_gt_43.1_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"    ))
#U0p25["ANY_10.75"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                   ))
#U0p25["RS_10.75"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_RS_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))
#U0p25["BC_10.75"] = Table.read(os.path.join(uchuu_dir_z0p25, "Ms_gt_10.75_BC_replication_-1.0_-1.0_0.0-wprp-pimax100-2pcf-bin0p1.fits"                ))

lgal_dir_z0p2 = '../data/lgal'
LGAL = {}
#LGAL["ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_Ms_1075_JK_2.fits" ))
#LGAL["LX43.1_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_43.1_Ms_1075_JK_2.fits" ))
#LGAL["LX43.2_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_43.2_Ms_1075_JK_2.fits" ))
#LGAL["LX43.3_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_43.3_Ms_1075_JK_2.fits" ))
#LGAL["LX43.4_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_43.4_Ms_1075_JK_2.fits" ))
#LGAL["LX43.5_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_43.5_Ms_1075_JK_2.fits" ))
#LGAL["LX43.6_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_43.6_Ms_1075_JK_2.fits" ))
#LGAL["LX43.7_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_43.7_Ms_1075_JK_2.fits" ))
#LGAL["LX43.8_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_43.8_Ms_1075_JK_2.fits" ))
#LGAL["LX43.9_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_43.9_Ms_1075_JK_2.fits" ))
#LGAL["LX44.0_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.0_Ms_1075_JK_2.fits" ))
#LGAL["LX44.1_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.1_Ms_1075_JK_2.fits" ))
#LGAL["LX44.2_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.2_Ms_1075_JK_2.fits" ))
#LGAL["LX44.3_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.3_Ms_1075_JK_2.fits" ))
#LGAL["LX44.4_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.4_Ms_1075_JK_2.fits" ))
#LGAL["LX44.5_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.5_Ms_1075_JK_2.fits" ))
#LGAL["LX44.6_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.6_Ms_1075_JK_2.fits" ))
#LGAL["LX44.7_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.7_Ms_1075_JK_2.fits" ))
#LGAL["LX44.8_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.8_Ms_1075_JK_2.fits" ))
#LGAL["LX44.9_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_44.9_Ms_1075_JK_2.fits" ))
#LGAL["LX45.0_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_45.0_Ms_1075_JK_2.fits" ))
#LGAL["LX45.1_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_45.1_Ms_1075_JK_2.fits" ))
#LGAL["LX45.2_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_45.2_Ms_1075_JK_2.fits" ))
#LGAL["LX45.3_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_45.3_Ms_1075_JK_2.fits" ))
#LGAL["LX45.4_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_45.4_Ms_1075_JK_2.fits" ))
#LGAL["LX45.5_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_LX_45.5_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.0_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.0_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.1_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.1_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.2_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.2_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.3_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.3_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.4_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.4_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.5_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.5_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.6_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.6_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.7_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.7_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.8_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.8_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.9_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.9_Ms_1075_JK_2.fits" ))
#LGAL["Mh15.0_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_15.0_Ms_1075_JK_2.fits" ))
#LGAL["Mh15.1_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_15.1_Ms_1075_JK_2.fits" ))
#LGAL["Mh15.2_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_15.2_Ms_1075_JK_2.fits" ))
#LGAL["Mh15.3_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_15.3_Ms_1075_JK_2.fits" ))
#LGAL["Mh15.4_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_15.4_Ms_1075_JK_2.fits" ))
#LGAL["Mh15.5_ANY_10.75"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_15.5_Ms_1075_JK_2.fits" ))
#LGAL["Mh14.0_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.0_Ms_1067_JK_100.fits" ))
#LGAL["Mh14.1_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.1_Ms_1067_JK_100.fits" ))
#LGAL["Mh14.2_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.2_Ms_1067_JK_100.fits" ))
#LGAL["Mh14.3_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.3_Ms_1067_JK_100.fits" ))
#LGAL["Mh14.4_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.4_Ms_1067_JK_100.fits" ))
#LGAL["Mh14.5_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.5_Ms_1067_JK_100.fits" ))
#LGAL["Mh14.6_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.6_Ms_1067_JK_100.fits" ))
#LGAL["Mh14.7_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.7_Ms_1067_JK_100.fits" ))
#LGAL["Mh14.8_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.8_Ms_1067_JK_100.fits" ))
#LGAL["Mh14.9_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_14.9_Ms_1067_JK_100.fits" ))
#LGAL["Mh15.0_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_15.0_Ms_1067_JK_100.fits" ))
#LGAL["Mh15.1_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_15.1_Ms_1067_JK_100.fits" ))
#LGAL["Mh15.2_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_M200c_15.2_Ms_1067_JK_100.fits" ))
LGAL["L0520_43.6_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.6_Ms_1067_JK_100.fits" ))
LGAL["L0520_43.8_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.8_Ms_1067_JK_100.fits" ))
LGAL["L0520_44.0_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.0_Ms_1067_JK_100.fits" ))
LGAL["L0520_44.2_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.2_Ms_1067_JK_100.fits" ))
LGAL["L0520_44.4_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.4_Ms_1067_JK_100.fits" ))

#
#
# S1 10.75
#
#

# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-S1-clusters-LGAL.png')
plt.figure(10, (5,5))

t_wp = CxG["S1_ANY_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='Galaxies x Clusters ', color='grey')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='grey', alpha=0.4)


for kk, cc in zip( np.array(list(LGAL.keys())), colors):
    t_wp = LGAL[kk]
    plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp'], lw=2,  ls='solid', label='LGAL LX>' +kk[6:10], color=cc)
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color=cc, alpha=0.3)

plt.ylim((1, 1500))
plt.xlim((0.05, 60))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$")#, $\pi_{max}=100$ Mpc$/h$")
plt.legend(loc=0, fontsize=8,ncol=1, title=r'$0.1<z<0.3$')#, title='LS10, r<19.5, 0.05<z<0.22')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


LGAL = {}
#LGAL["L0520_42.7_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.7_BC_1067_JK_100.fits" ))
#LGAL["L0520_42.7_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.7_Ms_1067_JK_100.fits" ))
#LGAL["L0520_42.7_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.7_RS_1067_JK_100.fits" ))
#LGAL["L0520_42.8_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.8_BC_1067_JK_100.fits" ))
#LGAL["L0520_42.8_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.8_Ms_1067_JK_100.fits" ))
#LGAL["L0520_42.8_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.8_RS_1067_JK_100.fits" ))
#LGAL["L0520_42.9_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.9_BC_1067_JK_100.fits" ))
#LGAL["L0520_42.9_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.9_Ms_1067_JK_100.fits" ))
#LGAL["L0520_42.9_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.9_RS_1067_JK_100.fits" ))
#LGAL["L0520_43.0_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.0_BC_1067_JK_100.fits" ))
#LGAL["L0520_43.0_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.0_Ms_1067_JK_100.fits" ))
#LGAL["L0520_43.0_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.0_RS_1067_JK_100.fits" ))

#LGAL["L0520_43.1_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.1_BC_1067_JK_100.fits" ))
#LGAL["L0520_43.1_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.1_Ms_1067_JK_100.fits" ))
#LGAL["L0520_43.1_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.1_RS_1067_JK_100.fits" ))

#LGAL["L0520_43.2_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.2_BC_1067_JK_100.fits" ))
#LGAL["L0520_43.2_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.2_Ms_1067_JK_100.fits" ))
#LGAL["L0520_43.2_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.2_RS_1067_JK_100.fits" ))
#LGAL["L0520_43.3_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.3_BC_1067_JK_100.fits" ))
#LGAL["L0520_43.3_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.3_Ms_1067_JK_100.fits" ))
#LGAL["L0520_43.3_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.3_RS_1067_JK_100.fits" ))
#LGAL["L0520_43.4_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.4_BC_1067_JK_100.fits" ))
#LGAL["L0520_43.4_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.4_Ms_1067_JK_100.fits" ))
#LGAL["L0520_43.4_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.4_RS_1067_JK_100.fits" ))
#LGAL["L0520_43.5_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.5_BC_1067_JK_100.fits" ))
#LGAL["L0520_43.5_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.5_Ms_1067_JK_100.fits" ))
#LGAL["L0520_43.5_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.5_RS_1067_JK_100.fits" ))
#LGAL["L0520_43.6_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.6_BC_1067_JK_100.fits" ))
#LGAL["L0520_43.6_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.6_Ms_1067_JK_100.fits" ))
#LGAL["L0520_43.6_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.6_RS_1067_JK_100.fits" ))
#LGAL["L0520_43.7_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.7_BC_1067_JK_100.fits" ))
#LGAL["L0520_43.7_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.7_Ms_1067_JK_100.fits" ))
#LGAL["L0520_43.7_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.7_RS_1067_JK_100.fits" ))
#LGAL["L0520_43.8_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.8_BC_1067_JK_100.fits" ))
#LGAL["L0520_43.8_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.8_Ms_1067_JK_100.fits" ))
#LGAL["L0520_43.8_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.8_RS_1067_JK_100.fits" ))
LGAL["L0520_43.9_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.9_BC_1067_JK_100.fits" ))
LGAL["L0520_43.9_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.9_Ms_1067_JK_100.fits" ))
LGAL["L0520_43.9_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.9_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.0_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.0_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.0_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.0_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.0_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.0_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.1_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.1_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.1_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.1_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.1_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.1_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.2_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.2_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.2_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.2_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.2_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.2_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.3_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.3_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.3_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.3_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.3_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.3_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.4_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.4_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.4_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.4_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.4_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.4_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.5_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.5_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.5_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.5_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.5_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.5_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.6_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.6_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.6_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.6_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.6_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.6_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.7_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.7_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.7_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.7_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.7_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.7_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.8_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.8_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.8_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.8_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.8_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.8_RS_1067_JK_100.fits" ))
#LGAL["L0520_44.9_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.9_BC_1067_JK_100.fits" ))
#LGAL["L0520_44.9_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.9_Ms_1067_JK_100.fits" ))
#LGAL["L0520_44.9_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.9_RS_1067_JK_100.fits" ))
#LGAL["L0520_45.0_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.0_BC_1067_JK_100.fits" ))
#LGAL["L0520_45.0_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.0_Ms_1067_JK_100.fits" ))
#LGAL["L0520_45.0_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.0_RS_1067_JK_100.fits" ))
#LGAL["L0520_45.1_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.1_BC_1067_JK_100.fits" ))
#LGAL["L0520_45.1_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.1_Ms_1067_JK_100.fits" ))
#LGAL["L0520_45.1_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.1_RS_1067_JK_100.fits" ))
#LGAL["L0520_45.2_BC_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.2_BC_1067_JK_100.fits" ))
#LGAL["L0520_45.2_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.2_Ms_1067_JK_100.fits" ))
#LGAL["L0520_45.2_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.2_RS_1067_JK_100.fits" ))
#LGAL["L0520_45.3_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.3_Ms_1067_JK_100.fits" ))
#LGAL["L0520_45.3_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.3_RS_1067_JK_100.fits" ))
#LGAL["L0520_45.4_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.4_Ms_1067_JK_100.fits" ))
#LGAL["L0520_45.4_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.4_RS_1067_JK_100.fits" ))
#LGAL["L0520_45.5_Ms_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.5_Ms_1067_JK_100.fits" ))
#LGAL["L0520_45.5_RS_1067"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_45.5_RS_1067_JK_100.fits" ))

#for lmin in np.arange(42.7, 45.3,0.1):
for lmin in np.arange(43.9, 43.91,0.1):
    l_str = str(np.round(lmin,1))
    # plot wprp
    p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-RSBC-S1-clusters-LGAL_LX_'+l_str+'.png')
    plt.figure(10, (5,5))

    t_wp = CxG["S1_ANY_10.75"]
    plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='Galaxies x Clusters ', color='grey')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='grey', alpha=0.4)
    t_wp = CxG["S1_BC_10.75"]
    plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=1,  ls='solid', label='blue-cloud', color='darkblue')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkblue', alpha=0.4)
    t_wp = CxG["S1_RS_10.75"]
    plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=1,  ls='solid', label='red-sequence', color='darkred')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkred', alpha=0.4)


    t_wp = LGAL["L0520_"+l_str+"_Ms_1067"]
    t_wp['wprp'] = t_wp['wprp']*1.5
    plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='k', label='LGAL')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']+0.05
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='k', alpha=0.3)
    t_wp = LGAL["L0520_"+l_str+"_BC_1067"]
    t_wp['wprp'] = t_wp['wprp']*1.5
    plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='magenta')#, label='LGAL sSFR>-11')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']+0.05
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='magenta', alpha=0.3)
    t_wp = LGAL["L0520_"+l_str+"_RS_1067"]
    t_wp['wprp'] = t_wp['wprp']*1.5
    plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='orange')#, label='LGAL sSFR<-11')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']+0.05
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='orange', alpha=0.3)

    plt.ylim((1, 1500))
    plt.xlim((0.05, 60))
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel(r"$r_p$ [Mpc/h]")
    plt.ylabel(r"$r_p\times w_p(r_p)$")#, $\pi_{max}=100$ Mpc$/h$")
    plt.legend(loc=0, fontsize=12,ncol=1, title='C1 x G1075')#, LX>'+l_str)#, title='LS10, r<19.5, 0.05<z<0.22')
    plt.tight_layout()
    plt.savefig(p2_fig)
    plt.clf()
    print(p2_fig)



LGAL = {}
#LGAL["L0520_42.7_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.7_BC_1023_JK_100.fits" ))
#LGAL["L0520_42.7_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.7_Ms_1023_JK_100.fits" ))
#LGAL["L0520_42.7_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.7_RS_1023_JK_100.fits" ))
#LGAL["L0520_42.8_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.8_BC_1023_JK_100.fits" ))
#LGAL["L0520_42.8_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.8_Ms_1023_JK_100.fits" ))
#LGAL["L0520_42.8_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.8_RS_1023_JK_100.fits" ))
#LGAL["L0520_42.9_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.9_BC_1023_JK_100.fits" ))
#LGAL["L0520_42.9_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.9_Ms_1023_JK_100.fits" ))
#LGAL["L0520_42.9_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_42.9_RS_1023_JK_100.fits" ))
#LGAL["L0520_43.0_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.0_BC_1023_JK_100.fits" ))
#LGAL["L0520_43.0_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.0_Ms_1023_JK_100.fits" ))
#LGAL["L0520_43.0_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.0_RS_1023_JK_100.fits" ))
#LGAL["L0520_43.1_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.1_BC_1023_JK_100.fits" ))
#LGAL["L0520_43.1_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.1_Ms_1023_JK_100.fits" ))
#LGAL["L0520_43.1_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.1_RS_1023_JK_100.fits" ))
#LGAL["L0520_43.2_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.2_BC_1023_JK_100.fits" ))
#LGAL["L0520_43.2_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.2_Ms_1023_JK_100.fits" ))
#LGAL["L0520_43.2_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.2_RS_1023_JK_100.fits" ))
#LGAL["L0520_43.3_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.3_BC_1023_JK_100.fits" ))
#LGAL["L0520_43.3_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.3_Ms_1023_JK_100.fits" ))
#LGAL["L0520_43.3_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.3_RS_1023_JK_100.fits" ))
#LGAL["L0520_43.4_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.4_BC_1023_JK_100.fits" ))
#LGAL["L0520_43.4_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.4_Ms_1023_JK_100.fits" ))
#LGAL["L0520_43.4_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.4_RS_1023_JK_100.fits" ))
#LGAL["L0520_43.5_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.5_BC_1023_JK_100.fits" ))
#LGAL["L0520_43.5_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.5_Ms_1023_JK_100.fits" ))
#LGAL["L0520_43.5_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.5_RS_1023_JK_100.fits" ))
#LGAL["L0520_43.6_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.6_BC_1023_JK_100.fits" ))
#LGAL["L0520_43.6_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.6_Ms_1023_JK_100.fits" ))
#LGAL["L0520_43.6_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.6_RS_1023_JK_100.fits" ))
#LGAL["L0520_43.7_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.7_BC_1023_JK_100.fits" ))
#LGAL["L0520_43.7_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.7_Ms_1023_JK_100.fits" ))
#LGAL["L0520_43.7_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.7_RS_1023_JK_100.fits" ))
#LGAL["L0520_43.8_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.8_BC_1023_JK_100.fits" ))
#LGAL["L0520_43.8_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.8_Ms_1023_JK_100.fits" ))
#LGAL["L0520_43.8_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.8_RS_1023_JK_100.fits" ))
LGAL["L0520_43.9_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.9_BC_1023_JK_100.fits" ))
LGAL["L0520_43.9_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.9_Ms_1023_JK_100.fits" ))
LGAL["L0520_43.9_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_43.9_RS_1023_JK_100.fits" ))
#LGAL["L0520_44.0_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.0_BC_1023_JK_100.fits" ))
#LGAL["L0520_44.0_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.0_Ms_1023_JK_100.fits" ))
#LGAL["L0520_44.0_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.0_RS_1023_JK_100.fits" ))
#LGAL["L0520_44.1_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.1_BC_1023_JK_100.fits" ))
#LGAL["L0520_44.1_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.1_Ms_1023_JK_100.fits" ))
#LGAL["L0520_44.1_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.1_RS_1023_JK_100.fits" ))
#LGAL["L0520_44.2_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.2_BC_1023_JK_100.fits" ))
#LGAL["L0520_44.2_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.2_Ms_1023_JK_100.fits" ))
#LGAL["L0520_44.2_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.2_RS_1023_JK_100.fits" ))
#LGAL["L0520_44.3_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.3_BC_1023_JK_100.fits" ))
#LGAL["L0520_44.3_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.3_Ms_1023_JK_100.fits" ))
#LGAL["L0520_44.3_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.3_RS_1023_JK_100.fits" ))
#LGAL["L0520_44.4_BC_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.4_BC_1023_JK_100.fits" ))
#LGAL["L0520_44.4_Ms_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.4_Ms_1023_JK_100.fits" ))
#LGAL["L0520_44.4_RS_1023"] = Table.read(os.path.join(lgal_dir_z0p2, "LGAL_01z03_L0520_44.4_RS_1023_JK_100.fits" ))

#for lmin in np.arange(42.7, 44.4,0.1):
for lmin in np.arange(43.9, 43.91,0.1):
    l_str = str(np.round(lmin,1))
    # plot wprp
    p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-RSBC-S1-clusters-LGAL_LX_'+l_str+'.png')
    plt.figure(10, (5,5))

    t_wp = CxG["S0_ANY_10.25"]
    plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='Galaxies x Clusters ', color='grey')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='grey', alpha=0.4)
    t_wp = CxG["S0_BC_10.25"]
    plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=1,  ls='solid', label='blue-cloud', color='darkblue')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkblue', alpha=0.4)
    t_wp = CxG["S0_RS_10.25"]
    plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=1,  ls='solid', label='red-sequence', color='darkred')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkred', alpha=0.4)


    t_wp = LGAL["L0520_"+l_str+"_Ms_1023"]
    plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='k', label='LGAL')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']+0.05
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='k', alpha=0.3)
    t_wp = LGAL["L0520_"+l_str+"_BC_1023"]
    plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='magenta')#, label='LGAL sSFR>-11')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']+0.05
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='magenta', alpha=0.3)
    t_wp = LGAL["L0520_"+l_str+"_RS_1023"]
    plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='orange')#, label='LGAL sSFR<-11')
    f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']+0.05
    plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='orange', alpha=0.3)

    plt.ylim((1, 1500))
    plt.xlim((0.05, 60))
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel(r"$r_p$ [Mpc/h]")
    plt.ylabel(r"$r_p\times w_p(r_p)$")#, $\pi_{max}=100$ Mpc$/h$")
    plt.legend(loc=0, fontsize=12,ncol=1, title='C0 x G1025')#, LX>'+l_str)#, title='LS10, r<19.5, 0.05<z<0.22')
    plt.tight_layout()
    plt.savefig(p2_fig)
    plt.clf()
    print(p2_fig)


## plot wprp
#p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-RSBC-S1-clusters-LGAL_LX438.png')
#plt.figure(10, (5,5))

#t_wp = CxG["S1_ANY_10.75"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='Galaxies x Clusters ', color='grey')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='grey', alpha=0.4)
#t_wp = CxG["S1_BC_10.75"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=1,  ls='solid', label='blue-cloud', color='darkblue')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkblue', alpha=0.4)
#t_wp = CxG["S1_RS_10.75"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=1,  ls='solid', label='red-sequence', color='darkred')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkred', alpha=0.4)


#t_wp = LGAL["L0520_44.0_Ms_1067"]
#plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='k', label='LGAL z=0.2')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='k', alpha=0.3)
#t_wp = LGAL["L0520_44.0_BC_1067"]
#plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='magenta', label='LGAL sSFR>-11')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='magenta', alpha=0.3)
#t_wp = LGAL["L0520_44.0_RS_1067"]
#plt.plot(t_wp['rp_mid'], t_wp['rp_mid']*t_wp['wprp_JK_mean'], lw=2,  ls='dashed', color='orange', label='LGAL sSFR<-11')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1-f_err), t_wp['wprp_JK_mean']*t_wp['rp_mid']*(1+f_err), color='orange', alpha=0.3)

#plt.ylim((1, 1500))
#plt.xlim((0.05, 60))
#plt.xscale('log')
##plt.yscale('log')
#plt.xlabel(r"$r_p$ [Mpc/h]")
#plt.ylabel(r"$r_p\times w_p(r_p)$")#, $\pi_{max}=100$ Mpc$/h$")
#plt.legend(loc=0, fontsize=8,ncol=1, title=r'$0.1<z<0.3$')#, title='LS10, r<19.5, 0.05<z<0.22')
#plt.tight_layout()
#plt.savefig(p2_fig)
#plt.clf()
#print(p2_fig)

sys.exit()

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
gal = open_hdf5_dict(filename)


#Millennium (Planck1) 0.315 0.049 0.685 67.3 0.96 0.826 21603 1.43 × 109 714 64 56
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
import astropy.units as u
cosmo = FlatLambdaCDM(H0=67.3 * u.km / u.s / u.Mpc, Om0=0.315)
h = 0.673
L_box = 714.0
z_array = np.arange(0, 7.5, 0.001)
dc_to_z = interp1d(cosmo.comoving_distance(z_array), z_array)


def get_radecz(x,y,z):
	rr = (x**2 + y**2 + z**2)**0.5
	theta = np.arccos(z / rr) * 180 / np.pi
	phi = np.arctan2(y, x) * 180 / np.pi
	ra = phi + 180.
	dec = theta - 90.
	redshift = dc_to_z(rr)
	return ra, dec, redshift



Xg=gal['Pos'].T[0]/1000.
Yg=gal['Pos'].T[1]/1000.
Zg=gal['Pos'].T[2]/1000.
ra_g, dec_g, redshift_g = get_radecz(Xg, Yg, Zg)


# plot wprp
p2_fig = os.path.join( fig_dir, 'M200c-LX-LGAL.png')
plt.figure(10, (5,5))
plt.plot(np.log10(gal['M_Crit200']), gal['XrayLum'], 'k,')
plt.xlabel(r"log10 M_Crit200")
plt.ylabel(r"XrayLum")
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)
p2_fig = os.path.join( fig_dir, 'M200c-LX-LGAL-type0.png')
plt.figure(10, (5,5))
plt.plot(np.log10(gal['M_Crit200'][gal['Type']==0]), gal['XrayLum'][gal['Type']==0], 'k,')
plt.xlabel(r"log10 M_Crit200")
plt.ylabel(r"XrayLum")
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

p2_fig = os.path.join( fig_dir, 'M200c-hist-LGAL-types.png')
plt.figure(10, (5,5))
plt.hist(np.log10(gal['M_Crit200'][gal['Type']==0]), histtype='step', label="type 0")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.hist(np.log10(gal['M_Crit200'][gal['Type']==1]), histtype='step', label="type 1")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.hist(np.log10(gal['M_Crit200'][gal['Type']==2]), histtype='step', label="type 2")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.xlabel(r"log10 M_Crit200")
plt.ylabel(r"XrayLum")
plt.legend(loc=0)
plt.yscale('log')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'M200c-hist-LGAL-type0.png')
plt.figure(10, (5,5))
plt.hist(np.log10(gal['M_Crit200'][gal['Type']==0]), histtype='step', label="type 0")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.hist(np.log10(gal['M_Crit200'][(gal['Type']==0)&(gal['XrayLum']==0)]), histtype='step', label="type 0 LX=0")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.hist(np.log10(gal['M_Crit200'][(gal['Type']==0)&(gal['XrayLum']>0)]), histtype='step', label="type 0, LX>0")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.xlabel(r"log10 M_Crit200")
plt.ylabel(r"XrayLum")
plt.legend(loc=0)
plt.yscale('log')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'M200c-hist-LGAL-type1.png')
plt.figure(10, (5,5))
plt.hist(np.log10(gal['M_Crit200'][gal['Type']==1]), histtype='step', label="type 1")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.hist(np.log10(gal['M_Crit200'][(gal['Type']==1)&(gal['XrayLum']==0)]), histtype='step', label="type 1 LX=0")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.hist(np.log10(gal['M_Crit200'][(gal['Type']==1)&(gal['XrayLum']>0)]), histtype='step', label="type 1, LX>0")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.xlabel(r"log10 M_Crit200")
plt.ylabel(r"XrayLum")
plt.legend(loc=0)
plt.yscale('log')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'M200c-hist-LGAL-type2.png')
plt.figure(10, (5,5))
plt.hist(np.log10(gal['M_Crit200'][gal['Type']==2]), histtype='step', label="type 2")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.hist(np.log10(gal['M_Crit200'][(gal['Type']==2)&(gal['XrayLum']==0)]), histtype='step', label="type 2 LX=0")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.hist(np.log10(gal['M_Crit200'][(gal['Type']==2)&(gal['XrayLum']>0)]), histtype='step', label="type 2, LX>0")#, gal['XrayLum'][gal['Type']==0], 'k,')
plt.xlabel(r"log10 M_Crit200")
plt.ylabel(r"XrayLum")
plt.legend(loc=0)
plt.yscale('log')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


import seaborn as sns
sns.set_theme(style="ticks")
import pandas as pd
df = pd.DataFrame()
df['M_Crit200'] = np.log10(gal['M_Crit200'])
df['XrayLum'] = gal['XrayLum']
df['Type'] = gal['Type']
# Load the penguins dataset
penguins = sns.load_dataset("penguins")

# Show the joint distribution using kernel density estimation
p2_fig = os.path.join( fig_dir, 'M200c-LX-LGAL-all_types.png')
plt.figure(10, (5,5))
g = sns.jointplot(
    data=df,
    x="M_Crit200", y="XrayLum", hue="Type",
    kind="kde",
)
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


