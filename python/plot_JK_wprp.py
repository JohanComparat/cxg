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

d1025 = '../data/_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100'
d1075 = '../data/_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100'
dC0xG1025 = '../data/C0xG1025'
dC1xG1075 = '../data/C1xG1075'

# "C:\Users\Johan Comparat\Documents\Shared\software\st_mod_data\data"
# export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
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
BGS["RS_10.25"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.25-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
BGS["RS_10.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-MstarMin_10.5-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["RS_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_10.75-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
BGS["RS_11.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["RS_11.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p05-HpxMask.fits" ) )

BGS["BC_10.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_10.0-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
BGS["BC_10.25"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.25-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
BGS["BC_10.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-MstarMin_10.5-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["BC_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_10.75-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
BGS["BC_11.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["BC_11.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p05-HpxMask.fits" ) )

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
CxG["S0_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S0_ANY_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S0_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S0_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S0_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )

CxG["S1_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S1_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S1_ANY_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S1_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S1_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )

CxG["S2_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S2_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S2_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S2_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )

CxG["S3_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask_JK.fits"  ) )
CxG["S3_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S3_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S3_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S3_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask.fits"  ) )

CxG["S4_ANY_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S4_ANY_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S4_ANY_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits"  ) )
CxG["S4_ANY_11.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S4_CROSS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p1-HpxMask.fits"  ) )

CxG["S2_BC_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_BC_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask.fits" ) )
CxG["S2_RS_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_RS_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

CxG["S1_BC_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_BC_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits" ) )
CxG["S1_RS_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_RS_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits" ) )

CxG["S0_BC_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_BC_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits" ) )
CxG["S0_RS_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_RS_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits" ) )



CxG["S1_BC_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
CxG["S1_RS_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

CxG["S0_BC_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
CxG["S0_RS_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

pcf_list_1025 = np.array( glob.glob( os.path.join( d1025, '*NSIDE_02*.fits' )))

z13_r, z13_wp, z13_wp_up = np.loadtxt('../data/zu-weinberg-fig13-red-datapoints-M1140-M1190.csv', unpack=True, delimiter=',')
z13_wp_ferr = z13_wp_up/z13_wp-1
#z13M1120_r, z13M1120_wp = np.loadtxt('../data/zu-weinberg-fig13-green-datapoints-M1120-M1140.csv', unpack=True, delimiter=',')
#z13_x = z13_r[np.argsort(z13_r)]
#z13_y = z13_xi[np.argsort(s06_r)]
#
#
# S1 10.75
#
#

# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS-S1-clusters-err-test.png')
plt.figure(13, (8., 6.))


pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_ANY_*NSIDE_04*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='solid', label='C1 x G1075, N4='+str(len(pcf_list_1075)), color='black')

pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_ANY_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='dotted', label='C1 x G1075, N8='+str(len(pcf_list_1075)), color='black')


pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_BC_*NSIDE_04*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='solid', label='C1 x G1075 BC, N4='+str(len(pcf_list_1075)), color='b')

pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_BC_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='dotted', label='C1 x G1075 BC, N8='+str(len(pcf_list_1075)), color='b')


pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_RS_*NSIDE_04*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='solid', label='C1 x G1075 RS, N4='+str(len(pcf_list_1075)), color='r')

pcf_list_1075 = np.array( glob.glob( os.path.join( dC1xG1075, '*_RS_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='dotted', label='C1 x G1075 RS, N8='+str(len(pcf_list_1075)), color='r')


t_wp = CxG["S1_ANY_10.75"]
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='dashed', color='black') # , label='C0 x G1025'
t_wp = CxG["S1_BC_10.75"]
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='blue-cloud', color='darkblue')
t_wp = CxG["S1_RS_10.75"]
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='red-sequence', color='darkred')

plt.plot(z13_r, z13_wp_ferr, lw=2, label='Zu 13 (M$^*\sim11.7$)', ls='dotted', color='green', zorder=0)

# plt.ylim((0.01, 0.3))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$\sigma(w_p)/w_p$")#")
plt.legend(loc=0, fontsize=10,ncol=2)
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)



# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS-S0-clusters-err-test.png')
plt.figure(13, (8., 6.))

pcf_list_1075 = np.array( glob.glob( os.path.join( dC0xG1025, '*_ANY_*NSIDE_04*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='solid', label='C0 x G1025, N4='+str(len(pcf_list_1075)), color='black')

pcf_list_1075 = np.array( glob.glob( os.path.join( dC0xG1025, '*_ANY_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='dotted', label='C0 x G1025, N8='+str(len(pcf_list_1075)), color='black')


pcf_list_1075 = np.array( glob.glob( os.path.join( dC0xG1025, '*_BC_*NSIDE_04*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='solid', label='C0 x G1025 BC, N4='+str(len(pcf_list_1075)), color='b')

pcf_list_1075 = np.array( glob.glob( os.path.join( dC0xG1025, '*_BC_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='dotted', label='C0 x G1025 BC, N8='+str(len(pcf_list_1075)), color='b')

pcf_list_1075 = np.array( glob.glob( os.path.join( dC0xG1025, '*_RS_*NSIDE_04*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='solid', label='C0 x G1025 RS, N4='+str(len(pcf_list_1075)), color='r')

pcf_list_1075 = np.array( glob.glob( os.path.join( dC0xG1025, '*_RS_*NSIDE_08*.fits' )))
all_wprp = []
for el in pcf_list_1075 :
	all_wprp.append(Table.read(el)['wprp'])
all_wprp = np.array(all_wprp)
t_ref = Table.read(el)
plt.plot(t_ref['rp_mid'], all_wprp.std(axis=0)/all_wprp.mean(axis=0), lw=2,  ls='dotted', label='C0 x G1025 RS, N8='+str(len(pcf_list_1075)), color='r')


t_wp = CxG["S0_ANY_10.25"]
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='dashed', color='black') # , label='C0 x G1025'
t_wp = CxG["S0_BC_10.25"]
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='blue-cloud', color='darkblue')
t_wp = CxG["S0_RS_10.25"]
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.plot(t_wp['rp_mid'], f_err, lw=2,  ls='solid', label='red-sequence', color='darkred')

plt.plot(z13_r, z13_wp_ferr, label='Zu 13 (SDSS)', ls='dotted', color='green', zorder=0)

# plt.ylim((0.01, 0.3))
plt.xlim((0.03, 60))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$\sigma(w_p)/w_p$")#")
plt.legend(loc=0, fontsize=10,ncol=2)#, title=r'$0.1<z<0.2$, ' '$\log_{10}(L_X\; [\mathrm{erg/s}])>42.7$ x'+ "\n" +r'$10.25<\log_{10}(M*[M_\odot])<12$')#, title='LS10, r<19.5, 0.05<z<0.22')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


