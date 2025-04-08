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

fig_dir  ='../figures/'

ZuMa_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/benchmark/zu-mandelbaum-1505.02781v1')

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

LS10_bgs_vlim_dir = os.path.join(os.environ['LSDR10'], 'sweep/BGS_VLIM_Mstar')
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

BGS["RS_10.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_10.0-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["RS_10.25"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.25-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["RS_10.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-MstarMin_10.5-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["RS_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_10.75-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["RS_11.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["RS_11.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "RS_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p05-HpxMask.fits" ) )

BGS["BC_10.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-MstarMin_10.0-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["BC_10.25"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-MstarMin_10.25-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["BC_10.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-MstarMin_10.5-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["BC_10.75"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-MstarMin_10.75-wprp-pimax100-bin0p1-HpxMask.fits" ) )
BGS["BC_11.0"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p05-HpxMask.fits" ) )
BGS["BC_11.5"]  = Table.read( os.path.join(LS10_bgs_vlim_dir, "BC_LS10_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_N_0120882-wprp-pimax100-bin0p05-HpxMask.fits" ) )


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M100-ANY-BC-RS.png')
plt.figure(12, (5,5))
t_wp = BGS["ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='10<M*<12, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='dashed', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='dashed', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)$ [Mpc/h]")
plt.legend(loc=3, fontsize=12, title='LS10, r<19.5, 0.05<z<0.18')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS.png')
plt.figure(11, (6,5))
t_wp = BGS["ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='All, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')
plt.ylim((0.5, 1e4))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)$ [Mpc/h]")
plt.legend(loc=3, fontsize=12, title='0.05<z<0.22,\n' +r'10.25$<\log_{10}(M*/M_\odot)<$12')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS.png')
plt.figure(11, (6,5))
t_wp = BGS["ANY_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='All, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')
plt.ylim((0.5, 1e4))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)$ [Mpc/h]")
plt.legend(loc=3, fontsize=12, title='0.05<z<0.31,\n'+r'10.75$<\log_{10}(M*/M_\odot)<$12')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)

# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M105-ANY-BC-RS.png')
plt.figure(12, (5,5))
t_wp = BGS["ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='10.5<M*<12, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='dashed', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='dashed', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)$ [Mpc/h]")
plt.legend(loc=3, fontsize=12, title='LS10, r<19.5, 0.05<z<0.26')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M110-ANY-BC-RS.png')
plt.figure(12, (5,5))
t_wp = BGS["ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='11.0<M*<12, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='dashed', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='dashed', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)$ [Mpc/h]")
plt.legend(loc=3, fontsize=12, title='LS10, r<19.5, 0.05<z<0.35')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M115-ANY-BC-RS.png')
plt.figure(12, (5,5))
t_wp = BGS["ANY_11.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='11.5<M*<12, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_11.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='dashed', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_11.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='dashed', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w_p(r_p)$ [Mpc/h]")
plt.legend(loc=3, fontsize=12, title='LS10, r<19.5, 0.05<z<0.35')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


#uchuu_dir = os.path.join(os.environ['UCHUU'], "FullSky/mockBGS/replication_0.0_0.0_0.0" )
#UCHUU = {}
##UCHUU['ALL_9.0'] = Table.read(os.path.join(uchuu_dir, "UCHUU_VLIM_ANY_9.5_Mstar_12.0_0.05_z_0.12_0err_glist-wprp-pimax100-bin0p1.fits" ))
##UCHUU['ALL_9.5'] = Table.read(os.path.join(uchuu_dir, "UCHUU_VLIM_ANY_9.0_Mstar_12.0_0.05_z_0.08_0err_glist-wprp-pimax100-bin0p1.fits" ))
#UCHUU['ALL_10.0'] = Table.read(os.path.join(uchuu_dir, "UCHUU_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_0err_glist-wprp-pimax100-bin0p1.fits" ))
#UCHUU['ALL_10.25'] = Table.read(os.path.join(uchuu_dir, "UCHUU_VLIM_ANY_10.25_Mstar_12.0_0.05_z_0.22_0err_glist-wprp-pimax100-bin0p1.fits" ))
#UCHUU['ALL_10.5'] = Table.read(os.path.join(uchuu_dir, "UCHUU_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_0err_glist-wprp-pimax100-bin0p1.fits" ))
#UCHUU['ALL_10.75'] = Table.read(os.path.join(uchuu_dir, "UCHUU_VLIM_ANY_10.75_Mstar_12.0_0.05_z_0.31_0err_glist-wprp-pimax100-bin0p1.fits" ))
#UCHUU['ALL_11.0'] = Table.read(os.path.join(uchuu_dir, "UCHUU_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_0err_glist-wprp-pimax100-bin0p1.fits" ))
#UCHUU['ALL_11.25'] = Table.read(os.path.join(uchuu_dir, "UCHUU_VLIM_ANY_11.25_Mstar_12.0_0.05_z_0.35_0err_glist-wprp-pimax100-bin0p1.fits" ))
#UCHUU['ALL_11.5'] = Table.read(os.path.join(uchuu_dir, "UCHUU_VLIM_ANY_11.5_Mstar_12.0_0.05_z_0.35_0err_glist-wprp-pimax100-bin0p1.fits" ))
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
CxG["S2_RS_11.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S2_CROSS_LS10_VLIM_RS_11.0_Mstar_12.0_0.05_z_0.35_N_1619838-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

CxG["S1_BC_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_BC_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits" ) )
CxG["S1_RS_10.5"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_RS_10.5_Mstar_12.0_0.05_z_0.26_N_3263228-wprp-pimax100-bin0p1-HpxMask.fits" ) )

CxG["S0_BC_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_BC_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits" ) )
CxG["S0_RS_10.0"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_RS_10.0_Mstar_12.0_0.05_z_0.18_N_2759238-wprp-pimax100-bin0p1-HpxMask.fits" ) )



CxG["S1_BC_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_BC_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
CxG["S1_RS_10.75"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S1_CROSS_LS10_VLIM_RS_10.75_Mstar_12.0_0.05_z_0.31_N_2802710-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )

CxG["S0_BC_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_BC_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )
CxG["S0_RS_10.25"]   = Table.read( os.path.join(LS10_bgs_vlim_dir, "eRASS1_VLIM_CLUSTERS_S0_CROSS_LS10_VLIM_RS_10.25_Mstar_12.0_0.05_z_0.22_N_3308841-wprp-pimax100-bin0p1-HpxMask_JK.fits" ) )



#
#
# S1 10.75
#
#


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS-S1-clusters.png')
plt.figure(13, (6,6))
t_wp = BGS["ANY_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed',  label=r'G1075', color='grey')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='grey', alpha=0.4)
t_wp = BGS["RS_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed', color='darkred')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkred', alpha=0.4)
t_wp = BGS["BC_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed', color='darkblue')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkblue', alpha=0.4)

#t_wp = CLU["S1_0.05_z_0.31"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=1,  ls='solid', label=r'Clusters, $\log_{10}(L_X\; [\mathrm{erg/s}])>43.1$', color='orange')
t_wp = CxG["S1_ANY_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='C1xG1075', color='grey')
t_wp = CxG["S1_BC_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='blue-cloud', color='darkblue')
t_wp = CxG["S1_RS_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='red-sequence', color='darkred')

plt.ylim((1, 2e3))
plt.xlim((0.03, 30))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$")#")
plt.legend(loc=4, fontsize=12,ncol=1, title=r'C1xG1075')#, title='LS10, r<19.5, 0.05<z<0.22')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)




p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS-S1-clusters-ratio.png')
plt.figure(13, (6,6))
t_ref = BGS["ANY_10.75"]
t_wp = BGS["RS_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='RS, Gal x Gal', color='darkred')
t_wp = BGS["BC_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='BC, Gal x Gal ', color='darkblue')

t_ref = CxG["S1_ANY_10.75"]
t_wp = CxG["S1_BC_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S1_RS_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')

plt.axhline(1,color='k', ls='dotted')
plt.ylim((0, 2))
plt.xlim((0.03, 30))
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w^{\rm RS\; or\; BC}_p(r_p)/w^{\rm All}_p$")
plt.legend(loc=1, fontsize=10, title='C1xG1075')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'wprp-obs-M1075-ANY-BC-RS-S1-clusters-ratioCROSSAUTO.png')
plt.figure(13, (6,6))
t_ref = BGS["ANY_10.75"]
t_wp = CxG["S1_ANY_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', color='grey', label='all galaxies')

t_ref = BGS["RS_10.75"]
t_wp = CxG["S1_RS_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', color='darkred', label='red-sequence')

t_ref = BGS["BC_10.75"]
t_wp = CxG["S1_BC_10.75"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', color='darkblue', label='blue-cloud')
b_gal = 1.43
b_clu = 3.35
b_clu_e = 0.23
xx = t_wp['rp_mid']
plt.fill_between(xx, y1=np.ones_like(xx)*( ((b_clu-b_clu_e)/b_gal)), y2=np.ones_like(xx)*( ((b_clu+b_clu_e)/b_gal)), color='darkgreen', alpha=0.2)
plt.axhline((b_clu/b_gal), color='darkgreen', ls='dashed', label='$b_{C}/b_{G}$')

plt.ylim((0.9, 30))
plt.xlim((0.03, 30))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w^{C1xG1075}_p/w^{G1075}_p$")
plt.legend(loc=3, fontsize=12,ncol=2, title='C1xG1075')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


#
#
# S0 10.25
#
#


# plot wprp
p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS-S0-clusters.png')
plt.figure(13, (6,6))
t_wp = BGS["ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed',  label=r'G1025', color='grey')
f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='grey', alpha=0.4)

t_wp = BGS["RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed', color='darkred')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkred', alpha=0.4)

t_wp = BGS["BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=2,  ls='dashed', color='darkblue')
#f_err = t_wp['wprp_JK_std']/t_wp['wprp_JK_mean']
#plt.fill_between(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid']*(1-f_err), t_wp['wprp']*t_wp['rp_mid']*(1+f_err), color='darkblue', alpha=0.4)

#t_wp = CLU["S0_0.05_z_0.22"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=1,  ls='solid', label='C1', color='orange')
t_wp = CxG["S0_ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='C1xG1025', color='grey')
t_wp = CxG["S0_BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='blue-cloud', color='darkblue')
t_wp = CxG["S0_RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']*t_wp['rp_mid'], lw=3,  ls='solid', label='red-sequence', color='darkred')

plt.ylim((1, 2e3))
plt.xlim((0.03, 30))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$r_p\times w_p(r_p)$")#")
plt.legend(loc=4, fontsize=12,ncol=1, title=r'C1xG1025')
plt.tight_layout()
plt.savefig(p2_fig)
plt.clf()
print(p2_fig)


p2_fig = os.path.join( fig_dir, 'wprp-obs-M1025-ANY-BC-RS-S0-clusters-ratio.png')
plt.figure(13, (6,6))
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
plt.xlim((0.03, 30))
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
plt.figure(13, (6,6))
t_ref = BGS["ANY_10.25"]
t_wp = CxG["S0_ANY_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', color='grey', label='all galaxies')

t_ref = BGS["RS_10.25"]
t_wp = CxG["S0_RS_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid',  color='darkred', label='red-sequence')

t_ref = BGS["BC_10.25"]
t_wp = CxG["S0_BC_10.25"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid',  color='darkblue', label='blue-cloud')
b_gal = 1.29
b_clu = 2.95
b_clu_e = 0.21
xx = t_wp['rp_mid']
plt.fill_between(xx, y1=np.ones_like(xx)*( ((b_clu-b_clu_e)/b_gal)), y2=np.ones_like(xx)*( ((b_clu+b_clu_e)/b_gal)), color='darkgreen', alpha=0.2)
plt.axhline((b_clu/b_gal), color='darkgreen', ls='dashed', label='$b_{C}/b_{G}$')

plt.ylim((0.9, 30))
plt.xlim((0.03, 30))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [Mpc/h]")
plt.ylabel(r"$w^{C0xG1025}_p/w^{G1025}_p$")
plt.legend(loc=3, fontsize=12,ncol=2, title=r'C0xG1025')
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
plt.figure(13, (6,6))
t_wp = BGS["ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='11.0<M*<12, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')

t_wp = CLU["S2_0.05_z_0.35"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='Clusters, S2 LX>43.3, N='+str(int(t_wp['N_data'][0])), color='orange')
t_wp = CxG["S2_ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='Galaxies x Clusters ', color='grey')
t_wp = CxG["S2_BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S2_RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')
plt.ylim((0.1, 3e4))
plt.xlim((0.03, 30))
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
plt.figure(13, (6,6))
t_ref = BGS["ANY_11.0"]
#plt.plot(t_ref['rp_mid'], t_ref['wprp'], lw=4,  ls='solid', label='11.0<M*<12, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', label='RS, Gal x Gal', color='darkred')
t_wp = BGS["BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=3,  ls='solid', label='BC, Gal x Gal', color='darkblue')

#t_wp = CLU["S2_0.05_z_0.35"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='Clusters, S2 LX>43.3, N='+str(int(t_wp['N_data'][0])), color='orange')
t_ref = CxG["S2_ANY_11.0"]
#plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='Galaxies x Clusters ', color='grey')
t_wp = CxG["S2_BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S2_RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')
plt.axhline(1,color='k', ls='dotted')
plt.ylim((0, 2))
plt.xlim((0.03, 30))
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
plt.figure(13, (6,6))
t_ref = BGS["ANY_11.0"]
t_wp = CxG["S2_ANY_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='Any', color='grey')

t_ref = BGS["RS_11.0"]
t_wp = CxG["S2_RS_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='RS', color='darkred')

t_ref = BGS["BC_11.0"]
t_wp = CxG["S2_BC_11.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='BC ', color='darkblue')

plt.ylim((0.1, 100))
plt.xlim((0.03, 30))
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
plt.figure(13, (6,6))
t_wp = BGS["ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='10.5<M*<12, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')

t_wp = CLU["S1_0.05_z_0.35"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='Clusters, S1 LX>43.1, N='+str(int(t_wp['N_data'][0])), color='orange')
t_wp = CxG["S1_ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='Galaxies x Clusters ', color='grey')
t_wp = CxG["S1_BC_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S1_RS_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')

plt.ylim((0.1, 3e4))
plt.xlim((0.03, 30))
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
plt.figure(13, (6,6))
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
plt.xlim((0.03, 30))
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
plt.figure(13, (6,6))
t_ref = BGS["ANY_10.5"]
t_wp = CxG["S1_ANY_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='Any', color='grey')

t_ref = BGS["RS_10.5"]
t_wp = CxG["S1_RS_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='RS', color='darkred')

t_ref = BGS["BC_10.5"]
t_wp = CxG["S1_BC_10.5"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='BC ', color='darkblue')

plt.ylim((0.1, 100))
plt.xlim((0.03, 30))
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
plt.figure(13, (6,6))
t_wp = BGS["ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=4,  ls='solid', label='10.0<M*<12, N='+str(int(t_wp['N_data'][0])), color='grey')
t_wp = BGS["RS_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='red-sequence, N='+str(int(t_wp['N_data'][0])), color='darkred')
t_wp = BGS["BC_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='blue-cloud, N='+str(int(t_wp['N_data'][0])), color='darkblue')

t_wp = CLU["S0_0.05_z_0.35"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=3,  ls='solid', label='Clusters, S0 LX>42.7, N='+str(int(t_wp['N_data'][0])), color='orange')
t_wp = CxG["S0_ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='Galaxies x Clusters ', color='grey')
t_wp = CxG["S0_BC_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='BC Gal x Clu ', color='darkblue')
t_wp = CxG["S0_RS_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp'], lw=2,  ls='dashed', label='RS Gal x Clu ', color='darkred')

plt.ylim((0.1, 3e4))
plt.xlim((0.03, 30))
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
plt.figure(13, (6,6))
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
plt.xlim((0.03, 30))
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
plt.figure(13, (6,6))
t_ref = BGS["ANY_10.0"]
t_wp = CxG["S0_ANY_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='Any', color='grey')

t_ref = BGS["RS_10.0"]
t_wp = CxG["S0_RS_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='RS', color='darkred')

t_ref = BGS["BC_10.0"]
t_wp = CxG["S0_BC_10.0"]
plt.plot(t_wp['rp_mid'], t_wp['wprp']/t_ref['wprp'][::2], lw=3,  ls='solid', label='BC ', color='darkblue')

plt.ylim((0.1, 100))
plt.xlim((0.03, 30))
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



