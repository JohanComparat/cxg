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
markersize = 20

fig_dir  ='../figures/CVCM'

d1025 = '../data/_ANY_10.25_Mstar_12.0_0.05_z_0.22_N_3308841_JK_wprp100'
d1075 = '../data/_ANY_10.75_Mstar_12.0_0.05_z_0.31_N_2802710_JK_wprp100'
dC0xG1025 = '../data/C0xG1025'
dC1xG1075 = '../data/C1xG1075'

pcf_list_1025_ANY_NSIDE_04 = np.array( glob.glob( os.path.join( dC0xG1025, '*_ANY_*NSIDE_04*.fits' )))
pcf_list_1025_ANY_NSIDE_08 = np.array( glob.glob( os.path.join( dC0xG1025, '*_ANY_*NSIDE_08*.fits' )))
pcf_list_1025_BC__NSIDE_04 = np.array( glob.glob( os.path.join( dC0xG1025,  '*_BC_*NSIDE_04*.fits' )))
pcf_list_1025_BC__NSIDE_08 = np.array( glob.glob( os.path.join( dC0xG1025,  '*_BC_*NSIDE_08*.fits' )))
pcf_list_1025_RS__NSIDE_04 = np.array( glob.glob( os.path.join( dC0xG1025,  '*_RS_*NSIDE_04*.fits' )))
pcf_list_1025_RS__NSIDE_08 = np.array( glob.glob( os.path.join( dC0xG1025,  '*_RS_*NSIDE_08*.fits' )))

pcf_list_1075_ANY_NSIDE_04 = np.array( glob.glob( os.path.join( dC1xG1075, '*_ANY_*NSIDE_04*.fits' )))
pcf_list_1075_ANY_NSIDE_08 = np.array( glob.glob( os.path.join( dC1xG1075, '*_ANY_*NSIDE_08*.fits' )))
pcf_list_1075_BC__NSIDE_04 = np.array( glob.glob( os.path.join( dC1xG1075,  '*_BC_*NSIDE_04*.fits' )))
pcf_list_1075_BC__NSIDE_08 = np.array( glob.glob( os.path.join( dC1xG1075,  '*_BC_*NSIDE_08*.fits' )))
pcf_list_1075_RS__NSIDE_04 = np.array( glob.glob( os.path.join( dC1xG1075,  '*_RS_*NSIDE_04*.fits' )))
pcf_list_1075_RS__NSIDE_08 = np.array( glob.glob( os.path.join( dC1xG1075,  '*_RS_*NSIDE_08*.fits' )))

basename = 'C0_1025_ANY_NSIDE_04'
p_2_2PCF_all = pcf_list_1025_ANY_NSIDE_04

def make_plots(p_2_2PCF_all, basename, title_name):
	print(basename, len(p_2_2PCF_all))
	mat = []
	for p_2_2PCF_all_i in p_2_2PCF_all: #[0]
		t = Table.read(p_2_2PCF_all_i)
		mat.append(t['wprp'])

	mat = np.transpose(mat)
	print(mat.shape)
	cv = np.cov(mat)
	cc = np.corrcoef(mat)

	rr = np.outer(t['rp_mid'], t['rp_mid'])
	xx = np.array([ t['rp_mid'] for jj in np.arange(len(t['rp_mid'])) ])
	yy = np.array([ np.ones(len(t['rp_mid']))*t['rp_mid'][jj] for jj in np.arange(len(t['rp_mid'])) ])
	x0 = np.hstack(( xx ))
	y0 = np.hstack(( yy ))
	z0 = np.hstack(( cc ))
	z1 = np.hstack(( cv ))

	p_2_2PCF_figure = os.path.join(fig_dir, basename +'_CC_wprp-pimax100.png')
	plt.figure(3, (6.5, 5))
	plt.scatter(x0, y0, c=z0, s=markersize, marker='s', vmin=-1, vmax=1, cmap='RdYlGn')
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$r_p$ [$h^{-1}$Mpc]',fontsize=14)
	plt.ylabel(r'$r_p$ [$h^{-1}$Mpc]  ',fontsize=14)
	plt.colorbar(label='correlation coefficient, N JK='+str(mat.shape[1]))
	#plt.title(title_name)
	plt.tight_layout()
	plt.savefig(p_2_2PCF_figure)
	plt.clf()
	print(p_2_2PCF_figure)


	p_2_2PCF_figure = os.path.join(fig_dir, basename +'_pCV_wprp-pimax100.png')
	plt.figure(3, (6.5, 5))
	plt.scatter(x0, y0, c=np.log10(z1), s=markersize, marker='s', vmin=-4, vmax=4., cmap='RdYlGn')
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.xscale('log')
	plt.yscale('log')
	#plt.ylim((1e-5, 2000))
	#plt.xlim((8e-4, 3))
	#plt.legend(loc='upper right', ncol=3, fontsize=12, title = str_params)
	plt.xlabel(r'$r_p$ [$h^{-1}$Mpc]',fontsize=14)
	plt.ylabel(r'$r_p$ [$h^{-1}$Mpc]  ',fontsize=14)
	plt.colorbar(label='log10 covariance matrix, N JK='+str(mat.shape[1]))
	#plt.title(title_name)
	plt.tight_layout()
	plt.savefig(p_2_2PCF_figure)
	plt.clf()
	print(p_2_2PCF_figure)

	p_2_2PCF_figure = os.path.join(fig_dir, basename +'_mCV_wprp-pimax100.png')
	plt.figure(3, (6.5, 5))
	plt.scatter(x0, y0, c=np.log10(-z1), s=markersize, marker='s', vmin=-4, vmax=4., cmap='RdYlGn')
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.xscale('log')
	plt.yscale('log')
	#plt.ylim((1e-5, 2000))
	#plt.xlim((8e-4, 3))
	#plt.legend(loc='upper right', ncol=3, fontsize=12, title = str_params)
	plt.xlabel(r'$r_p$ [$h^{-1}$Mpc]',fontsize=14)
	plt.ylabel(r'$r_p$ [$h^{-1}$Mpc]  ',fontsize=14)
	plt.colorbar(label='- log10 covariance matrix, N JK='+str(mat.shape[1]))
	#plt.title(title_name)
	plt.tight_layout()
	plt.savefig(p_2_2PCF_figure)
	plt.clf()
	print(p_2_2PCF_figure)


	p_2_2PCF_figure = os.path.join(fig_dir, basename +'_diagRelErr-pimax100.png')
	plt.figure(3, (6.5, 5))
	plt.plot(t['rp_mid'], 100*abs(mat.std(axis=1)/mat.mean(axis=1)))
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim((0.1, 20))
	#plt.xlim((8e-4, 3))
	#plt.legend(loc='upper right', ncol=3, fontsize=12, title = str_params)
	plt.xlabel(r'$r_p$ [$h^{-1}$Mpc]',fontsize=14)
	plt.ylabel(r'$\Delta w_p(r_p)$ [%]',fontsize=14)
	#plt.colorbar(label='- log10 covariance matrix, N JK='+str(mat.shape[1]))
	#plt.title(title_name)
	plt.tight_layout()
	plt.savefig(p_2_2PCF_figure)
	plt.clf()
	print(p_2_2PCF_figure)


make_plots(pcf_list_1025_ANY_NSIDE_04, 'C0_1025_ANY_NSIDE_04', 'C0 x G1025')
make_plots(pcf_list_1025_ANY_NSIDE_08, 'C0_1025_ANY_NSIDE_08', 'C0 x G1025')
make_plots(pcf_list_1025_BC__NSIDE_04, 'C0_1025_BC__NSIDE_04', 'C0 x G1025 BC')
make_plots(pcf_list_1025_BC__NSIDE_08, 'C0_1025_BC__NSIDE_08', 'C0 x G1025 BC')
make_plots(pcf_list_1025_RS__NSIDE_04, 'C0_1025_RS__NSIDE_04', 'C0 x G1025 RS')
make_plots(pcf_list_1025_RS__NSIDE_08, 'C0_1025_RS__NSIDE_08', 'C0 x G1025 RS')

make_plots(pcf_list_1075_ANY_NSIDE_04, 'C1_1075_ANY_NSIDE_04', 'C1 x 1075')
make_plots(pcf_list_1075_ANY_NSIDE_08, 'C1_1075_ANY_NSIDE_08', 'C1 x 1075')
make_plots(pcf_list_1075_BC__NSIDE_04, 'C1_1075_BC__NSIDE_04', 'C1 x 1075 BC')
make_plots(pcf_list_1075_BC__NSIDE_08, 'C1_1075_BC__NSIDE_08', 'C1 x 1075 BC')
make_plots(pcf_list_1075_RS__NSIDE_04, 'C1_1075_RS__NSIDE_04', 'C1 x 1075 RS')
make_plots(pcf_list_1075_RS__NSIDE_08, 'C1_1075_RS__NSIDE_08', 'C1 x 1075 RS')



