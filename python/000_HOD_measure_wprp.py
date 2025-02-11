import os, sys
import numpy as np
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
from astropy.table import Table
from astropy import constants as const
import pyccl as ccl
import matplotlib.pyplot as plt
from scipy import integrate
sys.path.append('/home/rseppi/eROclustering/src')
from utils import custom_plot

#p2data = '/home/rseppi/eROclustering/data/catalogs/eromapper_main_sample_v230303.fit'
p2data = '/home/rseppi/eROclustering/data/catalogs/erass1cl_main_v1.1_w_xrayresu_w_expbkg.fit'
p2rand = '/home/rseppi/eROclustering/data/catalogs/randoms-1-0-erass1sky-hod-cutselfunc20230620.fits'

zmin, zmax = 0.1, 0.8
pimax = 40.
h = 0.6776

data = Table.read(p2data)
sel = (data['Z_LAMBDA']>zmin) & (data['Z_LAMBDA']<zmax) & (data['PCONT']>=0) & (data['PCONT']<0.5) & (data['DEC']<32.375)
data = data[sel]
print('len data before cts cut',len(data))
sel = data['CTS500X']>20.
data = data[sel]
print('len data after cts cut',len(data))
print(np.min(data['CTS500X']))

randoms = Table.read(p2rand)
sel = (randoms['PSFDEPTH_G'] > 0) & (randoms['PSFDEPTH_R'] > 0) & (randoms['PSFDEPTH_Z'] > 0)
randoms = randoms[sel]
sel = (randoms['DEC'] < 32.375) & (randoms['redshift']>zmin) & (randoms['redshift']<zmax)
randoms = randoms[sel]
sel = (randoms['MASKBITS'] != 1) & (randoms['MASKBITS'] != 2) & (randoms['MASKBITS'] != 16) & \
      (randoms['MASKBITS'] != 128) & (randoms['MASKBITS'] != 1024) & (randoms['MASKBITS'] != 2048) & \
      (randoms['MASKBITS'] != 4096) & (randoms['MASKBITS'] != 8192)
randoms = randoms[sel]

edges1 = np.geomspace(0.5, 40., 10)
edges2 = np.geomspace(0.5, 40., 8)
edges3 = np.geomspace(0.5, 40., 6)
edges = [edges1, edges2, edges3]

def Wprp(data, randoms, p2outwp, edges):
    bins = (edges[1:] + edges[:-1]) / 2.
    print('redshift=',np.average(data['Z_LAMBDA']))
    print('Mass=',np.average(data['M500X'])*1e13*0.8)
    ra, dec = np.array(data['RA'], dtype=np.float64), np.array(data['DEC'], dtype=np.float64)
    rand_ra, rand_dec = np.array(randoms['RA'], dtype=np.float64), np.array(randoms['DEC'], dtype=np.float64)
    ra_plot = np.array([rand_ra, ra])
    dec_plot = np.array([rand_dec, dec])
    strin = p2outwp.split('.')[-2].split('exp')[-1]
    p2outsky = '/home/rseppi/eROclustering/figures/sky/eromapper_main_sample_v230303_wprp_z_'+str(zmin)+'_'+str(zmax)+'_exp'+strin+'.png'
    #if strin=='250':
    #    custom_plot.plot_sky(ra_arr=ra_plot, dec_arr=dec_plot, p2out=p2outsky)


    # Planck cosmology
    cosmo = ccl.Cosmology(Omega_c=0.26066, Omega_b=0.04897, h=h, sigma8=0.8228, n_s=0.96)
    dC = ccl.comoving_radial_distance(cosmo, 1. / (1. + np.array(data['Z_LAMBDA'], dtype=np.float64))) * h #Mpc/h
    dC_rand = ccl.comoving_radial_distance(cosmo, 1. / (1. + np.array(randoms['redshift'], dtype=np.float64))) * h

    N = len(ra)
    rand_N = len(rand_ra)
    print(N, 'sources')
    print(rand_N, 'randoms')

    nthreads = 2
    print('DD...')
    autocorr=1
    DD_counts = DDrppi_mocks(autocorr=autocorr, cosmology=2, nthreads=nthreads, binfile=edges,
                             pimax=pimax,
                             RA1=ra,
                             DEC1=dec,
                             CZ1=dC,
                             is_comoving_dist=True)
    print('DR...')
    autocorr = 0
    DR_counts = DDrppi_mocks(autocorr=autocorr, cosmology=2, nthreads=nthreads, pimax=pimax,
                             binfile=edges, RA1=ra, DEC1=dec, CZ1=dC,
                             RA2=rand_ra, DEC2=rand_dec, CZ2=dC_rand, is_comoving_dist=True)
    print('RR...')
    autocorr = 1
    RR_counts = DDrppi_mocks(autocorr=autocorr, cosmology=2, nthreads=nthreads, pimax=pimax,
                             binfile=edges, RA1=rand_ra, DEC1=rand_dec, CZ1=dC_rand, is_comoving_dist=True)
    N = len(ra)
    rand_N = len(rand_ra)
    nbins = len(edges) - 1
    print('Converting to wprp...')
    wprp = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
                                      DD_counts, DR_counts,
                                      DR_counts, RR_counts, nbins, pimax)

    np.savetxt(p2outwp, np.transpose([bins, wprp]))#, DD_counts['npairs'], DR_counts['npairs'], RR_counts['npairs']]), fmt='%s')
    print('done!', p2outwp)
    return wprp


exp = data['EXPOSURE']
cts = data['CTS500X']

#exp1, exp2, exp3, exp4 = 30, 120, 250
exp1, exp2, exp3 = 20, 150, 300
ctsmin = 20
selcts = cts>ctsmin
sel1 = (selcts)&(exp>exp1)
sel2 = (selcts)&(exp>exp2)
sel3 = (selcts)&(exp>exp3)
#sel4 = (selcts)&(exp>exp4)

exp = randoms['Texp']
sel1r = (exp>exp1)
sel2r = (exp>exp2)
sel3r = (exp>exp3)
#sel4r = (exp>exp4)

wprp0 = Wprp(data[sel1], randoms[sel1r], '/home/rseppi/eROclustering/data/2pcf/wprp/eRASS1_wprp_HOD_cts'+str(ctsmin)+'_exp'+str(exp1).zfill(3)+'.2pcf', edges1)
wprp1 = Wprp(data[sel2], randoms[sel2r], '/home/rseppi/eROclustering/data/2pcf/wprp/eRASS1_wprp_HOD_cts'+str(ctsmin)+'_exp'+str(exp2).zfill(3)+'.2pcf', edges2)
wprp2 = Wprp(data[sel3], randoms[sel3r], '/home/rseppi/eROclustering/data/2pcf/wprp/eRASS1_wprp_HOD_cts'+str(ctsmin)+'_exp'+str(exp3).zfill(3)+'.2pcf', edges3)
#wprp3 = Wprp(data[sel4], randoms[sel4r], '/home/rseppi/eROclustering/data/2pcf/wprp/eRASS1_wprp_HOD_cts'+str(ctsmin)+'_exp'+str(exp4).zfill(3)+'.2pcf')

wpl = [wprp0, wprp1, wprp2]#, wprp3]
expl = [exp1, exp2, exp3]#, exp4]

plt.figure(figsize=(9,7))
for jj,wp in enumerate(wpl):
    strin = '%d'%expl[jj]
    bins = (edges[jj][1:] + edges[jj][:-1]) / 2.
    plt.plot(bins, wp, label='> ' + strin + ' s', lw=4, c='C'+str(jj))
    plt.scatter(bins, wp, s=300, color='C'+str(jj))
    print(wp)
plt.ylabel(r'w$_p$(r$_p$)', fontsize=25)
plt.xlabel(r'r$_p$ [Mpc/h]', fontsize=25)
plt.yscale('log')
plt.xscale('log')
#plt.ylim(0.5)
plt.legend(fontsize=22)
plt.tick_params(labelsize=25, direction='in', which='both')
plt.tight_layout()
outf = '/home/rseppi/eROclustering/figures/HOD/wprp_eRASS1_model_HOD.png'
plt.savefig(outf)
plt.close()
print(outf)
sys.exit()


#Create 1 single catalog adding x-ray properties
p2main = '/home/rseppi/eROclustering/data/catalogs/eromapper_main_sample_v230303.fit'
p2resu = '/home/rseppi/eROclustering/data/catalogs/erass1_main_xrayresu_20230315.fits'

print('Reading...')
main = Table.read(p2main)
resu = Table.read(p2resu)

print('Cross matching...')
arr, index_main, index_resu = np.intersect1d(np.array(main['RA']),np.array(resu['RA']), return_indices=True)

newcols = np.array(resu.colnames)[~np.in1d(resu.colnames, main.colnames)]

for coln in newcols:
    if len(resu[coln].shape)>1:
        qty = np.ones((len(main),512)) * -99.
    else:
        qty = np.ones(len(main)) * -99.
    qty[index_main] = resu[coln]
    main[coln] = qty

p2out_main = '/home/rseppi/eROclustering/data/catalogs/eromapper_main_sample_v230303_w_xrayresu.fit'
main.write(p2out_main)