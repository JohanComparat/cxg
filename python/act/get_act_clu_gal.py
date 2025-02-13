import glob
import os, sys
import numpy as np

from sklearn.neighbors import BallTree
from astropy.table import Table, Column, vstack, hstack
nl = lambda selection : len(selection.nonzero()[0])
deg_to_rad = np.pi/180.

from astropy.cosmology import FlatLambdaCDM, Planck18
import astropy.units as u
cosmo = Planck18

all_files = np.array(glob.glob("/sps/euclid/OU-LE3/VMPZ-ID/REGREPROC1_R2/CatRed_CL/EUC_MER_PHZ_CA*fits.gz"))

# take the cluster catalog and get all catalog entries within 2 R lambda, ave all matches with the memmatchid set.
ACT = Table.read('/pbs/home/j/jcomparat/ACT_DR5_cluster-catalog_v1.1_eromapper_catalog.fit')
#ACT = Table.read('/media/sf_Shared/data/eromapper/data/cal/legacy_dr10_grz_z_v0.3_run/ACT_DR5_cluster-catalog_v1.1_eromapper_catalog.fit.gz')

ACT = ACT[( ACT['z_lambda']>0 ) & ( ACT['r_lambda'] > 0 )]

coord_S1 = deg_to_rad *np.transpose([ACT['dec'], ACT['ra']])
Tree_S1 = BallTree(coord_S1, metric='haversine')
ACT['z_lambda']
r_rad = ((ACT['r_lambda']*u.Mpc).to(u.kpc)/cosmo.kpc_proper_per_arcmin(ACT['z_lambda'])).to(u.rad)
ID_ACT = np.arange(len(ACT))
ACT['line_ID'] = ID_ACT
D_MAX = 0.4 # arc seconds

#distances_out_i, ids_out_i = Tree_S1.query_radius(coord_S1, r=r_rad*3, return_distance = True)

jj=0
for fname in all_files[::-1]:
    #fname = '/media/sf_Shared/data/Euclid/EUC_MER_PHZ_CAT_tresh0.8__TILE102165306-9535B5_20240718T180853.884342Z_00.00.fits.gz'
    print(fname, jj)
    jj+=1
    if os.path.isfile(fname):
        ff = Table.read(fname)
        coord_S2 = deg_to_rad *np.transpose([ff['DECLINATION'], ff['RIGHT_ASCENSION']])
        if len(coord_S2)>2:
            Tree_S2 = BallTree(coord_S2, metric='haversine')
            ids_out_i, distances_out_i = Tree_S2.query_radius(coord_S1, r=r_rad*3, return_distance = True)
            if len(np.hstack((ids_out_i)) )>=2:
                act_ids = []
                for jjj, el in zip(ACT['line_ID'], ids_out_i):
                    if len(el)>=1:
                        act_ids.append( ( np.ones_like(el)*jjj).astype('int') )

                act_ids = np.hstack(( act_ids ))
                distances_out = np.hstack(( distances_out_i ))
                ids_out = np.hstack(( ids_out_i )).astype('int')
                if len(ids_out)>0:
                    print('act_ids', act_ids)
                    print('ids_out', ids_out)
                    print('distances_out', distances_out)
                    t_out = ff[ids_out]
                    t_out['line_ID'] = act_ids
                    t_out['distance_cluster_rad'] = distances_out
                    print('ACT has', len(t_out), 'galaxies matching')
                    out_file = '/pbs/home/j/jcomparat/ACT_DR5_CLU_AND_' + os.path.basename(fname)
                    t_out.write(out_file, overwrite = True)
                    print(out_file, 'written')
