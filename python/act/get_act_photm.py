

import glob
import os, sys
import numpy as np

from sklearn.neighbors import BallTree
from astropy.table import Table, Column, vstack, hstack
nl = lambda selection : len(selection.nonzero()[0])
deg_to_rad = np.pi/180.


all_files = np.array(glob.glob("/sps/euclid/OU-LE3/VMPZ-ID/REGREPROC1_R2/CatRed_CL/EUC_MER_PHZ_CA*fits.gz"))
ACT = Table.read('/pbs/home/j/jcomparat/act-mem-euclid-test.fits.gz')
#ACT = Table.read('/media/sf_Shared/data/Euclid/act-mem-euclid-test.fits.gz')
ID_ACT = np.arange(len(ACT))
D_MAX = 0.4 # arc seconds
coord_S1 = deg_to_rad *np.transpose([ACT['dec'], ACT['ra']])
Tree_S1 = BallTree(coord_S1, metric='haversine')
jj=0
for fname in all_files[::-1][220:]:
    #fname = '/media/sf_Shared/data/Euclid/EUC_MER_PHZ_CAT_tresh0.8__TILE102165306-9535B5_20240718T180853.884342Z_00.00.fits.gz'
    print(fname, jj)
    jj+=1
    if os.path.isfile(fname):
        ff = Table.read(fname)
        coord_S2 = deg_to_rad *np.transpose([ff['DECLINATION'], ff['RIGHT_ASCENSION']])
        if len(coord_S2)>2:
            Tree_S2 = BallTree(coord_S2, metric='haversine')
            distances_out_i, ids_out_i = Tree_S2.query(coord_S1, k=1, return_distance = True)
            distances_out = np.transpose(distances_out_i)[0]
            ids_out = np.transpose(ids_out_i)[0]
            matching = np.hstack((distances_out))<D_MAX/3600*deg_to_rad
            print('ACT has', len(ids_out[matching]))
            if len(ids_out[matching])>0:
                t_out = hstack(( ACT[ID_ACT[matching]], ff[ids_out[matching]] ))
                out_file = '/pbs/home/j/jcomparat/act-mem-euclid-test_AND_' + os.path.basename(fname)
                t_out.write(out_file, overwrite = True)
                print(out_file, 'written')
