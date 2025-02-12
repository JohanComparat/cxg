dir_clusters='/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples'
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S3.fit         
randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits
randoms-1-0-erass1sky-hod-cutselfunc20230731_S4.fits
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S4.fit
randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits
erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit
randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits
randoms-1-0-erass1sky-hod-cutselfunc20230731_S3.fits

# stack euclid galaxy numbers as a function of separation
# in angular space
# around each cluster, retrieve the histogram of pair count of galaxies as a function of color
# around random points, retrieve the histogram of pair count of galaxies as a function of color
# up to 10'

# start with LS10
# color g-r

# S0
python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_010_z_011.npy" \
'g_mag' 'r_mag' 0.1 0.11

python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_011_z_012.npy" \
'g_mag' 'r_mag' 0.11 0.12

python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_012_z_013.npy" \
'g_mag' 'r_mag' 0.12 0.13

python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_013_z_014.npy" \
'g_mag' 'r_mag' 0.13 0.14

python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_014_z_015.npy" \
'g_mag' 'r_mag' 0.14 0.15


python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_015_z_016.npy" \
'g_mag' 'r_mag' 0.15 0.16


python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_016_z_017.npy" \
'g_mag' 'r_mag' 0.16 0.17


python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_017_z_018.npy" \
'g_mag' 'r_mag' 0.17 0.18


python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_018_z_019.npy" \
'g_mag' 'r_mag' 0.18 0.19

python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_DATA_S0_RAND_gr_019_z_020.npy" \
'g_mag' 'r_mag' 0.19 0.20


# S1 0.1-0.3
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"Counts_DATA_S1_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" \
"g_mag" "r_mag" """+str(z_lo)+' '+str(z_hi)

z_alls = np.arange(0.1, 0.3, 0.01)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.01,2))
    print(cmd)
    os.system(cmd)

# S2 0.1-0.4
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compute_wprp_Files_MaskHpx_MstarDEP_cross.py \
"/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" \
"MergeALL_BGSlike_LPH.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"Counts_DATA_S2_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" \
"g_mag" "r_mag" """+str(z_lo)+' '+str(z_hi)

z_alls = np.arange(0.1, 0.4, 0.01)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.01,2))
    print(cmd)
    os.system(cmd)

#
# MAKE A NICE FIGURE
#

#
# Do the same with Euclid galaxies (deep !)
#

python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_010_z_011.npy'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_011_z_012.npy'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_012_z_013.npy'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_013_z_014.npy'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_014_z_015.npy'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_015_z_016.npy'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_016_z_017.npy'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_017_z_018.npy'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_018_z_019.npy'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_019_z_020.npy'



# S1 0.1-0.3
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compare_ct.py "/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" "Counts_DATA_S1_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" """

z_alls = np.arange(0.1, 0.3, 0.01)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.01,2))
    print(cmd)
    os.system(cmd)

# S2 0.1-0.4
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compare_ct.py "/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep" "Counts_DATA_S2_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" """

z_alls = np.arange(0.1, 0.4, 0.01)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.01,2))
    print(cmd)
    os.system(cmd)



python plot_RS.py S0
python plot_RS.py S1
python plot_RS.py S2
