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

z_alls = np.arange(0.1, 0.39, 0.01)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.01,2))
    print(cmd)
    os.system(cmd)



python plot_RS.py S0
python plot_RS.py S1
python plot_RS.py S2



rsync -avz jcomparat@cca.in2p3.fr:~/color_catalog.fits /home/comparat/sf_Shared/data/Euclid/data/



rsync -avz jcomparat@cca.in2p3.fr:/sps/euclid/OU-LE3/CL/ial_workspace/Q1/data_Nov24/CatRed_filtered_noMagCut/CONCAT_CATRED_*.fits  /home/comparat/sf_Shared/data/Euclid/data/

rsync -avz jcomparat@cca.in2p3.fr:/sps/euclid/Users/sunayana/Q1_Nov24/no_mag/????_no_mag/galaxy_cats/CATRED_filtered_*.fits /home/comparat/sf_Shared/data/Euclid/data/


# S0 0.1-0.2
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compute_wprp_Euclid_cross.py \
"/home/comparat/sf_Shared/data/Euclid/data" \
"color_catalog.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_Euclid_S0_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" \
"gr" """+str(z_lo)+' '+str(z_hi)

z_alls = np.arange(0.1, 0.2, 0.02)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.02,2))
    print(cmd)

# S1 0.1-0.3
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compute_wprp_Euclid_cross.py \
"/home/comparat/sf_Shared/data/Euclid/data" \
"color_catalog.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"Counts_Euclid_S1_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" \
"gr" """+str(z_lo)+' '+str(z_hi)

z_alls = np.arange(0.1, 0.3, 0.02)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.02,2))
    print(cmd)

# S2 0.1-0.4
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compute_wprp_Euclid_cross.py \
"/home/comparat/sf_Shared/data/Euclid/data" \
"color_catalog.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"Counts_Euclid_S2_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" \
"gr" """+str(z_lo)+' '+str(z_hi)

z_alls = np.arange(0.1, 0.4, 0.02)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.02,2))
    print(cmd)

python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" "Counts_Euclid_S0_RAND_gr_10_z_12.npy" "gr" 0.1 0.12
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" "Counts_Euclid_S0_RAND_gr_12_z_14.npy" "gr" 0.12 0.14
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" "Counts_Euclid_S0_RAND_gr_14_z_16.npy" "gr" 0.14 0.16
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" "Counts_Euclid_S0_RAND_gr_16_z_18.npy" "gr" 0.16 0.18
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" "Counts_Euclid_S0_RAND_gr_18_z_20.npy" "gr" 0.18 0.2
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_10_z_12.npy" "gr" 0.1 0.12
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_12_z_14.npy" "gr" 0.12 0.14
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_14_z_16.npy" "gr" 0.14 0.16
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_16_z_18.npy" "gr" 0.16 0.18
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_18_z_20.npy" "gr" 0.18 0.2
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_20_z_22.npy" "gr" 0.2 0.22
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_22_z_24.npy" "gr" 0.22 0.24
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_24_z_26.npy" "gr" 0.24 0.26
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_26_z_28.npy" "gr" 0.26 0.28
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_gr_28_z_30.npy" "gr" 0.28 0.3
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_10_z_12.npy" "gr" 0.1 0.12
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_12_z_14.npy" "gr" 0.12 0.14
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_14_z_16.npy" "gr" 0.14 0.16
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_16_z_18.npy" "gr" 0.16 0.18
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_18_z_20.npy" "gr" 0.18 0.2
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_20_z_22.npy" "gr" 0.2 0.22
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_22_z_24.npy" "gr" 0.22 0.24
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_24_z_26.npy" "gr" 0.24 0.26
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_26_z_28.npy" "gr" 0.26 0.28
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_28_z_30.npy" "gr" 0.28 0.3
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_30_z_32.npy" "gr" 0.3 0.32
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_32_z_34.npy" "gr" 0.32 0.34
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_34_z_36.npy" "gr" 0.34 0.36
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_36_z_38.npy" "gr" 0.36 0.38
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_38_z_40.npy" "gr" 0.38 0.4
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_gr_40_z_42.npy" "gr" 0.4 0.42


# S0 0.1-0.2
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compare_ct.py "/home/comparat/sf_Shared/data/Euclid/data" "Counts_Euclid_S0_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" """

z_alls = np.arange(0.1, 0.2, 0.02)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.02,2))
    print(cmd)
    os.system(cmd)

# S1 0.1-0.3
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compare_ct.py "/home/comparat/sf_Shared/data/Euclid/data" "Counts_Euclid_S1_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" """

z_alls = np.arange(0.1, 0.3, 0.02)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.02,2))
    print(cmd)
    os.system(cmd)

# S2 0.1-0.4
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compare_ct.py "/home/comparat/sf_Shared/data/Euclid/data" "Counts_Euclid_S2_RAND_gr_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" """

z_alls = np.arange(0.1, 0.4, 0.02)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.02,2))
    print(cmd)
    os.system(cmd)




python plot_RS_Euclid.py S0
python plot_RS_Euclid.py S1
python plot_RS_Euclid.py S2


# VIS - y

# S0 0.1-0.2
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compute_wprp_Euclid_cross.py \
"/home/comparat/sf_Shared/data/Euclid/data" \
"color_catalog.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" \
"Counts_Euclid_S0_RAND_vy_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" \
"vy" """+str(z_lo)+' '+str(z_hi)

z_alls = np.arange(0.1, 0.2, 0.05)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.05,2))
    print(cmd)

# S1 0.1-0.3
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compute_wprp_Euclid_cross.py \
"/home/comparat/sf_Shared/data/Euclid/data" \
"color_catalog.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" \
"Counts_Euclid_S1_RAND_vy_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" \
"vy" """+str(z_lo)+' '+str(z_hi)

z_alls = np.arange(0.1, 0.3, 0.05)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.05,2))
    print(cmd)

# S2 0.1-0.4
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compute_wprp_Euclid_cross.py \
"/home/comparat/sf_Shared/data/Euclid/data" \
"color_catalog.fits" \
"/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" \
"erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" \
"randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" \
"Counts_Euclid_S2_RAND_vy_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" \
"vy" """+str(z_lo)+' '+str(z_hi)

z_alls = np.arange(0.1, 0.4, 0.05)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.05,2))
    print(cmd)

python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" "Counts_Euclid_S0_RAND_vy_10_z_15.npy" "vy" 0.1 0.15
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S0.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S0.fits" "Counts_Euclid_S0_RAND_vy_15_z_20.npy" "vy" 0.15 0.2
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_vy_10_z_15.npy" "vy" 0.1 0.15
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_vy_15_z_20.npy" "vy" 0.15 0.2
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_vy_20_z_25.npy" "vy" 0.2 0.25
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S1.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S1.fits" "Counts_Euclid_S1_RAND_vy_25_z_30.npy" "vy" 0.25 0.3
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_vy_10_z_15.npy" "vy" 0.1 0.15
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_vy_15_z_20.npy" "vy" 0.15 0.2
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_vy_20_z_25.npy" "vy" 0.2 0.25
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_vy_25_z_30.npy" "vy" 0.25 0.3
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_vy_30_z_35.npy" "vy" 0.3 0.35
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_vy_35_z_40.npy" "vy" 0.35 0.4
python compute_wprp_Euclid_cross.py "/home/comparat/sf_Shared/data/Euclid/data" "color_catalog.fits" "/home/comparat/sf_Shared/data/erosita/observations/eRASS/cluster_clustering/eRASS1_CLU_VolLimSamples" "erass1cl_main_v2.0_w_xrayresu_w_expbkg_S2.fit" "randoms-1-0-erass1sky-hod-cutselfunc20230731_S2.fits" "Counts_Euclid_S2_RAND_vy_40_z_45.npy" "vy" 0.4 0.45


# S0 0.1-0.2
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compare_ct_vy.py "/home/comparat/sf_Shared/data/Euclid/data" "vy" "Counts_Euclid_S0_RAND_vy_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" """

z_alls = np.arange(0.1, 0.2, 0.05)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.05,2))
    print(cmd)
    os.system(cmd)

# S1 0.1-0.3
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compare_ct_vy.py "/home/comparat/sf_Shared/data/Euclid/data" "vy" "Counts_Euclid_S1_RAND_vy_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" """

z_alls = np.arange(0.1, 0.3, 0.05)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.05,2))
    print(cmd)
    os.system(cmd)

# S2 0.1-0.4
import os
import numpy as np
init_command = lambda z_lo, z_hi : """python compare_ct_vy.py "/home/comparat/sf_Shared/data/Euclid/data" "vy" "Counts_Euclid_S2_RAND_vy_"""+str(int(z_lo*100))+"""_z_"""+str(int(z_hi*100))+""".npy" """

z_alls = np.arange(0.1, 0.4, 0.05)
for z_lo in z_alls:
    cmd = init_command(np.round(z_lo,2), np.round(z_lo+0.05,2))
    print(cmd)
    os.system(cmd)




python plot_RS_Euclid_vy.py S0
python plot_RS_Euclid_vy.py S1
python plot_RS_Euclid_vy.py S2
