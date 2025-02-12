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


#
# MAKE A NICE FIGURE
#

#
# Do the same with Euclid galaxies (deep !)
#

python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_010_z_011.npy' 'CT_s0_gr_010_z_011.png'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_011_z_012.npy' 'CT_s0_gr_011_z_012.png'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_012_z_013.npy' 'CT_s0_gr_012_z_013.png'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_013_z_014.npy' 'CT_s0_gr_013_z_014.png'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_014_z_015.npy' 'CT_s0_gr_014_z_015.png'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_015_z_016.npy' 'CT_s0_gr_015_z_016.png'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_016_z_017.npy' 'CT_s0_gr_016_z_017.png'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_017_z_018.npy' 'CT_s0_gr_017_z_018.png'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_018_z_019.npy' 'CT_s0_gr_018_z_019.png'
python compare_ct.py '/home/comparat/sf_Shared/data/legacysurvey/dr10/sweep' 'Counts_DATA_S0_RAND_gr_019_z_020.npy' 'CT_s0_gr_019_z_020.png'
