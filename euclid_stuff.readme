jcomparat
<!u9K7b;w1Q;uuYdjyT$

Analyse euclid

ssh avec l'alias

ccin2p3

data in :
/sps/euclid/OU-LE3/VMPZ-ID/Q1


50847 files
1.4TB

some scripts will be shared

import glob

cd, '/sps/euclid/OU-LE3/VMPZ-ID/Q1/'

phzcat = numpy.array(glob.glob('/sps/euclid/OU-LE3/VMPZ-ID/Q1/*PHZCAT*'))
idf_ph = numpy.array([el.split('.')[-4] for el in phzcat ])
mercat = numpy.array(glob.glob('/sps/euclid/OU-LE3/VMPZ-ID/Q1/*FINAL-CAT*'))
idf_fc = numpy.array([el.split('.')[-4] for el in mercat ])

belongs = numpy.isin(idf_fc, idf_ph)
alias ccin2p3='ssh -Y jcomparat@cca.in2p3.fr'

scp jcomparat@cca.in2p3.fr:/pbs/home/j/jcomparat/merging_MERPHZ_Msk_Selection_JM_MV.py .

$DATA/act-mem-euclid-test.fits.gz
scp /home/comparat/sf_Shared/data/act-mem-euclid-test.fits.gz jcomparat@cca.in2p3.fr:
scp /home/comparat/sf_Shared/data/act-mem-euclid-test.fits.gz jcomparat@cca.in2p3.fr:

# setup of anaconda

curl -O https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh
bash ~/Anaconda3-2024.10-1-Linux-x86_64.sh

get yml file from :
https://gitlab.mpcdf.mpg.de/joco/eroconda

# ongoing
screen -S conda_env

anaconda3/bin/conda env create --file eroconda.yml


anaconda3/bin/conda activate eroconda
pip install ipython


python -m venv --system-site-packages "joco"
source "joco"/bin/activate
python -m pip install --user ipython
python -m pip install --user astropy
/pbs/home/j/jcomparat/joco/bin/python -m pip install --upgrade pip

python -m pip install --user ipython
python -m pip install --user astropy
python -m pip install --user scikit


python get_act_photom.py
python get_act_clu_gal.py


rsync -avz jcomparat@cca.in2p3.fr:~/*.fits.gz .
ls ACT_DR5_CLU_* > list.around
ls act-mem-euclid-test_AND_EUC* > list.mem
stilts tcat in=@list.mem ifmt=fits omode=out ofmt=fits out=ACT-mem-Euclid.fits
stilts tcat in=@list.around ifmt=fits omode=out ofmt=fits out=ACT-gal2Rlambda-Euclid.fits



If you would like to perform some tests on your side regarding the latest Q1 cluster detection catalogues, you would be very welcome and please feel free to use the directories below:
/sps/euclid/Users/sunayana/Q1_Nov24
In each sub-directory EDFN, EDFS, EDFF you should be able to find masks_corr, galaxy_cats   directories contained the filtered coverage masks and galaxy catalogues. In addition you may also find the output products from DET-CL (AMICO and PZWAV) and RICH-CL (RICH-CL products are output_members_[field_name] and output_richness_[field_name] ).
Any questions or uncertainties please let me know. You may of course copy the files you wish to your own scratch spaces to look at them in more detail. Happy investigating


cd /home/comparat/sf_Shared/data/Euclid/sunayana/Q1_Nov24/EDFF
rsync -avz jcomparat@cca.in2p3.fr:/sps/euclid/Users/sunayana/Q1_Nov24/EDFF/output*fits .
rsync -avz jcomparat@cca.in2p3.fr:/sps/euclid/Users/sunayana/Q1_Nov24/EDFF/*AMICO-DETECTIONS* .


cd /home/comparat/sf_Shared/data/Euclid/sunayana/Q1_Nov24/EDFS
rsync -avz jcomparat@cca.in2p3.fr:/sps/euclid/Users/sunayana/Q1_Nov24/EDFS/output*fits .
rsync -avz jcomparat@cca.in2p3.fr:/sps/euclid/Users/sunayana/Q1_Nov24/EDFS/*AMICO-DETECTIONS* .


cd /home/comparat/sf_Shared/data/Euclid/sunayana/Q1_Nov24/EDFN
rsync -avz jcomparat@cca.in2p3.fr:/sps/euclid/Users/sunayana/Q1_Nov24/EDFN/output*fits .
rsync -avz jcomparat@cca.in2p3.fr:/sps/euclid/Users/sunayana/Q1_Nov24/EDFN/*AMICO-DETECTIONS* .

<!u9K7b;w1Q;uuYdjyT$


Choice of single cluster:

ACT-CL J0425.5-3742
RA=66.38478798255967
Dec=-37.71422445364507

RM LS10
lambda=79.48263
z_lambda=0.32163247

RM DES Yr3
lambda=96
z=0.32154

ACT
SNR = 9.23
y_c = 1.22 pm0.13
M500c 4.69 pm +.87 -0.73
M200m = 8.11 pm +1.62 - 1.35

eROSITA
sm05_066129_020_ML00018_001_c030
66.3858177891451
-37.713307274872825
z_lambda=0.32473436
Lambda=83.29637

In R500 = 1186.1782958608703 kpc [1171, 1200]
LX 0.5-2 = 1.86 \pm 0.1 10^{44} erg/s
FX = 5.942453822466271E-13

In 300kpc, 0.5-2 keV
LX = 9.916824510828491E43 [9.34-10.49]
FX=3.168632651080985E-13 [2.99, 3.32]
CR=0.3663 [0.3458, 0.3846]
CT=446.725 [421, 469]
RA_Xfit = 66.38478
DEC_Xfit=-37.71483

theta_cluster

skyDistanceRadians(66.3858177891451*DEGREE_RADIANS, -37.713307274872825*DEGREE_RADIANS, RIGHT_ASCENSION*DEGREE_RADIANS, DECLINATION*DEGREE_RADIANS)/ARC_SECOND_RADIANS

kpc per arcsecond at the mean redshift : 6.19347576


topcat -stilts plot2plane \
   xpix=942 ypix=574 \
   xlabel='theta_cluster*6.19347576' ylabel= \
   xmin=0 xmax=3825 ymin=0 ymax=1157 \
   legend=true \
   in=/media/sf_Shared/data/Euclid/ACT-gal2Rlambda-Euclid.fits x='theta_cluster*6.19347576' binsize=150.0 \
   layer_01=Histogram \
      icmd_01='select <cl>' \
      leglabel_01='9: cl' \
   layer_02=Histogram \
      icmd_02='select "cl && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=0.2 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<0.3"' \
      color_02=green \
      leglabel_02='9: c02' \
   layer_03=Histogram \
      icmd_03='select "cl && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=0.3 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<0.4"' \
      color_03=grey \
      leglabel_03='9: c03' \
   layer_04=Histogram \
      icmd_04='select "cl && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=0.4 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<0.5 "' \
      color_04=magenta \
      leglabel_04='9: c04' \
   layer_05=Histogram \
      icmd_05='select "cl && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=0.5 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<0.6"' \
      color_05=cyan \
      leglabel_05='9: c05' \
   layer_06=Histogram \
      icmd_06='select "cl&& FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=0.6 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<0.7"' \
      color_06=orange \
      leglabel_06='9: c06' \
   layer_07=Histogram \
      icmd_07='select "cl&&FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=0.1 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<0.2"' \
      color_07=pink \
      leglabel_07='9: c01' \
   layer_08=Histogram \
      icmd_08='select "cl&&FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=0.7 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<0.8"' \
      color_08=yellow \
      leglabel_08='9: c07' \
   layer_09=Histogram \
      icmd_09='select cl&&FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=0.9&&FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<1' \
      leglabel_09='9: c09' \
   layer_10=Histogram \
      icmd_10='select "cl && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=0.8 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<0.9"' \
      color_10=blue \
      leglabel_10='9: c08' \
   layer_11=Histogram \
      icmd_11='select "cl&&FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=1.0 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<1.1"' \
      color_11=green \
      leglabel_11='9: c10' \
   layer_12=Histogram \
      icmd_12='select "cl&&FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER>=1.1 && FLUX_G_EXT_DECAM_APER/FLUX_R_EXT_DECAM_APER<1.2"' \
      color_12=grey \
      leglabel_12='9: c11'

topcat -stilts plot2plane \
   xpix=1460 ypix=621 \
   xlabel='R [kpc]' ylabel='N galaxies [1/kpc2]' fontsize=17 \
   xmin=0 xmax=3000 ymin=-4.0E-11 ymax=3.61E-4 \
   legend=true legpos=1.0,1.0 \
   in=histogram x='(LOW+HIGH)/2.' xerrhi='(HIGH-low)/2.' shading=auto size=5 \
   layer_01=XYError \
      y_01='c01_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' yerrhi_01='sqrt(c01_COUNT)/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      thick_01=1 color_01=orange \
      leglabel_01='0.1<g-r<0.2' \
   layer_02=Line \
      y_02='c01_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_02=orange \
      leglabel_02='0.1<g-r<0.2' \
   layer_03=XYError \
      y_03='c02_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' yerrhi_03='sqrt(c02_COUNT)/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      thick_03=1 \
      leglabel_03='0.2<g-r<0.3' \
   layer_04=Line \
      y_04='c02_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      leglabel_04='0.2<g-r<0.3' \
   layer_05=XYError \
      y_05='c03_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' yerrhi_05='sqrt(c03_COUNT)/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      thick_05=1 color_05=blue \
      leglabel_05='0.3<g-r<0.4' \
   layer_06=Line \
      y_06='c03_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_06=blue \
      leglabel_06='0.3<g-r<0.4' \
   layer_07=Mark \
      y_07='c04_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_07=green \
      leglabel_07='0.4<g-r<0.5' \
   layer_08=Line \
      y_08='c04_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_08=green \
      leglabel_08='0.4<g-r<0.5' \
   layer_09=XYError \
      y_09='c05_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' yerrhi_09='sqrt(c05_COUNT)/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      thick_09=1 color_09=grey \
      leglabel_09='0.5<g-r<0.6' \
   layer_10=Line \
      y_10='c05_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_10=grey \
      leglabel_10='0.5<g-r<0.6' \
   layer_11=Mark \
      y_11='c06_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_11=magenta \
      leglabel_11='0.6<g-r<0.7' \
   layer_12=Line \
      y_12='c06_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_12=magenta \
      leglabel_12='0.6<g-r<0.7' \
   layer_13=Mark \
      y_13='c07_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_13=cyan \
      leglabel_13='0.7<g-r<0.8' \
   layer_14=Line \
      y_14='c07_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_14=cyan \
      leglabel_14='0.7<g-r<0.8' \
   layer_15=Mark \
      y_15='c08_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_15=pink \
      leglabel_15='0.8<g-r<0.9' \
   layer_16=Line \
      y_16='c08_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_16=pink \
      leglabel_16='0.8<g-r<0.9' \
   layer_17=Mark \
      y_17='c09_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_17=yellow \
      leglabel_17='0.9<g-r<1.0' \
   layer_18=Line \
      y_18='c09_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_18=yellow \
      leglabel_18='0.9<g-r<1.0' \
   layer_19=Mark \
      y_19='c10_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_19=black \
      leglabel_19='1.0<g-r<1.1' \
   layer_20=Line \
      y_20='c10_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_20=black \
      leglabel_20='1.0<g-r<1.1' \
   layer_21=Mark \
      y_21='c11_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_21=999900 \
      leglabel_21='1.1<g-r<1.2' \
   layer_22=Line \
      y_22='c11_COUNT/(PI*HIGH*HIGH - PI*LOW*LOW)' \
      color_22=999900 \
      leglabel_22='1.1<g-r<1.2' \
   layer_23=Function \
      axis_23=Vertical fexpr_23=1180 color_23=black thick_23=2 \
      leglabel_23='R500 ~1180 kpc'
