#!/bin/bash
#
# Challenge1_Euclid0.iptd.rebin.fits  Challenge2_Euclid1.iptd.rebin.fits  Challenge3_Euclid1.iptd.rebin.fits
# Challenge1_Euclid1.iptd.rebin.fits  Challenge2_Euclid2.iptd.rebin.fits  Challenge3_Euclid2.iptd.rebin.fits
# Challenge1_Euclid2.iptd.rebin.fits  Challenge2_Euclid3.iptd.rebin.fits  Challenge3_Euclid3.iptd.rebin.fits
# Challenge1_Euclid3.iptd.rebin.fits  Challenge2_Euclid4.iptd.rebin.fits  Challenge3_Euclid4.iptd.rebin.fits
# Challenge1_Euclid4.iptd.rebin.fits  Challenge2_Euclid5.iptd.rebin.fits  Challenge3_Euclid5.iptd.rebin.fits
# Challenge2_Euclid0.iptd.rebin.fits  Challenge3_Euclid0.iptd.rebin.fits
#  Challenge2_Euclid1.iptd.rebin.fits  Challenge3_Euclid1.iptd.rebin.fits Challenge1_Euclid1.iptd.rebin.fits  Challenge2_Euclid2.iptd.rebin.fits  Challenge3_Euclid2.iptd.rebin.fits Challenge1_Euclid2.iptd.rebin.fits  Challenge2_Euclid3.iptd.rebin.fits  Challenge3_Euclid3.iptd.rebin.fits Challenge1_Euclid3.iptd.rebin.fits  Challenge2_Euclid4.iptd.rebin.fits  Challenge3_Euclid4.iptd.rebin.fits Challenge1_Euclid4.iptd.rebin.fits  Challenge2_Euclid5.iptd.rebin.fits  Challenge3_Euclid5.iptd.rebin.fits Challenge2_Euclid0.iptd.rebin.fits  Challenge3_Euclid0.iptd.rebin.fits


for file in Challenge1_Euclid3.iptd.rebin.fits
do
      echo "Launch Dawis on file ${file}"
      qsub qsub_dawis.sh -v ncl=${file:0:29},nf=${file}
      sleep 2
done
