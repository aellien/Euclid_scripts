#!/bin/bash
# Crop image (hand made with box region)

path_simulations=/home/aellien/Euclid_ICL/simulations/out2/
path_scripts=/home/aellien/Euclid_ICL/Euclid_scripts

for cluster_num in 365000132000018 360009933000018 373000139000019
do
    cd ${path_simulations}${cluster_num}/vignets
    for full_mosaic in ${path_simulations}${cluster_num}/EUC_NIR_W-STK-IMAGE_H_z_*
    do
        echo $full_mosaic
        vignet_num=0
        fn="$(basename -- $full_mosaic)"
        while read -r ra dec
        do 
            
            vignet_num=$(($vignet_num+1))
            out=${fn:0:-8}_vignet_${vignet_num}.fits
            echo $out
            astcrop --mode=wcs --center=$ra,$dec \
                    --width=0.12,0.12 \
                    -h1 \
                    --output=$out \
                    --primaryimghdu \
                    $full_mosaic

        done < /home/aellien/Euclid_ICL/simulations/out2/cl_coo.txt
    done
done
cd $path_scripts
