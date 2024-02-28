#!/bin/bash
# Crop image (hand made with box region)

path_simulations=/home/aellien/Euclid_ICL/simulations/out1/
path_scripts=/home/aellien/Euclid_ICL/Euclid_scripts

#205000105000020 
#54000066009352
for cluster_num in 205000105000020 
do
    cd ${path_simulations}${cluster_num}/vignets
    for full_mosaic in ${path_simulations}${cluster_num}/EUC_NIR_W-STK-IMAGE_H_z_*
    do
        echo $full_mosaic
        z=${full_mosaic:82:3}
        vignet_num=0
        while read -r ra dec
        do 
           
            vignet_num=$(($vignet_num+1))
            echo z_${z}_vignet_${vignet_num}.fits
            astcrop --mode=wcs --center=$ra,$dec \
                    --width=0.12,0.12 \
                    -h1 \
                    --output=z_${z}_vignet_${vignet_num}.fits \
                    --primaryimghdu \
                    $full_mosaic

        done < /home/aellien/Euclid_ICL/simulations/out1/cl_coo.txt
    done
done
cd $path_scripts
