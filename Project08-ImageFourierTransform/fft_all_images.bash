#! /bin/bash

for filepath in images/*.png; 
do 
    filename=$(echo $filepath | sed "s/images\/\(.*\).png/\1/")
    build/img_dft $filepath build/output/{$filename\_dft,$filename\_phase}.png --log;
done
