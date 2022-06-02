#! /bin/bash

# for filepath in images/*.png; 
# do 
#     filename=$(echo $filepath | sed "s/images\/\(.*\).png/\1/")
#     build/img_dft $filepath build/output/{$filename\_dft,$filename\_phase}.png --log;
# done

build/img_dft images/double_slit_small_norm.png build/output/double_slit_small_norm_mag.png build/output/double_slit_small_norm_phase.png
build/img_dft images/double_slit_small_close.png build/output/double_slit_small_close_mag.png build/output/double_slit_small_close_phase.png
build/img_dft images/double_slit_small_widest.png build/output/double_slit_small_widest_mag.png build/output/double_slit_small_widest_phase.png
build/img_dft images/double_slit_norm_widest.png build/output/double_slit_norm_widest_mag.png build/output/double_slit_norm_widest_phase.png
build/img_dft images/hole.png build/output/hole_mag.png build/output/hole_phase.png
build/img_dft images/hole_tiny.png build/output/hole_tiny_mag.png build/output/hole_tiny_phase.png

build/img_dft images/webb.png build/output/webb_mag_n.png build/output/webb_phase_n.png --noshift
build/img_dft images/webb.png build/output/webb_mag_log.png build/output/webb_phase_log.png --log --noshift
build/img_dft images/webb.png build/output/webb_mag.png build/output/webb_phase.png --log

build/img_dft images/dune.png build/output/dune_mag.png build/output/dune_phase.png --log
