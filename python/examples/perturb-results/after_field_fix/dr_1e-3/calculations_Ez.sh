#!/bin/sh

DR=0.001
for RESOLUTION in 100 108 117 126 136 147 159 171 185 200
do
  python ../../../ring-cyl-perturbation-theory.py -res $RESOLUTION -dr $DR | grep -e component: -e res: -e dr: -e w: -e dwdR: > ring-cyl-perturbation-theory.Ez.perp_False_res_${RESOLUTION}_dr_${DR}.dat
  # python ../../../ring-cyl-perturbation-theory.py -perp True -res $RESOLUTION -dr $DR | grep -e component: -e res: -e dr: -e w: -e dwdR: > ring-cyl-perturbation-theory.Hz.perp_True_res_${RESOLUTION}_dr_${DR}.dat
done
