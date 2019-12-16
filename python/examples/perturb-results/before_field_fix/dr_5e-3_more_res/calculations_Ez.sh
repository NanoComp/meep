#!/bin/sh

DR=0.005
for RESOLUTION in 10 12 15 19 24 29 36 45 55 69 85 105 130 161 200
do
  python ../../ring-cyl-perturbation-theory.py -res $RESOLUTION -dr $DR | grep -e component: -e res: -e dr: -e w: -e dwdR: > ring-cyl-perturbation-theory.Ez.perp_False_res_${RESOLUTION}_dr_${DR}.dat
  # python ../../ring-cyl-perturbation-theory.py -perp True -res $RESOLUTION -dr $DR | grep -e component: -e res: -e dr: -e w: -e dwdR: > ring-cyl-perturbation-theory.Hz.perp_True_res_${RESOLUTION}_dr_${DR}.dat
done
