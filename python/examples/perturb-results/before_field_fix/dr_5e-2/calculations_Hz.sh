#!/bin/sh

DR=0.05
for RESOLUTION in 100 119 141 168 200
do
  # python ../../ring-cyl-perturbation-theory.py -res $RESOLUTION -dr $DR | grep -e component: -e res: -e dr: -e w: -e dwdR: > ring-cyl-perturbation-theory.Ez.perp_False_res_${RESOLUTION}_dr_${DR}.dat
  python ../../ring-cyl-perturbation-theory.py -perp True -res $RESOLUTION -dr $DR | grep -e component: -e res: -e dr: -e w: -e dwdR: > ring-cyl-perturbation-theory.Hz.perp_True_res_${RESOLUTION}_dr_${DR}.dat
done
