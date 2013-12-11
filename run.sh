#!/bin/bash

cd app_proiect

export OMP_NUM_THREADS=$1
export OMP_SCHEDULE=dynamic

time ./segment 0.6 300 300 landscape3.ppm out.ppm
