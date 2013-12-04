#!/bin/bash

cd app_proiect

export OMP_NUM_THREADS=$1

time ./segment 0.6 300 300 landscape2.ppm out.ppm
