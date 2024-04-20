#!/bin/bash

# Clean up first
rm src/*cxx.so
rm src/*cxx.d
rm src/*dict.pcm

# Now run code
input=$1

echo "Running analysis"
echo "Input edm4eic data is at [ ${input} ]"

# Extract the input filename without the extension
input_filename=$(ls "$input")
output_filename="${input_filename%.*}_output.root"

echo "Output at [ ${output_filename} ]"

root -b -q src/uchannelrho.cxx+\(\"${input}\",\"${output_filename}\"\)

# plotting benchmark figures.
# figures are stored under figures.
root -b -q macros/plot_rho_physics_benchmark.C\(\"${output_filename}\"\)

