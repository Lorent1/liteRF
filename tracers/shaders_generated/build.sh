#!/bin/sh
glslangValidator -V kernel2D_RayMarch.comp -o kernel2D_RayMarch.comp.spv -DGLSL -I.. -I/home/petr/dev/liteRF/external/json/include -I/home/petr/kernel_slicer/TINYSTL -I/home/petr/dev/liteRF/external/LiteMath 
