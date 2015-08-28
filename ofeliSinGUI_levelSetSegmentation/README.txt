Usage:

1- Level set algorithm configuration file is levelset.conf

2- Image to be segmented must be square, i. e. 500x500

3- Contour template must be os the same size of the image to be segmented

4- Compilation without GUI:

g++ main.cpp activecontour.cpp ac_withoutedges.cpp -o prueba -fopenmp -Dcimg_use_png -lpng -lz -Dcimg_display=0