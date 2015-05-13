#include "ac_withoutedges.hpp"
#include <stdio.h>
#include "CImg/CImg.h"

using namespace cimg_library;

void clean_boundaries(char* phi, char* phi_clean,int img_size,int img_width, int img_height, int* Lout, int* Lin);
bool isRedundantPointOfLin(char* phi, int x, int y, int img_width, int img_height);
bool isRedundantPointOfLout(char* phi, int x, int y, int img_width, int img_height);

static const int list_end = -9999999;

int main(int argc, char* argv[])
{
	
	cimg_library::CImg<unsigned char> image("imagenes/prueba.png"); //Imagen a segmentar
	int width1 = image.width(), height1=image.height(); //variables de la imagen a segmentar
	int size1 = width1*height1;//tamaño de la imagen a segmentar
	unsigned char* img = new unsigned char[size1]; //Imagen en forma de vector
	
	cimg_library::CImg<unsigned char> phi_img("imagenes/initial_phi.bmp"); //Imagen de la primera phi
	int width2 = phi_img.width(), height2=phi_img.height(); //variables de la imagen de la primera phi
	int size2 = width2*height2;//tamaño de la imagen phi
	char* phi = new char[size2]; //Imagen de la phi en forma de vector
	char* phi_clean = new char[size2]; //Imagen de la phi en forma de vector

	//Listas
	int* Lout = new int[size2+1];
	int* Lin = new int[size2+1];

	//Otras variables
	bool hasSmoothingCycle1 = true;
  	int Na1 = 30;
    int Ns1 = 3;
    int lambda_out1 = 1;
    int lambda_in1 = 1;
    int kernel_curve1 = 5;
    double std_curve1 = 2.0;
	
	//Coger los valores de la imagen a segmentar y se rellena el array
	for(int i=0; i < size1; i++)
	{
		int x,y;
		y = i/width1;
		x = i - (y*width1);
		
		img[i] = *image.data(x,y);
	}
	
	//Coger los valores de la imagen con el phi inicial y se rellena el array
	for(int i=0; i < size2; i++)
	{
		int x,y;
		y = i/width2;
		x = i - (y*width2);
		
		if(*phi_img.data(x,y) < 128)
		{
			phi[i] = 1;
		}
		else
		{
			phi[i] = -1;
		}
	}
	

	//Relleno de las listas 
	Lout[0] = list_end;
	Lin[0] = list_end;
	clean_boundaries(phi,phi_clean,size2,width2,height2, Lout,Lin);
	
	int n_out = 0;
	int n_in = 0;
	for( int offset = 0; offset < size2; offset++ )
	{
		if( phi_clean[offset] == 1 )
		{
			Lout[n_out++] = offset;
		}
		if( phi_clean[offset] == -1 )
		{
			Lin[n_in++] = offset;
		}
	}
	Lout[n_out] = list_end;
	Lin[n_in] = list_end;
	
	for(int k=0; k < img_size+1; k++)
	{
		if(Lout[k] == -9999999) break;
		printf(", %d", Lout[k]);

	}
	
	ofeli::ActiveContour* ac;

	//ac = new ofeli::ACwithoutEdges(img, width, height, phi_clean, hasSmoothingCycle1, kernel_curve1, std_curve1, Na1, Ns1, lambda_out1, lambda_in1);

	return 0;
	
}

void clean_boundaries(char* phi, char* phi_clean,int img_size,int img_width, int img_height, int* Lout, int* Lin)
{
    if( phi != NULL && phi_clean != NULL )
    {
        int offset, x, y; // position

        for( offset = 0; offset < img_size; offset++ )
        {
            if( phi[offset] == 1 )
            {
                // offset = x+y*img_width <==>
                y = offset/img_width;
                x = offset-y*img_width;

                // nettoyage si la condition de voisinage aux frontières n'est pas respectée
                if( isRedundantPointOfLout(phi,x,y,img_width,img_height) )
                {
                    phi_clean[offset] = 3;
                }
                else
                {
                    phi_clean[offset] = 1;
                }
            }

            if( phi[offset] == -1 )
            {
                // offset = x+y*img_width <==>
                y = offset/img_width;
                x = offset-y*img_width;

                if( isRedundantPointOfLin(phi,x,y,img_width,img_height ) )
                {
                    phi_clean[offset] = -3;
                }
                else
                {
                    phi_clean[offset] = -1;
                }
            }               

        int n_out = 0;
        int n_in = 0;
        for( int offset = 0; offset < img_size; offset++ )
        {
            if( phi_clean[offset] == 1 )
            {
                Lout[n_out++] = offset;
            }
            if( phi_clean[offset] == -1 )
            {
                Lin[n_in++] = offset;
            }
        }
        Lout[n_out] = list_end;
        Lin[n_in] = list_end;

    }
    return;
}
	
}
	
	
	
bool isRedundantPointOfLin(char* phi, int x, int y, int img_width, int img_height)
{
    // if ∃ a neighbor ∈ Lout | ∈ Rout

    if( y-1 >= 0 )
    {
        if( phi[ x+(y-1)*img_width ] >= 0 )
        {
            return false; // is not redundant point of Lin
        }
    }
    if( x-1 >= 0 )
    {
        if( phi[ (x-1)+y*img_width ] >= 0 )
        {
            return false; // is not redundant point of Lin
        }
    }
    if( x+1 < img_width )
    {
        if( phi[ (x+1)+y*img_width ] >= 0 )
        {
            return false; // is not redundant point of Lin
        }
    }
    if( y+1 < img_height )
    {
        if( phi[ x+(y+1)*img_width ] >= 0 )
        {
            return false; // is not redundant point of Lin
        }
    }

    // ==> ∀ neighbors ∈ Lin | ∈ Rin
    return true; // is redundant point of Lin
}

	
	
	
bool isRedundantPointOfLout(char* phi, int x, int y, int img_width, int img_height)
{
    // if ∃ a neighbor ∈ Lin | ∈ Rin

    if( y-1 >= 0 )
    {
        if( phi[ x+(y-1)*img_width ] <= 0 )
        {
            return false; // is not redundant point of Lout
        }
    }
    if( x-1 >= 0 )
    {
        if( phi[ (x-1)+y*img_width ] <= 0 )
        {
            return false; // is not redundant point of Lout
        }
    }
    if( x+1 < img_width )
    {
        if( phi[ (x+1)+y*img_width ] <= 0 )
        {
            return false; // is not redundant point of Lout
        }
    }
    if( y+1 < img_height )
    {
        if( phi[ x+(y+1)*img_width ] <= 0 )
        {
            return false; // is not redundant point of Lout
        }
    }

    // ==> ∀ neighbors ∈ Lout | ∈ Rout
    return true; // is redundant point of Lout
}