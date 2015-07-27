#include "ac_withoutedges.hpp"
#include <stdio.h>
#include "CImg/CImg.h"
#include <sys/time.h>
#include <locale>
#include <iostream>

using namespace cimg_library;

void clean_boundaries(char* phi, char* phi_clean,int img_size,int img_width, int img_height);
bool isRedundantPointOfLin(char* phi, int x, int y, int img_width, int img_height);
bool isRedundantPointOfLout(char* phi, int x, int y, int img_width, int img_height);
void isolateIslands(const char*  new_phi, int width, int height);

static const int list_end = -9999999;

template <class charT, charT sep>
class punct_facet: public std::numpunct<charT> {
protected:
    charT do_decimal_point() const { return sep; }
};

int main(int argc, char* argv[])
{
	CImg<unsigned char> image("imagenes/pruebaISLAS2.png" ); //Imagen a segmentar
	int width1 = image.width(), height1=image.height(); //variables de la imagen a segmentar
	int size1 = width1*height1;//tamaño de la imagen a segmentar
	unsigned char* img = new unsigned char[size1]; //Imagen en forma de vector
	
	CImg<unsigned char> phi_img("imagenes/initial_phi.png"); //Imagen de la primera phi
	int width2 = phi_img.width(), height2=phi_img.height(); //variables de la imagen de la primera phi
	int size2 = width2*height2;//tamaño de la imagen phi
	char* phi = new char[size2]; //Imagen de la phi en forma de vector
	char* phi_clean = new char[size2]; //Imagen de la phi en forma de vector

	//Listas
    const ofeli::list<int>* Lout1;
    const ofeli::list<int>* Lin1;
	
	//Otras variables
	bool hasSmoothingCycle1 = true;
  	int Na1 = 30;
    int Ns1 = 3;
    int lambda_out1 = 1;
    int lambda_in1 = 1;
    int kernel_curve1 = 5;
    double std_curve1 = 2.0;
	struct timeval time1;

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
			phi_clean[i]=1;
		}
		else
		{
			phi[i] = -1;
			phi_clean[i]=-1;
		}
	}
	clean_boundaries(phi,phi_clean,size2,width2,height2);		

	/******************************** ALGORITMO ********************************/
	ofeli::ActiveContour* ac;
	ac = new ofeli::ACwithoutEdges(img, width1, height1, phi_clean, hasSmoothingCycle1, kernel_curve1, std_curve1, Na1, Ns1, lambda_out1, lambda_in1);
	
	gettimeofday(&time1, NULL);	
	double t1 = time1.tv_sec * 1000000 + time1.tv_usec;

	ac->evolve(); //Evolución del contorno
	
	gettimeofday(&time1, NULL);	
	double t2 = time1.tv_sec * 1000000 + time1.tv_usec;
	
	double totalTime = (t2-t1)/1000000;	
	
	std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
	std::cout <<  totalTime  << std::endl;
	/***************************************************************************/
	
	//Comprobación de resultados
	Lout1 = &ac->get_Lout();
	Lin1 = &ac->get_Lin();
	const char*  new_phi = ac->get_phi();

	//DISPLAY DE LOS BORDES ENCONTRADOS
	CImg<float> imagen_a_color("imagenes/pruebaISLAS.png"); 
	for(int i=0; i < size1; i++)
	{
		int x,y;
		y = i/width1;
		x = i - (y*width1);
		if(new_phi[i] == -1 )
		{
			//Rojo
			imagen_a_color(x, y, 0,0) = 255;
			imagen_a_color(x, y, 0,1) = 0;
			imagen_a_color(x, y, 0,2) = 0;		
		}
		if(new_phi[i] == 1 )
		{
			//Azul
			imagen_a_color(x, y, 0,0) = 0;
			imagen_a_color(x, y, 0,1) = 26;
			imagen_a_color(x, y, 0,2) = 255;
		}
	}
	
	//imagen_a_color.display("RESULTADO");
	imagen_a_color.save("result.png");
	
	isolateIslands(new_phi, width1, height1);
	
/*	
	//IMPRIMIR LA PHI
	printf("TAM: %d\n", Lout1->size());
	std::cout << "PHI: "<< std::endl;
	int cont=0;
	for(int i=0; i < size1; i++)
	{
		if((i % width1 )== 0)
		{
			std::cout << "<-- fila " << cont << std::endl;
			cont++;
		} 
		printf(", %d", new_phi[i]);
	}
	std::cout << std::endl;
*/
	return 0;	
}












void clean_boundaries(char* phi, char* phi_clean,int img_size,int img_width, int img_height)
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




void isolateIslands(const char*  phi, int width, int height)
{
    int cont=0;
    int size=width*height, i=0;
	CImg<float> islasRecorridas("imagenes/pruebaISLAS2.png"); 
    bool* visitedIslands = new bool[size];
    for(int a=0; a < size; a++) visitedIslands[a] = false;
	
	//Find all the islands in the image
	for(int i=width+1; i < size; i++)
	{
	    ofeli::list<int>* islandPoints = new ofeli::list<int>();

		//Find a island 
		for(; i < size; i++)
		{
			if(phi[i] == -1 && !visitedIslands[i]) break;
		}
		int x,y,min_x=99999999,min_y=99999999,max_x=-1,max_y=-1;
		
		//Mark the first point of the island
		visitedIslands[i] = true;
		islandPoints->push_front(i);
		int islandPoint = i;
		int numberOfPoints=0;
		cont++;
		
		//Take the border of the island
		while(1)
		{	
			numberOfPoints++;
			
			y = islandPoint/width;
			x = islandPoint- (y*width);

			//Take borders to make the image later
			if(x < min_x) min_x = x;
			if(x > max_x) max_x = x;
			if(y < min_y) min_y = y;
			if(y > max_y) max_y = y;

			int offset =-1;
			//To find the next border point in current pixel neighbor prioritizing adyacent ones
			if(phi[(x-1) +  width * y] == -1 && !visitedIslands[(x-1) +  width * y])	offset = (x-1) +  width * y; 			//Adyacent neighbor
			else if(phi[(x+1) +  width * y] == -1 && !visitedIslands[(x+1) +  width * y])	offset = (x+1) +  width * y;		//Adyacent neighbor			
			else if(phi[x +  width * (y-1)] == -1 && !visitedIslands[x +  width * (y-1)])	offset = x +  width * (y-1);		//Adyacent neighbor					
			else if(phi[x +  width * (y+1)] == -1 && !visitedIslands[x +  width * (y+1)])	offset = x +  width * (y+1);		//Adyacent neighbor		
			else if(phi[(x-1) +  width * (y+1)] == -1 && !visitedIslands[(x-1) +  width * (y+1)])	offset = (x-1) +  width * (y+1);
			else if(phi[(x-1) +  width * (y-1)] == -1 && !visitedIslands[(x-1) +  width * (y-1)])	offset = (x-1) +  width * (y-1);
			else if(phi[(x+1) +  width * (y-1)] == -1 && !visitedIslands[(x+1) +  width * (y-1)])	offset = (x+1) +  width * (y-1);
			else if(phi[(x+1) +  width * (y+1)] == -1 && !visitedIslands[(x+1) +  width * (y+1)])	offset = (x+1) +  width * (y+1);
			
			//When we've make a round or the island
			if(offset == -1) break;
			
			//Save the point
			islandPoints->push_front(offset);
			visitedIslands[offset] = true;

			//Continue with the next point
			islandPoint = offset;			
		}//end while
		
		//To fix some errors
		if(numberOfPoints > 4)
		{
			//Create the island
			int tam1= max_x-min_x +1;
			int tam2 = max_y-min_y+1;
			CImg<float> island(tam1, tam2, 1,3,255); 
			
			//Sets the border of the island 
			for(ofeli::list<int>::iterator position = islandPoints->begin(); !position.end(); ++position)
			{
				int X,Y;
				Y = *position/width;
				X = *position-(Y*width);

				//Para ver el recorrido que se ha hecho
				islasRecorridas(X, Y, 0,0) = 255;
				islasRecorridas(X, Y, 0,1) = 0;
				islasRecorridas(X, Y, 0,2) = 0;
				
				//Reajust the positions
				X = X - min_x;
				Y = Y - min_y;

				island(X, Y, 0,0) = 0;
				island(X, Y, 0,1) = 0;
				island(X, Y, 0,2) = 0;

			}	
			/*
			//Fill up the rest of the island
			for(int island_y=0; island_y < tam2; island_y++)
			{		
				int border = -1, b_x= -1, b_y= -1, touchingBorder = -1;						
				for(int island_x=0; island_x < tam1; island_x++)
				{
					//std::cout << "island_x: " << island_x << std::endl;
					//std::cout << "island_y: " << island_y << std::endl;

					float f = *island.data(island_x, island_y, 0, 0);

					if(f == 0 && border == -1)
					{
						border = 1;
						b_x= island_x;
						b_y= island_y;
						touchingBorder=1;
					}
					//Founded other border
					else if(f == 0 && border == 1 && touchingBorder == -1)
					{	
						//Paint a line between the two borders
						for(int l=b_x; l < island_x; l++)
						{
							island(l, b_y, 0,0) = 0;
							island(l, b_y, 0,1) = 0;
							island(l, b_y, 0,2) = 0;
						}
						border = -1;
					}
					else if(f != 0)
					{
						touchingBorder=-1;
					}				
				}	
			}*/

			//Save the island as a image
			char nombre[25];
			sprintf (nombre, "islas/island%d.png", cont);
			island.save(nombre);
			islasRecorridas.save("recorrido.png");			
		}
		
	}//END FOR
}
