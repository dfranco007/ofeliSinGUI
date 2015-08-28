#include "ac_withoutedges.hpp"
#include <stdio.h>
#include <stdlib.h>
#include "CImg/CImg.h"
#include <sys/time.h>
#include <locale>
#include <iostream>
#include <fstream>
#include <string>

using namespace cimg_library;

void clean_boundaries(char* phi, char* phi_clean,int img_size,int img_width, int img_height);
bool isRedundantPointOfLin(char* phi, int x, int y, int img_width, int img_height);
bool isRedundantPointOfLout(char* phi, int x, int y, int img_width, int img_height);
void isolateIslands(char*  phi, int width, int height, int islandInnerValue, int islandInnerFarValue, int * innerIslands, int * allIslands, int* innerSize, int* maxWidth, int* maxHeight);
void fillUpIsland(char*  phi,int width, int min_x, int min_y, int islandInnerFarValue, CImg<float>* island, int innerPoint_x, int innerPoint_y);
bool goOverAnIsland(char*  phi, int islandPoint, int width, int height, int* min_x, int* min_y, int* max_x, int * max_y, int *inner_min_x, int *inner_min_y, int *inner_max_x, int *inner_max_y, ofeli::list<int>* islandPoints, bool* visitedIslands , int islandInnerValue);
void  adjustImagesAndCreatePSDFiles(int maxWidth, int maxHeight,int originalWidth, int originalHeight, int allIslands);

template <class charT, charT sep>
class punct_facet: public std::numpunct<charT> {
protected:
    charT do_decimal_point() const { return sep; }
};

int main(int argc, char* argv[])
{	
	//Create directories
	system("mkdir islands islands/rescaledAllIslands islands/allIslands_inOriginalImage islands/innerIslands_PSD_files islands/allIslands");

	//File variable configuration
	std::string line;
	std::ifstream file;
	file.open("levelset.conf");
	if(!file)
	{
		fprintf(stderr,"levelset.conf file not found.");
		return 1;
	}	
	//Skip first 5 lines cause are comments
	for(int i =1; i <6; i++)   getline (file,line);

	//Background and islands configuration
	bool inverse;
	int islandInnerValue,islandInnerFarValue;
	getline (file,line);
	line = line.substr(0,line.size()-1); //To remove the last CR character
	if(line.compare("true") ==0) inverse = true;
	else inverse = false;
	if(inverse)
	{
		islandInnerValue = 1;
		islandInnerFarValue = 3;
	}
	else
	{
		islandInnerValue = -1;
		islandInnerFarValue = -3;
	}

	getline (file,line);
	std::string imageName = line.substr(0,line.size()-1); //To remove the last CR character
	CImg<unsigned char> image(imageName.c_str()); //Image to be segmented and its variables
	int width = image.width(), height=image.height(); 
	int size = width*height;
	unsigned char* img = new unsigned char[size]; 
	
	getline (file,line);
	line = line.substr(0,line.size()-1); //To remove the last CR character
	CImg<unsigned char> phi_img(line.c_str()); //Contour template and its variables
	char* phi = new char[size];
	char* phi_clean = new char[size]; 

	//Test images size
	if( (image.width() != phi_img.width()) || (image.height() != phi_img.height()) )
	{
		fprintf(stderr,"The image to be segmented and its contour template have different sizes.");
		return 2;
	}
	if(image.width() != image.height())
	{
		fprintf(stderr,"The image to be segmented is not square");
		return 3;	
	}
	//Variables to make the segmentation
    const ofeli::list<int>* Lout1;
    const ofeli::list<int>* Lin1;
	getline (file,line);
	line = line.substr(0,line.size()-1); //To remove the last CR character
	bool hasSmoothingCycle1;
	if (line.compare("true") ==0)	hasSmoothingCycle1 = true;
	else	hasSmoothingCycle1 = false;
	getline (file,line);
	int Na1 = std::atoi(line.c_str()); 
	getline (file,line);
    int Ns1 = std::atoi(line.c_str()); 
	getline (file,line);	
    int lambda_out1 = std::atoi(line.c_str()); 
	getline (file,line);	
    int lambda_in1 = std::atoi(line.c_str()); 
	getline (file,line);	
    int kernel_curve1 = std::atoi(line.c_str()); 	
	getline (file,line);
    double std_curve1 = std::atof(line.c_str());
	struct timeval time1;

	//Fill up the image array
	for(int i=0; i < size; i++)
	{
		int x,y;
		y = i/width;
		x = i - (y*width);
		
		img[i] = *image.data(x,y);
	}
	
	//Fill up the phi of the level set with the contour template
	for(int i=0; i < size; i++)
	{
		int x,y;
		y = i/width;
		x = i - (y*width);
		
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
	clean_boundaries(phi,phi_clean,size,width,height);
	

/******************************** EVOLUTION ALGORITHM ********************************/
	
	//Create the contour 
	ofeli::ActiveContour* ac;
	ac = new ofeli::ACwithoutEdges(img, width, height, phi_clean, hasSmoothingCycle1, kernel_curve1, std_curve1, Na1, Ns1, lambda_out1, lambda_in1);
	
	gettimeofday(&time1, NULL);	
	double t1 = time1.tv_sec * 1000000 + time1.tv_usec;

	ac->evolve(); //Evolve the contour
	
	//Save the results
	char*  new_phi= ac->get_phi();
		
	gettimeofday(&time1, NULL);	
	double t2 = time1.tv_sec * 1000000 + time1.tv_usec;
	double totalTime = (t2-t1)/1000000;	
	
	std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
	std::cout <<  "Execution time: " << totalTime  << " s" << std::endl;
		
	//Calculate convering
	double covering = ac->calculateCovering(islandInnerValue);
	std::cout << "Covering: " << covering << "%" << std::endl;

	//Modify the border of the image to take in a easier way the islands of the border	
	for(int x=0; x < width; x++)//Top margin
		if(new_phi[x] == islandInnerFarValue) new_phi[x] = islandInnerValue;
	for(int y=0; y < height; y++)//Right margin
		if(new_phi[(width * y) -1] == islandInnerFarValue) new_phi[(width * y) -1] = islandInnerValue;
	for(int y=0; y < height; y++)//Left margin	
		if(new_phi[width * y] == islandInnerFarValue) new_phi[width * y] = islandInnerValue;
	for(int x=0; x < width; x++)//Bottom margin	
		if(new_phi[( (height-1) * width) + x +1] == islandInnerFarValue) new_phi[( (height-1) * width) + x +1] = islandInnerValue;	
	
/** SAVE THE RESULT IMAGE **/
	CImg<float> imagen_a_color(imageName.c_str()); //Must be the same image as the image to be segmented
	for(int i=0; i < size; i++)
	{
		int x,y;
		y = i/width;
		x = i - (y*width);
		if(new_phi[i] == -1 )
		{
			//Red represents Lin
			imagen_a_color(x, y, 0,0) = 255;
			imagen_a_color(x, y, 0,1) = 0;
			imagen_a_color(x, y, 0,2) = 0;		
		}
		if(new_phi[i] == 1 )
		{
			//Blue represents Lout
			imagen_a_color(x, y, 0,0) = 0;
			imagen_a_color(x, y, 0,1) = 26;
			imagen_a_color(x, y, 0,2) = 255;
		}
	}
	
	//imagen_a_color.display("RESULTADO");
	imagen_a_color.save("result.png");	
	
	//Isolate all islands and count them
	int *innerIslands = new int(0), *allIslands=new int(0), *innerSize =new int(0), *maxWidth =new int(0), *maxHeight =new int(0);
	isolateIslands(new_phi, width, height, islandInnerValue, islandInnerFarValue, innerIslands, allIslands, innerSize, maxWidth, maxHeight);
    std::cout << "Number of all islands: " << *allIslands << std::endl;
	std::cout << "Number of inner islands: " << *innerIslands << std::endl;
	
	float allIslandsDensity = *allIslands/(float)size;
	std::cout << "All islands density: " << allIslandsDensity << std::endl;
	float innerIslandsDensity = *innerIslands / (float)*innerSize;
	std::cout << "Inner islands density: " << innerIslandsDensity << std::endl;
	std::cout << "Average islands density(all and inner): " << (innerIslandsDensity + allIslandsDensity)/2 << std::endl;
	
    //Change the all images to the same size and squares them to calculate PSD correctly. Also creates PSD files.
	std::cout << "Preparing images to make PSD calculation... " << std::endl;
	adjustImagesAndCreatePSDFiles(*maxWidth, *maxHeight, width, height, *allIslands);
	system("cp islands/innerIslands_PSD_files/* islands/allIslands");

	//PSD
	//system("java -jar psd.jar islands/innerIslands_PSD_files/i" + *allIslands + "false true");
	//system("java -jar psd.jar islands/allIslands 1 false true");
	//system("java -jar psd.jar islands/innerIslands 1 false true");

    delete innerIslands; delete allIslands; delete innerSize;
    
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

    // ==> ∀ neighbours ∈ Lout | ∈ Rout
    return true; // is redundant point of Lout
}




void isolateIslands(char*  phi, int width, int height, int islandInnerValue, int islandInnerFarValue, int * innerIslands, int * allIslands, int* innerSize, int* maxWidth, int* maxHeight)
{
    int size=width*height, i=0, *inner_min_x =  new int(99999999), *inner_min_y = new int(99999999) , *inner_max_x = new int(-1) , *inner_max_y = new int(-1) ;
    bool* visitedIslands = new bool[size];
	bool innerIsland;
	CImg<float> innerIslandOnlyImage(width, height, 1,3,255);
	CImg<float> allIslandsInOriginalImage(width, height, 1,3,255);

	std::ofstream innerIslandsList("innerIslandsList.txt");

    for(int a=0; a < size; a++) visitedIslands[a] = false;

	//Find all the islands in the image
	for(int i=0; i < size; i++)
	{
	    ofeli::list<int>* islandPoints = new ofeli::list<int>();
		innerIsland = false;
		
		//Find a island 
		for(; i < size; i++)
		{
			if(phi[i] == islandInnerValue && !visitedIslands[i]) break;
		}
		int x,y,*min_x = new int(99999999),*min_y= new int(99999999),*max_x= new int(-1),*max_y = new int(-1);

		//Go over the island to be isolated
		if ( !goOverAnIsland(phi, i, width, height, min_x, min_y, max_x, max_y, inner_min_x, inner_min_y, inner_max_x, inner_max_y, islandPoints, visitedIslands, islandInnerValue ))
        {
			innerIsland = true;
        }				
		        
		//Create the island
		*innerSize = (*inner_max_x - *inner_min_x + 1) * (*inner_max_y - *inner_min_y + 1); //To calculate the inner density
		int tam1= *max_x-*min_x +1, tam2 = *max_y-*min_y+1, pointNumber=0;
        if(*maxWidth < tam1)  *maxWidth= tam1;
        if(*maxHeight < tam2) *maxHeight = tam2;
		CImg<float> island(tam1, tam2, 1,3,255);
		CImg<float> singleIslandInOriginalImage(width, height, 1,3,255);

		//Inner point control variables
		bool innerPointFounded = false;
		int innerPoint_x=-1, innerPoint_y=-1;
				
		//Sets the border of the island 
		for(ofeli::list<int>::iterator position = islandPoints->begin(); !position.end(); ++position)
		{					
			int X,Y;
			Y = *position/width;
			X = *position-(Y*width);
			
			//Take a inner point of the island 
			if(!innerPointFounded)
			{					
				if(phi[(X-1) +  width * Y] == islandInnerFarValue && (X -1 - *min_x) >= 0)
				{
					innerPoint_x = X -1 - *min_x;
					innerPoint_y= Y- *min_y;
					innerPointFounded=true;
				}
				else if(phi[(X+1) +  width * Y] == islandInnerFarValue &&  (X+1 -*min_x) < island.width() )
				{
					innerPoint_x = X+1 -*min_x;
					innerPoint_y= Y -*min_y;
					innerPointFounded=true;
				}
				else if(phi[X +  width * (Y-1)] == islandInnerFarValue && (Y-1 -*min_y) >= 0)
				{
					innerPoint_x = X -*min_x;
					innerPoint_y= Y-1 -*min_y;
					innerPointFounded=true;
				}
				else if(phi[X +  width * (Y+1)] == islandInnerFarValue && (Y+1 -*min_y) < island.height())
				{
					innerPoint_x = X -*min_x;
					innerPoint_y= Y+1 -*min_y;
					innerPointFounded=true;
				}
			}		
			
			//Readjust the positions
			X = X - *min_x;
			Y = Y - *min_y;

			island(X, Y, 0,0) = 0;
			island(X, Y, 0,1) = 0;
			island(X, Y, 0,2) = 0;
			
			pointNumber++;
		}

		//Take islands which have at least 5 pixels
		if(pointNumber > 5)
		{
			*allIslands = *allIslands +1;

			if(innerPointFounded)
			{
				//To fill up the islands
				fillUpIsland(phi,width, *min_x,*min_y,islandInnerFarValue,&island, innerPoint_x, innerPoint_y);
			
				//Print the inner island in a image with the size of the original image  			
				if(innerIsland)
				{	
					//Save the inner island number in a auxiliary file to make PSD file after								
					innerIslandsList << *allIslands << std::endl;
					
					*innerIslands= *innerIslands +1;

					for(int y = 0; y < island.height(); y++)
					{
						for(int x= 0; x < island.width(); x++)
						{
							if(*island.data(x, y, 0, 0) == 0)
							{
								innerIslandOnlyImage(x + *min_x, y + *min_y, 0,0) = 0;
								innerIslandOnlyImage(x + *min_x, y + *min_y, 0,1) = 0;
								innerIslandOnlyImage(x + *min_x, y + *min_y, 0,2) = 0;
							}	
						}
					}	
				}
				
				//Print the island in a image with the size of the original image and as a single image with the size of the original image			
				for(int y = 0; y < island.height(); y++)
				{
					for(int x= 0; x < island.width(); x++)
					{
						if(*island.data(x, y, 0, 0) == 0)
						{					
							singleIslandInOriginalImage(x + *min_x, y + *min_y, 0,0) = 0;
							singleIslandInOriginalImage(x + *min_x, y + *min_y, 0,1) = 0;
							singleIslandInOriginalImage(x + *min_x, y + *min_y, 0,2) = 0;
							
							allIslandsInOriginalImage(x + *min_x, y + *min_y, 0,0) = 0;
							allIslandsInOriginalImage(x + *min_x, y + *min_y, 0,1) = 0;
							allIslandsInOriginalImage(x + *min_x, y + *min_y, 0,2) = 0;
						}	
					}
				}	
			} 
			
			//Save the island as a single image
			char name[25], name2[50];
			sprintf (name, "islands/allIslands/island%d.png", *allIslands);
			island.save(name);
			sprintf (name2, "islands/allIslands_inOriginalImage/island%d.png", *allIslands);
			singleIslandInOriginalImage.save(name2);		
		}
		
        delete min_x; delete min_y; delete max_x; delete max_y;
	}//END FOR

	//Save the all inner islands as a separate image to make PSD
	innerIslandOnlyImage.save("islands/innerIslands1.png");
	innerIslandsList.close();
	
	allIslandsInOriginalImage.save("islands/allIslands1.png");
		
    delete inner_min_x; delete inner_min_y; delete inner_max_x; delete inner_max_y;
}


/** Used to fill up a island **/
void fillUpIsland(char*  phi,int width, int min_x, int min_y, int islandInnerFarValue, CImg<float>* island, int innerPoint_x, int innerPoint_y)
{
	if(*island->data(innerPoint_x, innerPoint_y, 0, 0) == 255 && phi[ ( (innerPoint_y + min_y) * width) + innerPoint_x + min_x] == islandInnerFarValue)
	{
		//Paint the pixel
		*island->data(innerPoint_x, innerPoint_y, 0,0) = 0;
		*island->data(innerPoint_x, innerPoint_y, 0,1) = 0;
		*island->data(innerPoint_x, innerPoint_y, 0,2) = 0;
		
		//Extends the paint recursively to its neighbours
		if(innerPoint_x +1 < island->width())	fillUpIsland(phi,width,min_x,min_y,islandInnerFarValue, island, innerPoint_x +1, innerPoint_y);
		if(innerPoint_x -1 >= 0)	fillUpIsland(phi,width,min_x,min_y,islandInnerFarValue,island, innerPoint_x -1, innerPoint_y);
		if(innerPoint_y +1 < island->height())	fillUpIsland(phi,width,min_x,min_y,islandInnerFarValue,island, innerPoint_x, innerPoint_y +1);
		if(innerPoint_y -1 >= 0)	fillUpIsland(phi,width,min_x,min_y,islandInnerFarValue,island, innerPoint_x, innerPoint_y -1);		
		if(innerPoint_y +1 < island->height() && innerPoint_x -1 >= 0)	fillUpIsland(phi,width,min_x,min_y,islandInnerFarValue,island, innerPoint_x -1, innerPoint_y +1);
		if(innerPoint_y +1 < island->height() && innerPoint_x +1 < island->width())	fillUpIsland(phi,width,min_x,min_y,islandInnerFarValue,island, innerPoint_x +1, innerPoint_y +1);
		if(innerPoint_y -1 >= 0 && innerPoint_x -1 >= 0)	fillUpIsland(phi,width,min_x,min_y,islandInnerFarValue,island, innerPoint_x -1, innerPoint_y -1);
		if(innerPoint_y -1 >= 0 && innerPoint_x +1 < island->width())	fillUpIsland(phi,width,min_x,min_y,islandInnerFarValue,island, innerPoint_x +1, innerPoint_y -1);		
	}
	
}


bool goOverAnIsland(char*  phi, int islandPoint, int width, int height, int* min_x, int* min_y, int* max_x, int * max_y, int *inner_min_x, int *inner_min_y, int *inner_max_x, int *inner_max_y, ofeli::list<int>* islandPoints, bool* visitedIslands , int islandInnerValue)
{
	bool touchingBorder = false;
	int x,y;
	
	y = islandPoint/width;
	x = islandPoint- (y*width);
	
	visitedIslands[x +  width * y] = true;	
	
	//Take borders to make the image later
	if(x < *min_x) *min_x = x;
	if(x > *max_x) *max_x = x;
	if(y < *min_y) *min_y = y;
	if(y > *max_y) *max_y = y;

	//Save the point
	islandPoints->push_front(islandPoint);

	//Test if the current point is in the the border of the image
	if( x == 0 || y ==0 || x == width -1 || y == height -1) touchingBorder = true;
	if(!touchingBorder)
	{
		if(x < *inner_min_x) *inner_min_x = x;
		if(x > *inner_max_x) *inner_max_x = x;
		if(y < *inner_min_y) *inner_min_y = y;
		if(y > *inner_max_y) *inner_max_y = y;
	}
   
    
    if( x-1 >= 0  )
    {
        if(phi[(x-1) +  width * y] == islandInnerValue && !visitedIslands[(x-1) +  width * y])
        {
            //Adjacent left neighbour
            if( goOverAnIsland(phi, (x-1) +  width * y , width,height, min_x, min_y, max_x, max_y, inner_min_x, inner_min_y, inner_max_x, inner_max_y, islandPoints, visitedIslands, islandInnerValue ) == true)
                touchingBorder = true;
        }	
    }
    if( x+1 < width)
    {
        if(phi[(x+1) +  width * y] == islandInnerValue && !visitedIslands[(x+1) +  width * y])
        {
            //Adjacent right neighbour
            if( goOverAnIsland(phi, (x+1) +  width * y , width,height, min_x, min_y, max_x, max_y, inner_min_x, inner_min_y, inner_max_x, inner_max_y, islandPoints, visitedIslands, islandInnerValue ) == true)
                touchingBorder = true;
        }				
    }
    if( y-1 >= 0 )
    {
        if(phi[x +  width * (y-1)] == islandInnerValue && !visitedIslands[x +  width * (y-1)])
        {
            //Adjacent bottom neighbour	
            if( goOverAnIsland(phi, x +  width * (y -1), width, height,min_x, min_y, max_x, max_y, inner_min_x, inner_min_y, inner_max_x, inner_max_y, islandPoints, visitedIslands, islandInnerValue ) == true)
                touchingBorder = true;
        }					
    }
    if(y+1 < height )
    {
        if(phi[x +  width * (y+1)] == islandInnerValue && !visitedIslands[x +  width * (y+1)])
        {
            //Adjacent top neighbour
            if( goOverAnIsland(phi, x +  width * (y +1), width,height, min_x, min_y, max_x, max_y, inner_min_x, inner_min_y, inner_max_x, inner_max_y, islandPoints, visitedIslands, islandInnerValue ) == true)
                touchingBorder = true;
        }			
    }
    if( x-1>= 0 && y+1 < height )
    {
        if(phi[(x-1) +  width * (y+1)] == islandInnerValue && !visitedIslands[(x-1) +  width * (y+1)])
        {
            if( goOverAnIsland(phi, (x-1) +  width * (y +1), width,height, min_x, min_y, max_x, max_y, inner_min_x, inner_min_y, inner_max_x, inner_max_y, islandPoints, visitedIslands, islandInnerValue ) == true)
                touchingBorder = true;
        }	
    }
    if( x-1>= 0 && y -1 >= 0 )
    {
        if(phi[(x-1) +  width * (y-1)] == islandInnerValue && !visitedIslands[(x-1) +  width * (y-1)])
        {
            if( goOverAnIsland(phi, (x -1) +  width * (y -1), width,height, min_x, min_y, max_x, max_y, inner_min_x, inner_min_y, inner_max_x, inner_max_y, islandPoints, visitedIslands, islandInnerValue ) == true)
                touchingBorder = true;
        }	
    }
    if( x+1 < width && y -1 >= 0 )
    {
        if(phi[(x+1) +  width * (y-1)] == islandInnerValue && !visitedIslands[(x+1) +  width * (y-1)])
        {
            if( goOverAnIsland(phi, (x+1) +  width * (y -1), width, height,min_x, min_y, max_x, max_y, inner_min_x, inner_min_y, inner_max_x, inner_max_y, islandPoints, visitedIslands, islandInnerValue ) == true)
                touchingBorder = true;
        }	
    }
    if( x+1 < width && y +1 < height )
    {
        if(phi[(x+1) +  width * (y+1)] == islandInnerValue && !visitedIslands[(x+1) +  width * (y+1)])
        {
            if( goOverAnIsland(phi, (x+1) +  width * (y +1), width, height,min_x, min_y, max_x, max_y, inner_min_x, inner_min_y, inner_max_x, inner_max_y, islandPoints, visitedIslands, islandInnerValue ) == true)
                touchingBorder = true;
        }	
    }

	return touchingBorder;
}


void  adjustImagesAndCreatePSDFiles(int maxWidth, int maxHeight,int originalWidth, int originalHeight, int allIslands)
{
	std::ifstream innerIslandsList("innerIslandsList.txt");
	
	//To the square images 
	if(maxWidth > maxHeight)	maxHeight = maxWidth;
	else	maxWidth = maxHeight;

	CImg<float> innerIslands("islands/innerIslands1.png"); 
	innerIslands.resize(maxWidth, maxHeight, 1, 3, 1);
	innerIslands.save("islands/scaledInnerIslands1.png");	

	CImg<float> allIslands1("islands/allIslands1.png"); 
	allIslands1.resize(maxWidth, maxHeight, 1, 3, 1);
	allIslands1.save("islands/scaledAllIslands1.png");
	
	//Island images
	for(int i=1; i <= allIslands; i++)
    {	
		char name1[50];
		sprintf (name1, "islands/allIslands/island%d.png", i);
		CImg<float> island(name1); 
		
		int posX = (maxWidth - island.width())/2;
		int posY = (maxHeight - island.height())/2;
		
		CImg<float> newIsland(maxWidth, maxHeight, 1,3,255);
		
		//Complete the new island centred and resized 
		for(int y = 0; y < island.height(); y++)
		{
			for(int x= 0; x < island.width(); x++)
			{
				if(*island.data(x, y, 0, 0) == 0)
				{
					newIsland(x+ posX, y+ posY, 0,0) = 0;
					newIsland(x+ posX, y+ posY, 0,1) = 0;
					newIsland(x+ posX, y+ posY, 0,2) = 0;
				}
				
			}
		}
		
		char name2[25];
		sprintf (name2, "islands/allIslands/island%d.png", i);
		newIsland.save(name2);	
		
		//Save the island with the size of the original image
		sprintf (name1, "islands/rescaledAllIslands/island%d.png", i);
		CImg<float> newIsland2(originalWidth, originalHeight, 1,3,255);

		int aux_x = (originalWidth - newIsland.width())/2;
		int aux_y = (originalHeight - newIsland.height())/2;
		
		for(int y = 0; y < newIsland.height(); y++)
		{
			for(int x= 0; x < newIsland.width(); x++)
			{
				if(*newIsland.data(x, y, 0, 0) == 0)
				{
					newIsland2(x + aux_x, y + aux_y, 0,0) = 0;
					newIsland2(x + aux_x, y + aux_y, 0,1) = 0;
					newIsland2(x + aux_x, y + aux_y, 0,2) = 0;
				}
				
			}
		}		
		newIsland2.save(name1);
    }
	
	int number;
	innerIslandsList >> number;
	
	//Create island PSD files 
	for(int i=1; i <= allIslands; i++)
    {
		char name1[50], name2[50];
		sprintf (name1, "islands/allIslands/island%d.png", i);
		CImg<float> island(name1); 
		
		if(number == i)
		{
			sprintf (name2, "islands/innerIslands_PSD_files/i%d.txt", i);
			innerIslandsList >> number;
		}						
		else	
		{
			sprintf (name2, "islands/allIslands/i%d.txt", i);
		}

		std::ofstream outfile(name2);

		//Set header of the PSD file	
		outfile << "#	" << island.width() << "	" << island.height() << std::endl;
						
		//Complete the new island centred and resized 
		for(int y = 0; y < island.height(); y++)
		{
			for(int x= 0; x < island.width(); x++)
			{
				if(*island.data(x, y, 0, 0) == 0)
					outfile << x << "	" << y << "		" << 0 << std::endl;
				else 
					outfile << x << "	" << y << "		" << -1 << std::endl;
			}
		}
		outfile.close();
		
		//Make another PSD file of the island with the size of the original image
		sprintf (name1, "islands/rescaledAllIslands/island%d.png", i);		
		CImg<float> island2(name1); 
		sprintf (name2, "islands/rescaledAllIslands/i%d.txt", i);			
		std::ofstream outfile2(name2);

		//Set header of the PSD file	
		outfile2 << "#	" << island2.width() << "	" << island2.height() << std::endl;
						
		//Complete the new island centred and resized 
		for(int y = 0; y < island2.height(); y++)
		{
			for(int x= 0; x < island2.width(); x++)
			{
				if(*island2.data(x, y, 0, 0) == 0)
					outfile2 << x << "	" << y << "		" << 0<< std::endl;
				else 
					outfile2 << x << "	" << y << "		" << -1 << std::endl;
			}
		}
		outfile2.close();	
		
		//Make another PSD file of the isolated island in the original image
		sprintf (name1, "islands/allIslands_inOriginalImage/island%d.png", i);		
		CImg<float> island3(name1); 
		sprintf (name2, "islands/allIslands_inOriginalImage/i%d.txt", i);			
		std::ofstream outfile3(name2);

		//Set header of the PSD file	
		outfile3 << "#	" << island3.width() << "	" << island3.height() << std::endl;
						
		//Complete the new island centred and resized 
		for(int y = 0; y < island3.height(); y++)
		{
			for(int x= 0; x < island3.width(); x++)
			{
				if(*island3.data(x, y, 0, 0) == 0)
					outfile3 << x << "	" << y << "		" << 0 << std::endl;
				else 
					outfile3 << x << "	" << y << "		" << -1 << std::endl;
			}
		}
		outfile3.close();	
    }	
	
	CImg<float> island("islands/scaledInnerIslands1.png"); 
	std::ofstream outfile("islands/scaledInnerIslands1.txt");
	
	//Set header of the PSD file	
	outfile << "#	" << island.width() << "	" << island.height() << std::endl;
									
	//Create scaled inner island image PSD file
	for(int y = 0; y < island.height(); y++)
	{
		for(int x= 0; x < island.width(); x++)
		{
			if(*island.data(x, y, 0, 0) == 0)
				outfile << x << "	" << y << "		" << 0 << std::endl;
			else 
				outfile << x << "	" << y << "		" << -1 << std::endl;
		}
	}	
	outfile.close();

	CImg<float> island2("islands/scaledAllIslands1.png"); 
	std::ofstream outfile2("islands/scaledAllIslands1.txt");
	
	//Set header of the PSD file	
	outfile2 << "#	" << island2.width() << "	" << island2.height() << std::endl;
															
	//Create scaled all islands image PSD file
	for(int y = 0; y < island2.height(); y++)
	{
		for(int x= 0; x < island2.width(); x++)
		{
			if(*island2.data(x, y, 0, 0) == 0)
				outfile2 << x << "	" << y << "		" << 0 << std::endl;
			else 
				outfile2 << x << "	" << y << "		" << -1 << std::endl;
		}
	}	
	outfile2.close();

	CImg<float> island3("islands/innerIslands1.png"); 
	std::ofstream outfile3("islands/innerIslands1.txt");
	
	//Set header of the PSD file	
	outfile3 << "#	" << island3.width() << "	" << island3.height() << std::endl;
															
	//Create inner island image PSD file
	for(int y = 0; y < island3.height(); y++)
	{
		for(int x= 0; x < island3.width(); x++)
		{
			if(*island3.data(x, y, 0, 0) == 0)
				outfile3 << x << "	" << y << "		" << 0 << std::endl;
			else 
				outfile3 << x << "	" << y << "		" << -1 << std::endl;
		}
	}	
	outfile3.close();
	
	CImg<float> island4("islands/allIslands1.png"); 
	std::ofstream outfile4("islands/allIslands1.txt");
	
	//Set header of the PSD file	
	outfile4 << "#	" << island4.width() << "	" << island4.height() << std::endl;
															
	//Create all islands image PSD file
	for(int y = 0; y < island4.height(); y++)
	{
		for(int x= 0; x < island4.width(); x++)
		{
			if(*island4.data(x, y, 0, 0) == 0)
				outfile4 << x << "	" << y << "		" << 0 << std::endl;
			else 
				outfile4 << x << "	" << y << "		" << -1 << std::endl;
		}
	}	
	outfile4.close();
	
	innerIslandsList.close();
	system("rm innerIslandsList.txt");
}
