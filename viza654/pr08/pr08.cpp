// =============================================================================
// VIZA654/CSCE646 at Texas A&M University
// Project 8
// Created by Ashwath Rajendran 
// This file is supplied with an associated makefile. Put both files in the same
// directory, navigate to that directory from the Linux shell, and type 'make'.
// This will create a program called 'pr08' that you can run by entering
// =============================================================================

#include <cstdlib>
#include <iostream>
#include <GL/glut.h>
#include<iostream>
#include<cstdio>
#include<cmath>
#include<vector>
#include<stack>
#include<queue>
#include<map>
#include<sstream>
#include<algorithm>
#include<string>
#include<limits.h>
#include <fstream>
#include <cassert>
#include <sstream>
#include <string>
#include "cyPoint.h"
#include "Color.h"


using namespace cy;
using namespace std;
#define print(p) printf("(%f,%f) \n",p.x,p.y);
#define pb(a) push_back(a)

// =============================================================================
// These variables will store the input ppm image's width, height, and color
// =============================================================================
int width, height, maxcolor=255;
unsigned char *pixmap;


int PIXEL(int x,int y) {
	return (y * width + x) * 3;
}

// =============================================================================
// OpenGL Display and Mouse Processing Functions.
//
// You can read up on OpenGL and modify these functions, as well as the commands
// in main(), to perform more sophisticated display or GUI behavior. This code
// will service the bare minimum display needs for most assignments.
// =============================================================================
static void windowResize(int w, int h)
{   
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0,(w/2),0,(h/2),0,1); 
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity() ;
}
static void windowDisplay(void)
{
  glClear(GL_COLOR_BUFFER_BIT);
  glRasterPos2i(0,0);
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, pixmap);
  glFlush();
}
static void processMouse(int button, int state, int x, int y)
{
  if(state == GLUT_UP)
  exit(0);               // Exit on mouse click.
}
static void init(void)
{
  glClearColor(1,1,1,1); // Set background color.
}

double distanceBetweenPoints(Point3f A,Point3f B) {	return sqrt((A-B)%(A-B));}
double distanceBetweenPoints(Point2f A,Point2f B) {	return sqrt((A-B)%(A-B));}

unsigned char* readPPM(const string filename,int  &width,int  &height,int  &maxcolor)
{ 
	  int ch, bit, comment;
	    FILE *fp;
	    fp=fopen(filename.c_str(),"r");				
	    //fp=fopen("red.ppm","r");			
	    unsigned char *readpixmap;		
	    if(fp == NULL)
	    {
		    printf("\n File cannot be read!! Error opening file!\n");
		    exit(0);
	    }
	    char magic[10];	
	    fscanf(fp, "%s", magic);
	    // string magic;
	    // fp>>magic;	// wrong method!!!! 
	    if(magic[0]!='P'|| (magic[1]!='6' &&  magic[1]!='5'))			
	    {
		    printf("\n Magic number %s of the input file is not P6 or P5\n", magic);
		    exit(0);
	    }
	    if(magic[1]=='6')
	    {
		    // continue with p6 algo

		    //skip all comment lines 
		    ch=fgetc(fp);	
		    do {
			    if (ch == '#');
			    ch = fgetc(fp);
		    } while (ch == '\n');
		    ungetc(ch, fp);
		    ch= fgetc(fp);
		    while (ch == '#') 
		    {
			    while (fgetc(fp) != '\n') ;
			    ch = fgetc(fp);
		    }
		    ungetc(ch, fp);
		    fscanf (fp, "%d %d %d", &width, &height, &maxcolor);	
		    fgetc(fp);
		    readpixmap = new unsigned char[width * height * 3];         
		    int y, x, pixel;
		    unsigned char red, green, blue;
		    for(y = height-1; y >= 0; y--) 
		    {
			    for(x = 0; x < width; x++) 
			    {
				    fscanf(fp, "%c%c%c", &red, &green, &blue);
				    pixel = (y * width + x) * 3; 
				    readpixmap[pixel] = red;
				    pixel++;
				    readpixmap[pixel] = green;
				    pixel++;
				    readpixmap[pixel] = blue;
			    }
		    }
	    }
	    else if(magic[1]=='5')
	    {
		    // continue with p5 algo
		    //
		    ch=fgetc(fp);	
		    do {
			    if (ch == '#');
			    ch = fgetc(fp);
		    } while (ch == '\n');
		    ungetc(ch, fp);
		    ch= fgetc(fp);
		    while (ch == '#') 
		    {
			    while (fgetc(fp) != '\n') ;
			    ch = fgetc(fp);
		    }
		    ungetc(ch, fp);
		    fscanf (fp, "%d %d %d", &width, &height, &maxcolor);	
		    fgetc(fp);
		    readpixmap = new unsigned char[width * height * 1];         
		    int y, x, pixel;
		    unsigned char red, green, blue;
		    for(y = height-1; y >= 0; y--) 
		    {
			    for(x = 0; x < width; x++) 
			    {
				    fscanf(fp, "%c", &red);
				    pixel = (y * width + x); 
				    readpixmap[pixel] = red;
				   
			    }
		    }

	    }
	    fclose(fp);
	    return readpixmap;

}			
	
void writePPM(const string filename)
{ 
	  int ch, bit, comment;
	    FILE *fp;
	    fp=fopen(filename.c_str(),"w");				
	    //fp=fopen("red.ppm","r");					
	    if(fp == NULL)
	    {
		    printf("\n File cannot be written!! Error opening file!\n");
		    exit(0);
	    }
	    char *magic = "P6";	
	    fprintf(fp, "%s", magic);
	    	    
	    fputc('\n',fp);

	    fprintf (fp, "%d %d", width, height);	
	    fputc('\n',fp);
	    fprintf (fp, "%d", maxcolor);	
	    fputc('\n',fp);
	    // pixmap = new unsigned char[width * height * 3];         
	    int y, x, pixel;
	    unsigned char red, green, blue;
	    for(y = height-1; y >= 0; y--) 
	    {
		    for(x = 0; x < width; x++) 
		    {
			    pixel = (y * width + x) * 3; 
			    red = pixmap[pixel];
			    pixel++;
			    green = pixmap[pixel];
			    pixel++;
			    blue = pixmap[pixel];
			    fprintf(fp, "%c%c%c", red, green, blue);
		    }
	    }
	    fclose(fp);

}			
	
double h0(double t){	return 2*t*t*t - 3*t*t + 1;}
double h1(double t){	return (-2)*t*t*t + 3*t*t ;}
double h2(double t){	return t*t*t - 2*t*t + t;}
double h3(double t){	return t*t*t - t*t ;}

const int filterSize = 5;

double gausFnValue(Point2f P, Point2f P0,  Point2f n0, Point2f n1, double S0, double S1)
{   return exp((-1)*( ((n0%(P-P0))*(n0%(P-P0)))/(S0*S0) + ((n1%(P-P0))*(n1%(P-P0)))/(S1*S1) ));}
unsigned char* applyFilterBlur(unsigned char *pixmap, double filter[filterSize][filterSize])
{
	int m0 = filterSize/2,n0 = m0, indx , indy, pixel;
	unsigned char *pixmapret = new unsigned char[width * height * 3];
	double sumR = 0, sumG,sumB;
	for(int y = height-1; y >= 0; y--) 
	{
		for(int x = 0; x < width; x++) 
		{
			sumR = 0, sumB=0,sumG=0;
			pixel = (y * width + x) * 3; 
			for(int m = 0; m<filterSize;m++)
			{
				for(int n = 0; n<filterSize;n++)
				{
					indx = x-n0+n; if(indx<0)	indx = 0; if(indx>=width)	indx = width - 1;
					indy = y-m0+m; if(indy<0)	indy = 0; if(indy>=height)	indy = height - 1;
					sumR += pixmap[PIXEL(indx,indy)+0]*filter[m][n];
					sumG += pixmap[PIXEL(indx,indy)+1]*filter[m][n];
					sumB += pixmap[PIXEL(indx,indy)+2]*filter[m][n];
				}
			}
			//*
			if(sumR < 0)	sumR =0;
			if(sumG < 0)	sumG =0;
			if(sumB < 0)	sumB =0;
			if(sumR > 255)	sumR = 255;
			if(sumG > 255)	sumG = 255;
			if(sumB > 255)	sumB = 255;
			//*/

			pixmapret[pixel]   = sumR; 
			pixmapret[pixel+1] = sumG;
			pixmapret[pixel+2] = sumB;
		}
	}


	return pixmapret;
}
unsigned char* applyFilterEmboss(unsigned char *pixmap, double filter[filterSize][filterSize])
{
	int m0 = filterSize/2,n0 = m0, indx , indy, pixel;
	unsigned char *pixmapret = new unsigned char[width * height * 3];
	double sumR = 0, sumG,sumB;
	for(int y = height-1; y >= 0; y--) 
	{
		for(int x = 0; x < width; x++) 
		{
			sumR = 0, sumB=0,sumG=0;
			pixel = (y * width + x) * 3; 
			for(int m = 0; m<filterSize;m++)
			{
				for(int n = 0; n<filterSize;n++)
				{
					indx = x-n0+n; if(indx<0)	indx = 0; if(indx>=width)	indx = width - 1;
					indy = y-m0+m; if(indy<0)	indy = 0; if(indy>=height)	indy = height - 1;
					sumR += pixmap[PIXEL(indx,indy)+0]*filter[m][n];
					sumG += pixmap[PIXEL(indx,indy)+1]*filter[m][n];
					sumB += pixmap[PIXEL(indx,indy)+2]*filter[m][n];
				}
			}
			//printf("sumr,g,b: %f,%f,%f \n",sumR,sumG,sumB);
			/*
			if(sumR < 0)	sumR =0;
			if(sumG < 0)	sumG =0;
			if(sumB < 0)	sumB =0;
			if(sumR > 255)	sumR = 255;
			if(sumG > 255)	sumG = 255;
			if(sumB > 255)	sumB = 255;
			*/
			sumR = (sumR+255*3)/6;
                        sumG = (sumR+255*3)/6;
                        sumB = (sumR+255*3)/6;
			if(sumR < 0)	sumR =0;
			if(sumG < 0)	sumG =0;
			if(sumB < 0)	sumB =0;
			if(sumR > 255)	sumR = 255;
			if(sumG > 255)	sumG = 255;
			if(sumB > 255)	sumB = 255;

			pixmapret[pixel]   = sumR; 
			pixmapret[pixel+1] = sumG;
			pixmapret[pixel+2] = sumB;
		}
	}


	return pixmapret;
}
unsigned char* applyFilterEdge(unsigned char *pixmap, double filter[3][3])
{
	int filterSize = 3;
	int m0 = filterSize/2,n0 = m0, indx , indy, pixel;
	unsigned char *pixmapret = new unsigned char[width * height * 3];
	double sumR = 0, sumG,sumB;
	for(int y = height-1; y >= 0; y--) 
	{
		for(int x = 0; x < width; x++) 
		{
			sumR = 0, sumB=0,sumG=0;
			pixel = (y * width + x) * 3; 
			for(int m = 0; m<filterSize;m++)
			{
				for(int n = 0; n<filterSize;n++)
				{
					indx = x-n0+n; if(indx<0)	indx = 0; if(indx>=width)	indx = width - 1;
					indy = y-m0+m; if(indy<0)	indy = 0; if(indy>=height)	indy = height - 1;
					sumR += pixmap[PIXEL(indx,indy)+0]*filter[m][n];
					sumG += pixmap[PIXEL(indx,indy)+1]*filter[m][n];
					sumB += pixmap[PIXEL(indx,indy)+2]*filter[m][n];
				}
			}
			
			sumR = (sumR+255*3)/6;
                        sumG = (sumR+255*3)/6;
                        sumB = (sumR+255*3)/6;
			if(sumR < 0)	sumR =0;
			if(sumG < 0)	sumG =0;
			if(sumB < 0)	sumB =0;
			if(sumR > 255)	sumR = 255;
			if(sumG > 255)	sumG = 255;
			if(sumB > 255)	sumB = 255;

			pixmapret[pixel]   = sumR; 
			pixmapret[pixel+1] = sumG;
			pixmapret[pixel+2] = sumB;
		}
	}


	return pixmapret;
}
double* applyFilterErode(double *pixmap, double filter[7][7])
{
	int filterSize = 7;
	int m0 = filterSize/2,n0 = m0, indx , indy, pixel;
	double*pixmapret = new double[width * height * 3];
	double sumR = 0, sumG,sumB,minColR = 500, minColG=500,minColB=500;
	for(int y = height-1; y >= 0; y--) 
	{
		for(int x = 0; x < width; x++) 
		{
			sumR = 500, sumB=500,sumG=500;
			pixel = (y * width + x) * 3; 
			for(int m = 0; m<filterSize;m++)
			{
				for(int n = 0; n<filterSize;n++)
				{
					if(filter[m][n]==0)
						continue;
					indx = x-n0+n; if(indx<0)	indx = 0; if(indx>=width)	indx = width - 1;
					indy = y-m0+m; if(indy<0)	indy = 0; if(indy>=height)	indy = height - 1;
					sumR = min(sumR,pixmap[PIXEL(indx,indy)+0]*filter[m][n]);
					sumG = min(sumG,pixmap[PIXEL(indx,indy)+1]*filter[m][n]);
					sumB = min(sumB,pixmap[PIXEL(indx,indy)+2]*filter[m][n]);
				}
			}
			
			pixmapret[pixel]   = sumR; 
			pixmapret[pixel+1] = sumG;
			pixmapret[pixel+2] = sumB;
		}
	}


	return pixmapret;
}
 double* applyFilterDilate(double *pixmap, double filter[7][7])
{
	int filterSize = 7;
	int m0 = filterSize/2,n0 = m0, indx , indy, pixel;
	double *pixmapret = new double[width * height * 3];
	double sumR = 0, sumG,sumB,minColR = 500, minColG=500,minColB=500;
	for(int y = height-1; y >= 0; y--) 
	{
		for(int x = 0; x < width; x++) 
		{
			sumR = 0, sumB=0,sumG=0;
			pixel = (y * width + x) * 3; 
			for(int m = 0; m<filterSize;m++)
			{
				for(int n = 0; n<filterSize;n++)
				{
					if(filter[m][n]==0)
						continue;
					indx = x-n0+n; if(indx<0)	indx = 0; if(indx>=width)	indx = width - 1;
					indy = y-m0+m; if(indy<0)	indy = 0; if(indy>=height)	indy = height - 1;
					sumR = max(sumR,pixmap[PIXEL(indx,indy)+0]*filter[m][n]);
					sumG = max(sumG,pixmap[PIXEL(indx,indy)+1]*filter[m][n]);
					sumB = max(sumB,pixmap[PIXEL(indx,indy)+2]*filter[m][n]);
				}
			}
			
			pixmapret[pixel]   = sumR; 
			pixmapret[pixel+1] = sumG;
			pixmapret[pixel+2] = sumB;
		}
	}


	return pixmapret;
}
void displayFilter(double filter[filterSize][filterSize])
{
	cout<<"displaying the filter: \n";
	for(int i = 0; i<filterSize;i++)
	{
		for(int j = 0; j<filterSize;j++)
		{
			printf("%f ", filter[i][j]);
		}
		cout<<endl;
	}
}
double embossFunction(Point2f p, Point2f p0 , Point2f n0)
{	return n0%(p-p0);}
/*
 *   Input RGB color primary values: r, g, and b on scale 0 - 255
 *     Output HSV colors: h on scale 0-360, s and v on scale 0-1
 *     */

#define maximum(x, y, z) ((x) > (y)? ((x) > (z)? (x) : (z)) : ((y) > (z)? (y) : (z)))
#define minimum(x, y, z) ((x) < (y)? ((x) < (z)? (x) : (z)) : ((y) < (z)? (y) : (z)))

void RGBtoHSV(int r, int g, int b, double &h, double &s, double &v)
{
	double red, green, blue;
	double max, min, delta;

	red = r / 255.0; green = g / 255.0; blue = b / 255.0;  /* r, g, b to 0 - 1 scale */

	max = maximum(red, green, blue);
	min = minimum(red, green, blue);

	v = max;        /* value is maximum of r, g, b */

	if(max == 0){    /* saturation and hue 0 if value is 0 */
		s = 0;
		h = 0;
	}
	else
	{
		s = (max - min) / max;           /* saturation is color purity on scale 0 - 1 */

		delta = max - min;
		if(delta == 0)                    /* hue doesn't matter if saturation is 0 */
			h = 0;
		else{
			if(red == max)                  /* otherwise, determine hue on scale 0 - 360 */
				h = (green - blue) / delta;
			else if(green == max)
				h = 2.0 + (blue - red) / delta;
			else /* (blue == max) */
				h = 4.0 + (red - green) / delta;
			h = h * 60.0;
			if(h < 0)
				h = h + 360.0;
		}
	}
}

unsigned char* ImageRGBtoHSV(unsigned char* rgbImage)
{
	unsigned char* hsvImage = new unsigned char[width * height * 3 ];
	for(int y = height-1; y >= 0; y--) 
	{
		for(int x = 0; x < width; x++) 
		{
			int pixel = (y * width + x) * 3; 
			double h=0,s=0,v=0;
			RGBtoHSV(rgbImage[pixel], rgbImage[pixel+1],rgbImage[pixel+2],h,s,v);
			
			//if(rgbImage[pixel+1]>200 && rgbImage[pixel+2]<200)
			//	printf("r,g,b,h,s,v,%d,%d,%d, %f,%f,%f\n",rgbImage[pixel], rgbImage[pixel+1],rgbImage[pixel+2],h,s,v);
			hsvImage[pixel]   = h;
			hsvImage[pixel+1] = s;
			hsvImage[pixel+2] = v;
		}
	}
	return hsvImage;

}
double* generateAlphaMap(unsigned char* hsvImage,double minHue, double maxHue,double minSaturation, double maxSaturation, double minValue, double maxValue)
{	
	// given an image returns the alphamap
	double* alphaMap = new double[width * height * 3];
	for(int y = height-1; y >= 0; y--) 
	{
		for(int x = 0; x < width; x++) 
		{
			int pixel = (y * width + x) * 3; 
				// printf("h,s,v,%d,%d,%d, \n",hsvImage[pixel], hsvImage[pixel+1],hsvImage[pixel+2]);
			if(hsvImage[pixel]<=maxHue && hsvImage[pixel]>=minHue && hsvImage[pixel+1]<=maxSaturation && hsvImage[pixel+1]>=minSaturation && hsvImage[pixel+2]<=maxValue && hsvImage[pixel+2]>=minValue )
			{
				// printf("h,s,v,%d,%d,%d, \n",hsvImage[pixel], hsvImage[pixel+1],hsvImage[pixel+2]);
				alphaMap[pixel] = 1;
			}
			else 
			{
				alphaMap[pixel] = 0;
				// cout<<(int)alphaMap[pixel];
			}
		}
	}
	return alphaMap;
}
double trans[3][3];
double transInverse[3][3]; 
void copyMatToTrans(double a[3][3])
{

	for (int row = 0; row < 3; row++) {
		for (int col = 0; col < 3; col++) {
			trans[row][col] = a[row][col];
		}
	}
           
}
void copyMatToTransInv(double a[3][3])
{

	for (int row = 0; row < 3; row++) {
		for (int col = 0; col < 3; col++) {
			transInverse[row][col] = a[row][col];
		}
	}
           
}
Point2f transPoint(double x, double y)
{
    // assuming matA = mxn; matB = nxp
    // multiply matrices 
    double product[3][1]={0}; 
    double bMatrix[3][1];
    bMatrix[0][0]=x;
    bMatrix[1][0]=y;
    bMatrix[2][0]=1;

        for (int row = 0; row < 3; row++) {
            //for (int col = 0; col < 3; col++) {
            int col=0;
            // Multiply the row of A by the column of B to get the row, column of product.
                for (int inner = 0; inner < 3; inner++) {
                    product[row][col] += trans[row][inner] * bMatrix[inner][col];
                }
                //std::cout << product[row][col] << "  ";
            //}
            //std::cout << "\n";
        }
    return Point2f(product[0][0],product[1][0]);
}

Point2f transInvPoint(double x, double y)
{
    // assuming matA = mxn; matB = nxp
    // multiply matrices 
    double product[3][1]={0}; 
    double bMatrix[3][1];
    bMatrix[0][0]=x;
    bMatrix[1][0]=y;
    bMatrix[2][0]=1;

        for (int row = 0; row < 3; row++) {
            //for (int col = 0; col < 3; col++) {
            int col=0;
            // Multiply the row of A by the column of B to get the row, column of product.
                for (int inner = 0; inner < 3; inner++) {
                    product[row][col] += transInverse[row][inner] * bMatrix[inner][col];
                }
                //std::cout << product[row][col] << "  ";
            //}
            //std::cout << "\n";
        }
    return Point2f(product[0][0],product[1][0]);
}
int numToReduce = 8;
	
// =============================================================================
// main() Program Entry
// =============================================================================
int main(int argc, char *argv[]) 
{

	// Check the number of parameters
	if (argc < 2) {
		// Tell the user how to run the program
		std::cout<< "Usage: " << argv[0] << " option" << std::endl;
		return 1;
	}
	
	int widthFore=0,widthBack=0,heightFore=0,heightBack=0;
	string input = argv[1];
	int N = 1, M = N;
	float I = 0, J = 0;
	cout<< "input arg is : " <<input<<endl;
	//initialize the global variables
	width = 640;
	height = 480;
	pixmap = new unsigned char[width * height * 3]; 
	if(input == "floyd")
	{
		if (argc < 2) {
			std::cout<< "Usage: " << argv[0] << " floyd <image> [numToReduce]"  << std::endl;
			return 1;
		}


		string Image = argv[2];
		if(argc == 4)	numToReduce=atof(argv[3]);

		unsigned char* pixmapInput = readPPM(Image,widthFore,heightFore,maxcolor);
		width = widthFore; height = heightFore;
		unsigned char* pixmapRet = new unsigned char[width * height * 3]; 
		printf("W %d, H %d\n", width ,height);
		int level = numToReduce;
		int step = 255 / level;
		// cout<<step<<endl<<round(200*1.0/255*3)*step<<endl;
		for(int y=height-1;y>=0;y--){
			for(int x=0;x<width;x++){
				
				int pixel = ( y * width + x) * 3; 
				pixmapRet[pixel]  = round(pixmapInput[pixel]*1.0/255*level)*step   ;
				pixmapRet[pixel+1]= round(pixmapInput[pixel+1]*1.0/255*level)*step ;
				pixmapRet[pixel+2]= round(pixmapInput[pixel+2]*1.0/255*level)*step ;
				double errRed  = abs(pixmapInput[pixel] - pixmapRet[pixel]);
				double errBlue = abs(pixmapInput[pixel+1] - pixmapRet[pixel+1]);
				double errGreen= abs(pixmapInput[pixel+2] - pixmapRet[pixel+2]);


				//printf("W %d, H %d\n",PIXEL(x+1,y),height);
				//*
				if(x+1 < width )		
				{
					pixmapInput[PIXEL(x+1,y)]  =min(pixmapInput[PIXEL(x+1,y)]+(int)(errRed*7/16.0),255);
					pixmapInput[PIXEL(x+1,y)+1]=min(pixmapInput[PIXEL(x+1,y)+1]+(int)(errGreen*7/16.0),255);
					pixmapInput[PIXEL(x+1,y)+2]=min(pixmapInput[PIXEL(x+1,y)+2]+(int)(errBlue*7/16.0),255);
				}
				if(y-1 >= 0)			
				{
					pixmapInput[PIXEL(x,y-1)]  =min(pixmapInput[PIXEL(x,y-1)]+(int)(errRed*5/16.0),255);    
					pixmapInput[PIXEL(x,y-1)+1]=min(pixmapInput[PIXEL(x,y-1)+1]+(int)(errGreen*5/16.0),255);
					pixmapInput[PIXEL(x,y-1)+2]=min(pixmapInput[PIXEL(x,y-1)+2]+(int)(errBlue*5/16.0),255); 
				}
				if(x+1 < width && y-1>=0 )	
				{
					pixmapInput[PIXEL(x+1,y-1)]  =min(pixmapInput[PIXEL(x+1,y-1)]+(int)(errRed/16.0),255);    
					pixmapInput[PIXEL(x+1,y-1)+1]=min(pixmapInput[PIXEL(x+1,y-1)+1]+(int)(errGreen/16.0),255);
					pixmapInput[PIXEL(x+1,y-1)+2]=min(pixmapInput[PIXEL(x+1,y-1)+2]+(int)(errBlue/16.0),255); 
				}
				if(x-1 >=0 && y-1>=0 )	
				{
					pixmapInput[PIXEL(x-1,y-1)]  =min(pixmapInput[PIXEL(x-1,y-1)]+(int)(errRed*3/16.0),255);    
					pixmapInput[PIXEL(x-1,y-1)+1]=min(pixmapInput[PIXEL(x-1,y-1)+1]+(int)(errGreen*3/16.0),255);
					pixmapInput[PIXEL(x-1,y-1)+2]=min(pixmapInput[PIXEL(x-1,y-1)+2]+(int)(errBlue*3/16.0),255); 
				}

				//*/

				if(pixmapRet[pixel]<0 || pixmapRet[pixel+1]<0 || pixmapRet[pixel+2]<0||pixmapRet[pixel]>255 || pixmapRet[pixel+1]>255 || pixmapRet[pixel+2]>255  )
					cout<<"WHAT!!!";
				pixmapRet[pixel]  = min((int)pixmapRet[pixel], 255)  ;
				pixmapRet[pixel+1]= min((int)pixmapRet[pixel+1],255) ;
				pixmapRet[pixel+2]= min((int)pixmapRet[pixel+2],255) ;

			}
		}
		printf("W %d, H %d\n", width ,height);

		width = widthFore;
		height = heightFore;

		pixmap = pixmapRet;
	//exit(0);
		
	}
	else if(input == "random")
	{
		if (argc < 2) {
			std::cout<< "Usage: " << argv[0] << " random <image> [numToReduce]"  << std::endl;
			return 1;
		}


		string Image = argv[2];
		if(argc == 4)	numToReduce=atof(argv[3]);

		unsigned char* pixmapInput = readPPM(Image,widthFore,heightFore,maxcolor);
		width = widthFore; height = heightFore;
		unsigned char* pixmapRet = new unsigned char[width * height * 3]; 
		printf("W %d, H %d\n", width ,height);
		int level = numToReduce;
		int step = 255 / level;
		// cout<<step<<endl<<round(200*1.0/255*3)*step<<endl;
		for(int y=height-1;y>=0;y--){
			for(int x=0;x<width;x++){
				
				int pixel = ( y * width + x) * 3; 
				pixmapInput[pixel]  +=  ((float)((float)(rand()%100)/(100.0)) - 0.5)*step;
				pixmapInput[pixel+1]+=  ((float)((float)(rand()%100)/(100.0)) - 0.5)*step;		
				pixmapInput[pixel+2]+=  ((float)((float)(rand()%100)/(100.0)) - 0.5)*step;
				pixmapRet[pixel]  = round(pixmapInput[pixel]  *1.0/255*level)*step   ;
				pixmapRet[pixel+1]= round(pixmapInput[pixel+1]*1.0/255*level)*step ;
				pixmapRet[pixel+2]= round(pixmapInput[pixel+2]*1.0/255*level)*step ;
				double errRed  = abs(pixmapInput[pixel] - pixmapRet[pixel]);
				double errBlue = abs(pixmapInput[pixel+1] - pixmapRet[pixel+1]);
				double errGreen= abs(pixmapInput[pixel+2] - pixmapRet[pixel+2]);


				//printf("W %d, H %d\n",PIXEL(x+1,y),height);
				//*
				if(x+1 < width )		
				{
					errRed = 
					pixmapInput[PIXEL(x+1,y)]  =min(pixmapInput[PIXEL(x+1,y)]+(int)(errRed/2.0),255);
					pixmapInput[PIXEL(x+1,y)+1]=min(pixmapInput[PIXEL(x+1,y)+1]+(int)(errGreen/2.0),255);
					pixmapInput[PIXEL(x+1,y)+2]=min(pixmapInput[PIXEL(x+1,y)+2]+(int)(errBlue/2.0),255);
				}
				if(y-1 >= 0)			
				{
					pixmapInput[PIXEL(x,y-1)]  =min(pixmapInput[PIXEL(x,y-1)]+(int)(errRed/4.0),255);    
					pixmapInput[PIXEL(x,y-1)+1]=min(pixmapInput[PIXEL(x,y-1)+1]+(int)(errGreen/4.0),255);
					pixmapInput[PIXEL(x,y-1)+2]=min(pixmapInput[PIXEL(x,y-1)+2]+(int)(errBlue/4.0),255); 
				}
				if(x+1 < width && y-1>=0 )	
				{
					pixmapInput[PIXEL(x+1,y-1)]  =min(pixmapInput[PIXEL(x+1,y-1)]+(int)(errRed/4.0),255);    
					pixmapInput[PIXEL(x+1,y-1)+1]=min(pixmapInput[PIXEL(x+1,y-1)+1]+(int)(errGreen/4.0),255);
					pixmapInput[PIXEL(x+1,y-1)+2]=min(pixmapInput[PIXEL(x+1,y-1)+2]+(int)(errBlue/4.0),255); 
				}
				//*/

				if(pixmapRet[pixel]<0 || pixmapRet[pixel+1]<0 || pixmapRet[pixel+2]<0||pixmapRet[pixel]>255 || pixmapRet[pixel+1]>255 || pixmapRet[pixel+2]>255  )
					cout<<"WHAT!!!";
				pixmapRet[pixel]  = min((int)pixmapRet[pixel], 255)  ;
				pixmapRet[pixel+1]= min((int)pixmapRet[pixel+1],255) ;
				pixmapRet[pixel+2]= min((int)pixmapRet[pixel+2],255) ;

			}
		}
		printf("W %d, H %d\n", width ,height);

		width = widthFore;
		height = heightFore;

		pixmap = pixmapRet;
	//exit(0);
		
	}

	else if(input == "uniform")
	{
		if (argc < 2) {
			std::cout<< "Usage: " << argv[0] << " uniform <image> [numToReduce]"  << std::endl;
			return 1;
		}


		string Image = argv[2];
		if(argc == 4)	numToReduce=atof(argv[3]);

		unsigned char* pixmapInput = readPPM(Image,widthFore,heightFore,maxcolor);
		width = widthFore; height = heightFore;
		unsigned char* pixmapRet = new unsigned char[width * height * 3]; 
		printf("W %d, H %d\n", width ,height);
		int level = numToReduce;
		int step = 255 / level;
		// cout<<step<<endl<<round(200*1.0/255*3)*step<<endl;
		for(int y=height-1;y>=0;y--){
			for(int x=0;x<width;x++){
				
				int pixel = ( y * width + x) * 3; 
				pixmapRet[pixel]  = round(pixmapInput[pixel]*1.0/255*level)*step   ;
				pixmapRet[pixel+1]= round(pixmapInput[pixel+1]*1.0/255*level)*step ;
				pixmapRet[pixel+2]= round(pixmapInput[pixel+2]*1.0/255*level)*step ;
				double errRed  = abs(pixmapInput[pixel] - pixmapRet[pixel]);
				double errBlue = abs(pixmapInput[pixel+1] - pixmapRet[pixel+1]);
				double errGreen= abs(pixmapInput[pixel+2] - pixmapRet[pixel+2]);


				//printf("W %d, H %d\n",PIXEL(x+1,y),height);
				//*
				if(x+1 < width )		
				{
					pixmapInput[PIXEL(x+1,y)]  =min(pixmapInput[PIXEL(x+1,y)]+(int)(errRed/2.0),255);
					pixmapInput[PIXEL(x+1,y)+1]=min(pixmapInput[PIXEL(x+1,y)+1]+(int)(errGreen/2.0),255);
					pixmapInput[PIXEL(x+1,y)+2]=min(pixmapInput[PIXEL(x+1,y)+2]+(int)(errBlue/2.0),255);
				}
				if(y-1 >= 0)			
				{
					pixmapInput[PIXEL(x,y-1)]  =min(pixmapInput[PIXEL(x,y-1)]+(int)(errRed/4.0),255);    
					pixmapInput[PIXEL(x,y-1)+1]=min(pixmapInput[PIXEL(x,y-1)+1]+(int)(errGreen/4.0),255);
					pixmapInput[PIXEL(x,y-1)+2]=min(pixmapInput[PIXEL(x,y-1)+2]+(int)(errBlue/4.0),255); 
				}
				if(x+1 < width && y-1>=0 )	
				{
					pixmapInput[PIXEL(x+1,y-1)]  =min(pixmapInput[PIXEL(x+1,y-1)]+(int)(errRed/4.0),255);    
					pixmapInput[PIXEL(x+1,y-1)+1]=min(pixmapInput[PIXEL(x+1,y-1)+1]+(int)(errGreen/4.0),255);
					pixmapInput[PIXEL(x+1,y-1)+2]=min(pixmapInput[PIXEL(x+1,y-1)+2]+(int)(errBlue/4.0),255); 
				}
				//*/

				if(pixmapRet[pixel]<0 || pixmapRet[pixel+1]<0 || pixmapRet[pixel+2]<0||pixmapRet[pixel]>255 || pixmapRet[pixel+1]>255 || pixmapRet[pixel+2]>255  )
					cout<<"WHAT!!!";
				pixmapRet[pixel]  = min((int)pixmapRet[pixel], 255)  ;
				pixmapRet[pixel+1]= min((int)pixmapRet[pixel+1],255) ;
				pixmapRet[pixel+2]= min((int)pixmapRet[pixel+2],255) ;

			}
		}
		printf("W %d, H %d\n", width ,height);

		width = widthFore;
		height = heightFore;

		pixmap = pixmapRet;
	//exit(0);
		
	}
	else
	{
		
		cout<<"invalid argument\n";
		return 1;
	}

	// OpenGL Commands:
	// Once "glutMainLoop" is executed, the program loops indefinitely to all
	// glut functions.  
	glutInit(&argc, argv);
	glutInitWindowPosition(100, 100); // Where the window will display on-screen.
	glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutCreateWindow("pr08");
	init();
	glutReshapeFunc(windowResize);
	glutDisplayFunc(windowDisplay);
	glutMouseFunc(processMouse);
	glutMainLoop();

	return 0; //This line never gets reached. 
}

