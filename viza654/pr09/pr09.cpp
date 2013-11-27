// =============================================================================
// VIZA654/CSCE646 at Texas A&M University
// Project 9
// Created by Ashwath Rajendran 
// This file is supplied with an associated makefile. Put both files in the same
// directory, navigate to that directory from the Linux shell, and type 'make'.
// This will create a program called 'pr09' that you can run by entering
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
#define PIXEL(x,y) (y * width + x) * 3

// =============================================================================
// These variables will store the input ppm image's width, height, and color
// =============================================================================
int width, height, maxcolor=255;
unsigned char *pixmap;
string format = "P6";
double tx=0, ty = 0, sx=0,sy=0;



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

	    if(magic[1]=='6')
	    {
		    format = "P6";

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
		    format = "P5";

		    for(y = height-1; y >= 0; y--) 
		    {
			    for(x = 0; x < width; x++) 
			    {
				    fscanf(fp, "%c", &red);
				    pixel = (y * width + x) * 3; 
				    readpixmap[pixel] = red;
				    pixel++;
				    readpixmap[pixel] = red;
				    pixel++;
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
	    char *magic;	
	    if(format == "P6")
	    {
		    magic = "P6";
	    }else if(format == "P5")
	    {
		    magic = "P5";
	    }


	    fprintf(fp, "%s", magic);
	    	    
	    fputc('\n',fp);

	    fprintf (fp, "%d %d", width, height);	
	    fputc('\n',fp);
	    fprintf (fp, "%d", maxcolor);	
	    fputc('\n',fp);
	    // pixmap = new unsigned char[width * height * 3];         
	    int y, x, pixel;
	    unsigned char red, green, blue;
	    if(format == "P6")
	    {
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
	    }
	    else if(format == "P5")
	    {
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
				    fprintf(fp, "%c", red);
			    }
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
			
			//if(rgbImage[pixel]>10)
			//    printf("r,g,b,h,s,v,%d,%d,%d, %f,%f,%f\n",rgbImage[pixel], rgbImage[pixel+1],rgbImage[pixel+2],h,s,v);
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
				//printf("h,s,v,%d,%d,%d, \n",hsvImage[pixel], hsvImage[pixel+1],hsvImage[pixel+2]);
				alphaMap[pixel] = 1;
			}
			else 
			{
				alphaMap[pixel] = 0;
				//printf("alpha is zero here : h,s,v,%d,%d,%d, \n",hsvImage[pixel], hsvImage[pixel+1],hsvImage[pixel+2]);
				//cout<<(int)alphaMap[pixel];
			}
		}
	}
	return alphaMap;
}

double* correctAlphaMap(unsigned char* hsvImage,double minHue, double maxHue,double minSaturation, double maxSaturation, double minValue, double maxValue, double* alpha)
{	
	// given an image returns the alphamap
	double* alphaMap = new double[width * height * 3];
	for(int y = height-1; y >= 0; y--) 
	{
		for(int x = 0; x < width; x++) 
		{
			int pixel = (y * width + x) * 3; 
				// printf("h,s,v,%d,%d,%d, \n",hsvImage[pixel], hsvImage[pixel+1],hsvImage[pixel+2]);
			if(alpha[pixel]==1 && hsvImage[pixel]<=maxHue && hsvImage[pixel]>=minHue && hsvImage[pixel+1]<=maxSaturation && hsvImage[pixel+1]>=minSaturation && hsvImage[pixel+2]<=maxValue && hsvImage[pixel+2]>=minValue )
			{
				//printf("h,s,v,%d,%d,%d, \n",hsvImage[pixel], hsvImage[pixel+1],hsvImage[pixel+2]);
				alphaMap[pixel] = 1;
			}
			else if(alpha[pixel]==0 && hsvImage[pixel]<=maxHue && hsvImage[pixel]>=minHue && hsvImage[pixel+1]<=maxSaturation && hsvImage[pixel+1]>=minSaturation && hsvImage[pixel+2]<=maxValue && hsvImage[pixel+2]>=minValue )
			{
				//printf("h,s,v,%d,%d,%d, \n",hsvImage[pixel], hsvImage[pixel+1],hsvImage[pixel+2]);
				alphaMap[pixel] = 1;
			}

			else 
			{
				alphaMap[pixel] = alpha[pixel];
				//printf("alpha is zero here : h,s,v,%d,%d,%d, \n",hsvImage[pixel], hsvImage[pixel+1],hsvImage[pixel+2]);
				//cout<<(int)alphaMap[pixel];
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


unsigned char* warpOperation(int argc, char *argv[], string input, float T,unsigned char* pixmapToComposite)
{
	T = T/2.0;
int widthFore=0,widthBack=0,heightFore=0,heightBack=0;
	int N = 1, M = N;
	float I = 0, J = 0;
	//cout<< "input arg is : " <<input<<endl;
	//initialize the global variables
	//width = 120;
	//height = 160;
	string Img = argv[2];
	pixmap = readPPM(Img,width,height,maxcolor);
	unsigned char *pixmapRet = new unsigned char[width * height * 3]; 
	if(input == "over")
	{
		if (argc < 9) {
			std::cout<< "Usage: " << argv[0] << " over <background image> <foreground or masking image2> min_hue max_hue min_saturation max_saturation min_value max_value"  << std::endl;
			return pixmapRet;
		}


		string backgroundImage = argv[2],foregroundImage = argv[3];
		double min_hue=atof(argv[4]), max_hue =atof(argv[5]), min_saturation=atof(argv[6]),max_saturation=atof(argv[7]), min_value=atof(argv[8]), max_value=atof(argv[9]);

		unsigned char* pixmapFore = readPPM(foregroundImage,widthFore,heightFore,maxcolor);
		pixmapFore = pixmapToComposite;


		unsigned char* pixmapBack = readPPM(backgroundImage,widthBack,heightBack,maxcolor);

		unsigned char* hsvImage = ImageRGBtoHSV(pixmapFore);
		cout<<"\nhsv image created..\n";

		double* alphaMap = generateAlphaMap(hsvImage,min_hue,max_hue,min_saturation,max_saturation,min_value,max_value);  
		double filter[7][7]=	{	
			{0,0,0,1,0,0,0},
			{0,0,1,1,1,0,0},
			{0,1,1,1,1,1,0},
			{1,1,1,1,1,1,1},
			{0,1,1,1,1,1,0},
			{0,0,1,1,1,0,0},
			{0,0,0,1,0,0,0}, };
		alphaMap = applyFilterErode(alphaMap, filter);
		alphaMap = applyFilterDilate(alphaMap, filter);

		//alphaMap = correctAlphaMap(hsvImage,min_hue,max_hue,min_saturation,max_saturation,min_value,max_value, alphaMap); 

		printf("mH= %f,mxH= %f,mS= %f,mxS= %f,mV= %f,mxV= %f\n",min_hue,max_hue,min_saturation,max_saturation,min_value,max_value);
		if(widthFore != widthBack || heightFore != heightBack)
		{
			printf("height or width mismatch between fore and background images\n wF %d,hF %d,wB %d,hB %d\n",widthFore,widthBack,heightFore,heightBack);
			return pixmapRet;
		}
		width = widthFore;height = heightFore;

		for(int y = height-1; y >= 0; y--) 
		{
			for(int x = 0; x < width; x++) 
			{
				int pixel = (y * width + x) * 3; 
				// set the following line for normal operation in photoshop
				// alphaMap[pixel]=0.5;
				//printf("alphaMap of pixel = %f\n",alphaMap[pixel]);
				pixmapRet[pixel]   = alphaMap[pixel]*pixmapBack[pixel]   + (1-alphaMap[pixel])*pixmapFore[pixel]; 
				pixmapRet[pixel+1] = alphaMap[pixel]*pixmapBack[pixel+1] + (1-alphaMap[pixel])*pixmapFore[pixel+1];
				pixmapRet[pixel+2] = alphaMap[pixel]*pixmapBack[pixel+2] + (1-alphaMap[pixel])*pixmapFore[pixel+2];
			}
		}


	}
	else if(input == "translate")
	{
		if (argc < 3) {
			std::cout<< "Usage: " << argv[0] << " scale <image> <tx> <ty>"  << std::endl;
			return pixmapRet;
		}


		string Image = argv[3];//argv[2];
		tx=+T*10;
		ty=-T*20;
		cout<<"scale x,y: "<<tx<<", "<<ty<<endl;

		unsigned char* pixmapInput = readPPM(Image,widthFore,heightFore,maxcolor);

		width = widthFore;
		height = heightFore;

		printf("W %d, H %d\n", width ,height);
		double trans[3][3]={
			{1 , 0, tx},
			{0 , 1, ty},
			{0 , 0,  1}
		};
		double transInverse[3][3] ={
			{1 , 0, -tx},
			{0 , 1, -ty},
			{0 , 0,  1}
		}; 
		copyMatToTrans(trans);
		copyMatToTransInv(transInverse);

        // calculating new width and height
        int maxWidth = INT_MIN,minWidth=INT_MAX,maxHeight=INT_MIN,minHeight=INT_MAX;
        /*
	for(int y=height-1;y>=0;y--){
			for(int x=0;x<width;x++){
				Point2f point=transPoint(x,y);
			if (point.x > maxWidth) maxWidth = point.x;
      			if (point.y > maxHeight) maxHeight = point.y;
      			if (point.x < minWidth) minWidth = point.x;
      			if (point.y < minHeight) minHeight = point.y;
      			
      		}
      	}
	*/
	printf("maxW %d, minW %d, maxH %d, minH %d\n", maxWidth ,minWidth,maxHeight,minHeight);
	// int offsetWidth = abs(minWidth), offsetHeight = abs(minHeight);
	maxWidth = width ;// + abs(tx);
	minWidth = 0;
	maxHeight = height ;//+ abs(ty);
	minHeight = 0;
	printf("maxW %d, minW %d, maxH %d, minH %d\n", maxWidth ,minWidth,maxHeight,minHeight);
	// pixmapRet = new unsigned char[maxWidth*maxHeight*3];

	// Actual translation
	for(int y=maxHeight-1;y>=0;y--){
		for(int x=0;x<maxWidth;x++){
				Point2f point=transInvPoint(x,y);
				int pixel = ( (y) * maxWidth + (x)) * 3; 
				if(point.x >= width || point.x < 0 || point.y >= height || point.y < 0)
				{
					// set color as black
					pixmapRet[pixel]   = 0; 
					pixmapRet[pixel+1] = 0;
					pixmapRet[pixel+2] = 0;
				}
				else
				{
					int pixelInput = ((int)point.y * width + (int)point.x)*3;
					pixmapRet[pixel]   = pixmapInput[pixelInput];
                                        pixmapRet[pixel+1] = pixmapInput[pixelInput+1];
				        pixmapRet[pixel+2] = pixmapInput[pixelInput+2];
				}
		}
	}

	width = maxWidth;
	height= maxHeight;

	//exit(0);
		
	}
	else if(input == "scale")
	{
		if (argc < 5) {
			std::cout<< "Usage: " << argv[0] << " scale <image> <sx> <sy>"  << std::endl;
			return NULL;
		}


		string Image = argv[3];
		//double sx;
		//double sy;
		sx=1+T*0.15;
		sy=1+T*0.15;

		cout<<"scale x,y: "<<sx<<", "<<sy<<endl;

		unsigned char* pixmapInput = readPPM(Image,widthFore,heightFore,maxcolor);
		pixmapInput = pixmapToComposite;
		width = widthFore;
		height = heightFore;

		printf("W %d, H %d\n", width ,height);
		double trans[3][3]={
			{sx, 0, 0},
			{0 ,sy, 0},
			{0 , 0, 1}
		};
		double transInverse[3][3] ={
			{1.0/sx, 0, 0},
			{0 ,1.0/sy, 0},
			{0 , 0, 1}
		}; 
		copyMatToTrans(trans);
		copyMatToTransInv(transInverse);

        // calculating new width and height
        int maxWidth = INT_MIN,minWidth=INT_MAX,maxHeight=INT_MIN,minHeight=INT_MAX;
        /*
	 * for(int y=height-1;y>=0;y--){
			for(int x=0;x<width;x++){
				Point2f point=transPoint(x,y);
			if (point.x > maxWidth) maxWidth = point.x;
      			if (point.y > maxHeight) maxHeight = point.y;
      			if (point.x < minWidth) minWidth = point.x;
      			if (point.y < minHeight) minHeight = point.y;
      			
      		}
      	}
	*/
	printf("maxW %d, minW %d, maxH %d, minH %d\n", maxWidth ,minWidth,maxHeight,minHeight);
	int offsetWidth = abs(minWidth), offsetHeight = abs(minHeight);
	maxWidth = (int)(width);
	minWidth = 0;
	maxHeight = (int)(height);
	minHeight = 0;
	printf("maxW %d, minW %d, maxH %d, minH %d\n", maxWidth ,minWidth,maxHeight,minHeight);
	// pixmapRet = new unsigned char[maxWidth*maxHeight*3];

	// Actual scaling
	for(int y=maxHeight-1;y>=0;y--){
		for(int x=0;x<maxWidth;x++){
				Point2f point=transInvPoint(x,y);
				int pixel = ( (y) * maxWidth + (x)) * 3; 
				if(point.x >= width || point.x < 0 || point.y >= height || point.y < 0)
				{
					// set color as black
					pixmapRet[pixel]   = 0; 
					pixmapRet[pixel+1] = 0;
					pixmapRet[pixel+2] = 0;
				}
				else
				{
					int pixelInput = ((int)point.y * width + (int)point.x)*3;
					pixmapRet[pixel]   = pixmapInput[pixelInput];
                                        pixmapRet[pixel+1] = pixmapInput[pixelInput+1];
				        pixmapRet[pixel+2] = pixmapInput[pixelInput+2];
				}
		}
	}

	width = maxWidth;
	height= maxHeight;

	//exit(0);
		
	}
	else
	{
		
		cout<<"invalid argument\n";
		return pixmapRet;
	}
	return pixmapRet;
}


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
	int numTime = 30;
	string input = argv[1],outputFileName="output",ext = ".ppm";
	unsigned char *pixmap1,*pixmapPrev;

	//width = 1200;
	//height = 1600;
	for(int T = 0; T<numTime;T++)
	{
	    cout<<"start\n"<<argc<<input<<endl;
		input = "translate";
		pixmap1 = warpOperation(argc,argv,input, T,NULL);
		input = "scale";
		pixmap1 = warpOperation(argc,argv,input, T,pixmap1);
		

		cout<<"\n NEXT....\n";
		input = "over";
		pixmap = warpOperation(argc,argv,input, T,pixmap1);
		cout<<"end\n";

		if(pixmap == NULL)
		{
			cout<<"null detected\n";
			continue;
		}
		
		pixmapPrev = pixmap;

		// Write the output file that is being displayed
		//*
		time_t newTime;
		time(&newTime);

		std::string number;std::stringstream strstream;strstream <<(newTime+T);strstream >> number;
		string fileName = outputFileName + number + ext;cout<<fileName<<endl;

		if(format == "P6")	writePPM(fileName);
		else if(format == "P5")	writePPM(fileName);
		// break;
		
		// */
		// pixmap = NULL;

	}
	// OpenGL Commands:
	// Once "glutMainLoop" is executed, the program loops indefinitely to all
	// glut functions.  
	glutInit(&argc, argv);
	glutInitWindowPosition(100, 100); // Where the window will display on-screen.
	glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutCreateWindow("pr09");
	init();
	glutReshapeFunc(windowResize);
	glutDisplayFunc(windowDisplay);
	glutMouseFunc(processMouse);
	glutMainLoop();

	return 0; //This line never gets reached. 
}

