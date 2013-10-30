// =============================================================================
// VIZA654/CSCE646 at Texas A&M University
// Project 1
// Created by Ashwath Rajendran 
// This file is supplied with an associated makefile. Put both files in the same
// directory, navigate to that directory from the Linux shell, and type 'make'.
// This will create a program called 'pr01' that you can run by entering
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

unsigned char* readPPM(const string filename)
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
	
vector<Point3f> swapPoints(vector<Point3f> v, int src, int dest )
{
	int N = v.size();
	if(src < 0 || dest < 0 || src >= N || dest >=N)
	{
		cout<<"Entered src dest is larger than maximum number of elements in vector\n";
		exit(0);
	}

	Point3f t = v[src];
	v[src] = v[dest];
	v[dest] = t;
	return v;
}
void UniformClusterImage(const string inputImage,int k)
{}
void KmeansClusterImage(const string inputImage,int k)
{

	if(k>maxcolor)
	{
		cout<<"Entered number is larger than the maximum number of colors in the image\n";
		exit(0);
	}
	unsigned char *image = readPPM(inputImage.c_str());
	pixmap = image;
	vector<vector<Point3f> > AllClusters;
	AllClusters.resize(k);
    int minInd = 0, maxInd = 1;
    double thresh = 5.0;
    int rndIndex = ceil(rand()%(maxInd-minInd)+minInd);
    int pixel;

    vector<Point3f> currentMean, newMean, AllPoints;
    for(int y = height-1; y >= 0; y--) 
    {
	    for(int x = 0; x < width; x++) 
	    {
		    pixel = (y * width + x) * 3; 
		    AllPoints.push_back(Point3f(image[pixel],image[pixel+1],image[pixel+2]));
	    }
    }

    AllPoints = swapPoints(AllPoints,0,100);
    int numPoints = AllPoints.size();

    /*
    for(int i = 0;i<k;i++)
	    print(AllPoints[numPoints-i-1]);
	
	*/

    for(int i = 0;i<k;i++)
    {
	    maxInd = numPoints-i-1; minInd = 0;
	    rndIndex = ceil(rand()%(maxInd-minInd)+minInd);
	    AllPoints = swapPoints(AllPoints,rndIndex,numPoints-i-1);
    }
    //*
    cout<<"Initial random mean..\n";
    for(int i = 0;i<k;i++)
    {

	    newMean.push_back(AllPoints[numPoints-i-1]);
	    currentMean.push_back(AllPoints[numPoints-i-1]);
    }

    //*/

    double diff = 10;
    double dist = 0;int kInd = 0;

    int iterNum = 0;
    while(diff > thresh)
    {
	    iterNum++;
	    cout<<"iteration number: "<<iterNum<<" diff value: "<<diff<<endl;

	    diff = 0;
	    // calculate dist b/n all points and means
	    for(int i = 0;i<numPoints;i++)
	    {
		    dist = INT_MAX;
		    kInd = 0;
		    for(int mean = 0; mean<k;mean++)
		    {
			    int d = distanceBetweenPoints(AllPoints[i],currentMean[mean]);
			    if( dist > d)
			    {
				    dist = d;
				    kInd = mean;
			    }
		    }
		    AllClusters[kInd].push_back(AllPoints[i]);
	    }

	    // recalculate means
	    for(int i = 0; i < k; i++)
	    {
		    Point3f sum(0,0,0) ;
		    // sum of all the cluster points
		    for(int m = 0; m<AllClusters[i].size();m++)
			    sum = sum + AllClusters[i][m];
		    sum = sum/AllClusters[i].size();
		    Point3f diffVec = sum - currentMean[i];
		    // diff between old mean and new mean
		    diff+= diffVec.Length();
		    newMean[i] = sum;
		    // reassing current mean as the newly calculated mean
		    currentMean[i] = sum;
	    }

    }

    // reduce the colors in image
    for(int y = height-1; y >= 0; y--) 
    {
	    for(int x = 0; x < width; x++) 
	    {
		    pixel = (y * width + x) * 3; 

		    dist = INT_MAX;
		    kInd = 0;
		    for(int mean = 0; mean<k;mean++)
		    {
			    int d = distanceBetweenPoints(Point3f(image[pixel],image[pixel+1],image[pixel+2]),currentMean[mean]);
			    if( dist > d)
			    {
				    dist = d;
				    kInd = mean;
			    }
		    }
		    image[pixel] = ceil(currentMean[kInd].x);
		    image[pixel+1] = ceil(currentMean[kInd].y);
		    image[pixel+2] = ceil(currentMean[kInd].z);

	    }
    }



	pixmap = image;





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
unsigned char* applyFilterErode(unsigned char *pixmap, double filter[7][7])
{
	int filterSize = 7;
	int m0 = filterSize/2,n0 = m0, indx , indy, pixel;
	unsigned char *pixmapret = new unsigned char[width * height * 3];
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
			//printf("sumr,g,b: %f,%f,%f \n",sumR,sumG,sumB);
			/*
			if(sumR < 0)	sumR =0;
			if(sumG < 0)	sumG =0;
			if(sumB < 0)	sumB =0;
			if(sumR > 255)	sumR = 255;
			if(sumG > 255)	sumG = 255;
			if(sumB > 255)	sumB = 255;
			

			sumR = (sumR+255*3)/6;
                        sumG = (sumR+255*3)/6;
                        sumB = (sumR+255*3)/6;
			if(sumR < 0)	sumR =0;
			if(sumG < 0)	sumG =0;
			if(sumB < 0)	sumB =0;
			if(sumR > 255)	sumR = 255;
			if(sumG > 255)	sumG = 255;
			if(sumB > 255)	sumB = 255;

			*/

			pixmapret[pixel]   = sumR; 
			pixmapret[pixel+1] = sumG;
			pixmapret[pixel+2] = sumB;
		}
	}


	return pixmapret;
}
unsigned char* applyFilterDilate(unsigned char *pixmap, double filter[7][7])
{
	int filterSize = 7;
	int m0 = filterSize/2,n0 = m0, indx , indy, pixel;
	unsigned char *pixmapret = new unsigned char[width * height * 3];
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
			//printf("sumr,g,b: %f,%f,%f \n",sumR,sumG,sumB);
			/*
			if(sumR < 0)	sumR =0;
			if(sumG < 0)	sumG =0;
			if(sumB < 0)	sumB =0;
			if(sumR > 255)	sumR = 255;
			if(sumG > 255)	sumG = 255;
			if(sumB > 255)	sumB = 255;
			

			sumR = (sumR+255*3)/6;
                        sumG = (sumR+255*3)/6;
                        sumB = (sumR+255*3)/6;
			if(sumR < 0)	sumR =0;
			if(sumG < 0)	sumG =0;
			if(sumB < 0)	sumB =0;
			if(sumR > 255)	sumR = 255;
			if(sumG > 255)	sumG = 255;
			if(sumB > 255)	sumB = 255;

			*/

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
	string input = argv[1];
	int N = 1, M = N;
	float I = 0, J = 0;
	cout<< "input arg is : " <<input<<endl;
	//initialize the global variables
	width = 640;
	height = 480;
	pixmap = new unsigned char[width * height * 3]; 
	if(input == "blur")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " blur <input image> <direction theta in radians>" << std::endl;
			return 1;
		}

		double filter[filterSize][filterSize]={0};
		string inputImage = argv[2];
		double theta = atof(argv[3]), I ,J;
		Point2f n0(cos(theta), sin(theta)), n1(-sin(theta),cos(theta)), P0((int)(filterSize/2),(int)(filterSize/2));
		pixmap = readPPM(inputImage);
		// cout<<"theta = "<<theta<<endl<<"atof value: "<<atof(string("3.14").c_str())<<endl;
		// define filter parameters S0 and S1 for gausian filter
		double S0filter = 1, S1filter = 100, sumFilter = 0;
		for(int i = 0; i<filterSize;i++)
		{
			for(int j = 0; j<filterSize;j++)
			{
				double rndx =(float)((float)(rand()%100)/(100.0));
				double rndy=(float)((float)(rand()%100)/(100.0));
				// sampling the filter
				double sum = 0;
				for (int p = 0; p < M; p++)
				{
					for (int q = 0; q < N; q++)
					{
						I=i+(((float)p/M)+((rndx)/M));
						J=j+(((float)q/N)+((rndy)/N));
						sum+=gausFnValue(Point2f(I,J),P0,n0,n1,S0filter,S1filter);
					}
				}
				filter[i][j] = sum/(M*N);
				sumFilter+=filter[i][j];
			}
		}

		cout<<"displaying the filter: \n";
		for(int i = 0; i<filterSize;i++)
		{
			for(int j = 0; j<filterSize;j++)
			{
				filter[i][j]/=sumFilter;
				printf("%f ", filter[i][j]);
			}
			cout<<endl;
		}
		pixmap = applyFilterBlur(pixmap , filter);

		
	}
	else  if(input == "emboss")
	{
		if (argc < 3) {
			std::cout<< "Usage: " << argv[0] << " emboss <input image> <direction theta in radians>" << std::endl;
			return 1;
		}
		string inputImage = argv[2];
		double theta = atof(argv[3]), I ,J, col=0;
		double filter[filterSize][filterSize]={0};
		Point2f n0(cos(theta), sin(theta)), n1(-sin(theta),cos(theta)), P0((int)(filterSize/2),(int)(filterSize/2));
		pixmap = readPPM(inputImage);
		cout<<"theta = "<<theta<<endl<<"atof value: "<<atof(string("3.14").c_str())<<endl;
		// emboss parameters
		double D = 1;
		double sumFilterPlus= 0,numPlus=0 , sumFilterNeg = 0,numNeg=0;
		double embossFilterValue = 0.5;
		// creating the filter programatically at the given angle
		for(int i = 0; i<filterSize;i++)
		{
			for(int j = 0; j<filterSize;j++)
			{
				double rndx =(float)((float)(rand()%100)/(100.0));
				double rndy=(float)((float)(rand()%100)/(100.0));
				// sampling the filter
				double sum = 0;
				for (int p = 0; p < M; p++)
				{
					for (int q = 0; q < N; q++)
					{
						I=i+(((float)p/M)+((rndx)/M));
						J=j+(((float)q/N)+((rndy)/N));
						// print(n0);
						double s = embossFunction(Point2f(I,J),P0,n0); 	
						if(abs(s)<=D)
						{
							sum+=0;
							cout<<":";
						}
						else if(s>D)
							sum+=1;
						else if(s<-1*D)
							sum-=1;
						//printf("s: %f, sum= %f \n",s,sum);
					}
				}
				// filter[i][j] = sum/(M*N);
				if(sum>0)	
				{
					if(sum> ((float)(M*N)/2))
						sum = embossFilterValue;
					else
						sum = 0;
					filter[i][j] = sum;
					sumFilterPlus+=filter[i][j];numPlus++;	
				}
				else if(sum<0.0)
				{
					if(sum< -1*(((float)(M*N)/2)))
						sum = -1*embossFilterValue;
					else
						sum = 0;
					
					filter[i][j] = sum;
					
					sumFilterNeg+=filter[i][j];
					numNeg++;	
				}
				else if (sum ==0)
				{
					filter[i][j] = sum;
				}


			}
		}
		//printf("numPlus: %f, numNeg: %f, sumFilterPlus %f, sumFilterNeg: %f \n",numPlus, numNeg,sumFilterPlus,sumFilterNeg);
		// displayFilter(filter);

		// sum of all values = 0
		double c = sumFilterPlus - abs(sumFilterNeg); 
		if(c > 0)
		{
			c/=numPlus;
			for(int i = 0; i<filterSize;i++)
				for(int j = 0; j<filterSize;j++)
					if(filter[i][j] > 0)
							filter[i][j]-= c;
		}
		else if(c < 0)
		{
			c/=numNeg;
			for(int i = 0; i<filterSize;i++)
				for(int j = 0; j<filterSize;j++)
					if(filter[i][j]<0)
						filter[i][j]+= c;
		}


		displayFilter(filter);
		pixmap = applyFilterEmboss(pixmap , filter);


	}
	else  if(input == "edge")
	{
		if (argc < 3) {
			std::cout<< "Usage: " << argv[0] << " edge <input image>" << std::endl;
			return 1;
		}
		string inputImage = argv[2];
		double theta = atof(argv[3]), I ,J, col=0;
		double filter[3][3]=		{	{-1,-1,-1},
							{-1,+8,-1},
							{-1,-1,-1} };
		Point2f n0(cos(theta), sin(theta)), n1(-sin(theta),cos(theta)), P0((int)(filterSize/2),(int)(filterSize/2));
		pixmap = readPPM(inputImage);
		cout<<"theta = "<<theta<<endl<<"atof value: "<<atof(string("3.14").c_str())<<endl;

		pixmap = applyFilterEdge(pixmap, filter);
	}
	else  if(input == "erode")
	{
		if (argc < 3) {
			std::cout<< "Usage: " << argv[0] << " erode <input image> " << std::endl;
			return 1;
		}
		string inputImage = argv[2];
		double filter[7][7]=		{	{0,0,0,1,0,0,0},
							{0,0,1,1,1,0,0},
							{0,1,1,1,1,1,0},
							{1,1,1,1,1,1,1},
							{0,1,1,1,1,1,0},
							{0,0,1,1,1,0,0},
							{0,0,0,1,0,0,0}, };
		pixmap = readPPM(inputImage);

		pixmap = applyFilterErode(pixmap, filter);
	}
	else  if(input == "dilate")
	{
		if (argc < 3) {
			std::cout<< "Usage: " << argv[0] << " dilate <input image> " << std::endl;
			return 1;
		}
		string inputImage = argv[2];
		double filter[7][7]=		{	{0,0,0,1,0,0,0},
							{0,0,1,1,1,0,0},
							{0,1,1,1,1,1,0},
							{1,1,1,1,1,1,1},
							{0,1,1,1,1,1,0},
							{0,0,1,1,1,0,0},
							{0,0,0,1,0,0,0}, };
		pixmap = readPPM(inputImage);

		pixmap = applyFilterDilate(pixmap, filter);

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
	glutCreateWindow("pr05");
	init();
	glutReshapeFunc(windowResize);
	glutDisplayFunc(windowDisplay);
	glutMouseFunc(processMouse);
	glutMainLoop();

	return 0; //This line never gets reached. 
}

