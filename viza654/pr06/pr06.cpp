// =============================================================================
// VIZA654/CSCE646 at Texas A&M University
// Project 6
// Created by Ashwath Rajendran 
// This file is supplied with an associated makefile. Put both files in the same
// directory, navigate to that directory from the Linux shell, and type 'make'.
// This will create a program called 'pr06' that you can run by entering
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
	
	int widthFore=0,widthBack=0,heightFore=0,heightBack=0;
	string input = argv[1];
	int N = 1, M = N;
	float I = 0, J = 0;
	cout<< "input arg is : " <<input<<endl;
	//initialize the global variables
	width = 640;
	height = 480;
	pixmap = new unsigned char[width * height * 3]; 
	if(input == "difference")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " subtract <background image> <foreground image2>" << std::endl;
			return 1;
		}

		double filter[filterSize][filterSize]={0};
		string backgroundImage = argv[2],foregroundImage = argv[3];
		unsigned char* pixmapBack;
		unsigned char* pixmapFore;  
		pixmapFore = readPPM(foregroundImage,widthFore,heightFore,maxcolor);
		pixmapBack = readPPM(backgroundImage,widthBack,heightBack,maxcolor);
		
		if(widthFore != widthBack || heightFore != heightBack)
		{
			printf("height or width mismatch between fore and background images\n wF %d,hF %d,wB %d,hB %d\n",widthFore,widthBack,heightFore,heightBack);
			return 1;
		}
		width = widthFore;height = heightFore;
		
		for(int y = height-1; y >= 0; y--) 
		    {
			    for(int x = 0; x < width; x++) 
			    {
				    int pixel = (y * width + x) * 3; 
				    pixmap[pixel] = abs(pixmapBack[pixel] - pixmapFore[pixel++]);	// red
				    //printf("pix back %d, fore %d",int(pixmapBack[pixel]),int( pixmapFore[pixel]));pixel++;
				    pixmap[pixel] = abs(pixmapBack[pixel] - pixmapFore[pixel++]);	// blue
				    pixmap[pixel] = abs(pixmapBack[pixel] - pixmapFore[pixel++]);	// green
			    }
		    }

		
	}
	else if(input == "subtract")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " subtract <background image> <foreground image2>" << std::endl;
			return 1;
		}

		double filter[filterSize][filterSize]={0};
		string backgroundImage = argv[2],foregroundImage = argv[3];
		unsigned char* pixmapBack;
		unsigned char* pixmapFore;  
		pixmapFore = readPPM(foregroundImage,widthFore,heightFore,maxcolor);
		pixmapBack = readPPM(backgroundImage,widthBack,heightBack,maxcolor);
		
		if(widthFore != widthBack || heightFore != heightBack)
		{
			printf("height or width mismatch between fore and background images\n wF %d,hF %d,wB %d,hB %d\n",widthFore,widthBack,heightFore,heightBack);
			return 1;
		}
		width = widthFore;height = heightFore;
		
		for(int y = height-1; y >= 0; y--) 
		    {
			    for(int x = 0; x < width; x++) 
			    {
				    int pixel = (y * width + x) * 3; 
				    pixmap[pixel] = max(pixmapBack[pixel] - pixmapFore[pixel++],0);	// red
				    pixmap[pixel] = max(pixmapBack[pixel] - pixmapFore[pixel++],0);	// blue
				    pixmap[pixel] = max(pixmapBack[pixel] - pixmapFore[pixel++],0);	// green
			    }
		    }

		
	}

	else if(input == "addition")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " subtract <background image> <foreground image2>" << std::endl;
			return 1;
		}

		double filter[filterSize][filterSize]={0};
		string backgroundImage = argv[2],foregroundImage = argv[3];
		unsigned char* pixmapBack;
		unsigned char* pixmapFore;  
		pixmapFore = readPPM(foregroundImage,widthFore,heightFore,maxcolor);
		pixmapBack = readPPM(backgroundImage,widthBack,heightBack,maxcolor);
		
		if(widthFore != widthBack || heightFore != heightBack)
		{
			printf("height or width mismatch between fore and background images\n wF %d,hF %d,wB %d,hB %d\n",widthFore,widthBack,heightFore,heightBack);
			return 1;
		}
		width = widthFore;height = heightFore;
		
		for(int y = height-1; y >= 0; y--) 
		    {
			    for(int x = 0; x < width; x++) 
			    {
				    int pixel = (y * width + x) * 3; 
				    pixmap[pixel] = min(pixmapBack[pixel] + pixmapFore[pixel++],255);	// red
				    pixmap[pixel] = min(pixmapBack[pixel] + pixmapFore[pixel++],255);	// blue
				    pixmap[pixel] = min(pixmapBack[pixel] + pixmapFore[pixel++],255);	// green
			    }
		    }

		
	}

	else  if(input == "min")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " subtract <background image> <foreground image2>" << std::endl;
			return 1;
		}

		double filter[filterSize][filterSize]={0};
		string backgroundImage = argv[2],foregroundImage = argv[3];
		unsigned char* pixmapBack;
		unsigned char* pixmapFore;  
		pixmapFore = readPPM(foregroundImage,widthFore,heightFore,maxcolor);
		pixmapBack = readPPM(backgroundImage,widthBack,heightBack,maxcolor);
		
		if(widthFore != widthBack || heightFore != heightBack)
		{
			printf("height or width mismatch between fore and background images\n wF %d,hF %d,wB %d,hB %d\n",widthFore,widthBack,heightFore,heightBack);
			return 1;
		}
		width = widthFore;height = heightFore;
		
		for(int y = height-1; y >= 0; y--) 
		    {
			    for(int x = 0; x < width; x++) 
			    {
				    int pixel = (y * width + x) * 3; 
				    pixmap[pixel] = min(pixmapBack[pixel],pixmapFore[pixel++]);	// red
				    pixmap[pixel] = min(pixmapBack[pixel],pixmapFore[pixel++]);	// blue
				    pixmap[pixel] = min(pixmapBack[pixel],pixmapFore[pixel++]);		// green
			    }
		    }

	}
	else  if(input == "max")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " subtract <background image> <foreground image2>" << std::endl;
			return 1;
		}

		double filter[filterSize][filterSize]={0};
		string backgroundImage = argv[2],foregroundImage = argv[3];
		unsigned char* pixmapBack;
		unsigned char* pixmapFore;  
		pixmapFore = readPPM(foregroundImage,widthFore,heightFore,maxcolor);
		pixmapBack = readPPM(backgroundImage,widthBack,heightBack,maxcolor);
		
		if(widthFore != widthBack || heightFore != heightBack)
		{
			printf("height or width mismatch between fore and background images\n wF %d,hF %d,wB %d,hB %d\n",widthFore,widthBack,heightFore,heightBack);
			return 1;
		}
		width = widthFore;height = heightFore;
		
		for(int y = height-1; y >= 0; y--) 
		    {
			    for(int x = 0; x < width; x++) 
			    {
				    int pixel = (y * width + x) * 3; 
				    pixmap[pixel] = max(pixmapBack[pixel],pixmapFore[pixel++]);	// red
				    pixmap[pixel] = max(pixmapBack[pixel],pixmapFore[pixel++]);	// blue
				    pixmap[pixel] = max(pixmapBack[pixel],pixmapFore[pixel++]);		// green
			    }
		    }

	}
	else  if(input == "multiply")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " subtract <background image> <foreground image2>" << std::endl;
			return 1;
		}

		double filter[filterSize][filterSize]={0};
		string backgroundImage = argv[2],foregroundImage = argv[3];
		unsigned char* pixmapBack;
		unsigned char* pixmapFore;  
		pixmapFore = readPPM(foregroundImage,widthFore,heightFore,maxcolor);
		pixmapBack = readPPM(backgroundImage,widthBack,heightBack,maxcolor);
		
		if(widthFore != widthBack || heightFore != heightBack)
		{
			printf("height or width mismatch between fore and background images\n wF %d,hF %d,wB %d,hB %d\n",widthFore,widthBack,heightFore,heightBack);
			return 1;
		}
		width = widthFore;height = heightFore;
		
		for(int y = height-1; y >= 0; y--) 
		    {
			    for(int x = 0; x < width; x++) 
			    {
				    int pixel = (y * width + x) * 3; 
				    pixmap[pixel] = pixmapBack[pixel]*pixmapFore[pixel++]/255;	// red
				    pixmap[pixel] = pixmapBack[pixel]*pixmapFore[pixel++]/255;	// blue
				    pixmap[pixel] = pixmapBack[pixel]*pixmapFore[pixel++]/255;	// green
			    }
		    }

	}
	else  if(input == "burn")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " subtract <background image> <foreground image2>" << std::endl;
			return 1;
		}

		double filter[filterSize][filterSize]={0};
		string backgroundImage = argv[2],foregroundImage = argv[3];
		unsigned char* pixmapBack;
		unsigned char* pixmapFore;  
		pixmapFore = readPPM(foregroundImage,widthFore,heightFore,maxcolor);
		pixmapBack = readPPM(backgroundImage,widthBack,heightBack,maxcolor);
		
		if(widthFore != widthBack || heightFore != heightBack)
		{
			printf("height or width mismatch between fore and background images\n wF %d,hF %d,wB %d,hB %d\n",widthFore,widthBack,heightFore,heightBack);
			return 1;
		}
		width = widthFore;height = heightFore;
		
		for(int y = height-1; y >= 0; y--) 
		    {
			    for(int x = 0; x < width; x++) 
			    {
				    int pixel = (y * width + x) * 3; 

				    //printf("pix back %d\n",abs((int)((int) (255 - (256.0*(255.0 - pixmapBack[pixel]))/(pixmapFore[pixel++]+1))%255)));
				    pixmap[pixel] =min(255 - (256+(255 - pixmapBack[pixel]))/(pixmapFore[pixel++]+1),0);
				    pixmap[pixel] =min(255 - (256+(255 - pixmapBack[pixel]))/(pixmapFore[pixel++]+1),0);
				    pixmap[pixel] =min(255 - (256+(255 - pixmapBack[pixel]))/(pixmapFore[pixel++]+1),0);
			    }
		    }

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
	glutCreateWindow("pr06");
	init();
	glutReshapeFunc(windowResize);
	glutDisplayFunc(windowDisplay);
	glutMouseFunc(processMouse);
	glutMainLoop();

	return 0; //This line never gets reached. 
}
