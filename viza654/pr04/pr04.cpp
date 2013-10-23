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

// =============================================================================
// main() Program Entry
// =============================================================================
int main(int argc, char *argv[]) 
{

	// Check the number of parameters
	if (argc < 2) {
		// Tell the user how to run the program
		std::cout<< "Usage: " << argv[0] << " option" << std::endl;
		/* "Usage messages" are a conventional way of telling the user
		 * how to run a program if they enter the command incorrectly.
		 */
		return 1;
	}
	string input = argv[1];
	int N = 4, M = N;
	float I = 0, J = 0;
	cout<< "input arg is : " <<input<<endl;
	//initialize the global variables
	width = 640;
	height = 480;
	pixmap = new unsigned char[width * height * 3]; 
	if(input == "reduce")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " reduce <input image> <number of colors to reduce>" << std::endl;
			return 1;
		}
		string inputImage = argv[2];
		int k = atoi(argv[3]);
		cout<<"number of colors = "<<k<<endl;
		KmeansClusterImage(inputImage,k);




	}
	else  if(input == "replace")
	{
		if (argc < 4) {
			std::cout<< "Usage: " << argv[0] << " replace <input image> <R G B of src color> <R G B of dest color> <threshold>" << std::endl;
			return 1;
		}
		string inputImage = argv[2];
		pixmap = readPPM(inputImage);
		int sR = atoi(argv[3]),sG = atoi(argv[4]),sB = atoi(argv[5]);
		int dR = atoi(argv[6]),dG = atoi(argv[7]),dB = atoi(argv[8]);
		int thresh = atoi(argv[9]),pixel,red,blue,green;
		for(int y = height-1; y >= 0; y--) 
		{
			for(int x = 0; x < width; x++) 
			{
				pixel = (y * width + x) * 3; 
				red = pixmap[pixel];green = pixmap[pixel+1];blue = pixmap[pixel+2];
				if( (red >sR-thresh && red < sR+thresh) &&
					(green >sG-thresh && green< sG+thresh) &&
					(blue >sB-thresh && blue < sB+thresh) 	)
				{
					pixmap[pixel] = dR;
					pixmap[pixel+1] = dG;
					pixmap[pixel+2] = dB;
				}

				
			}
		}




	}
	else  if(input == "curve")
	{

		
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
	glutCreateWindow("pr04");
	init();
	glutReshapeFunc(windowResize);
	glutDisplayFunc(windowDisplay);
	glutMouseFunc(processMouse);
	glutMainLoop();

	return 0; //This line never gets reached. 
}

