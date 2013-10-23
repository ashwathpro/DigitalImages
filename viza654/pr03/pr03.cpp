// =============================================================================
// VIZA654/CSCE646 at Texas A&M University
// Project 1
// Created by Ashwath Rajendran 
// This file is supplied with an associated makefile. Put both files in the same
// directory, navigate to that directory from the Linux shell, and type 'make'.
// This will create a program called 'pr03' that you can run by entering
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
int width, height, maxcolor;
unsigned char *pixmap;
// =============================================================================
// fillColor(int r, int g, int b, int stx = 0, int sty = 0, int endx = width , int endy = height) 
//
// This function stores the RGB values of each pixel to "pixmap."
// It starts filling the given r,g,b values from stx,sty to endx,endy
// Then, "glutDisplayFunc" below will use pixmap to display the pixel colors.
// =============================================================================
void fillColor(int r, int g, int b, int stx = 0, int sty = 0, int endx = width , int endy = height)
{
	if( r < 0 || b < 0 || g < 0 || r > 255 || g > 255 || b > 255)
		cout<<"RGB should be between 0 to 255\n";

   for(int y = sty; y < endy ; y++) {
     for(int x = stx; x < endx; x++) {
       int i = (y * width + x) * 3; 
       pixmap[i++] = r;
       pixmap[i++] = g; 
       pixmap[i] = b;     
     }
   }
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
static void windowDisplayGrey(void)
{
  glClear(GL_COLOR_BUFFER_BIT);
  glRasterPos2i(0,0);
  glPixelStorei(GL_UNPACK_ALIGNMENT,1);
  glDrawPixels(width, height, GL_LUMINANCE, GL_UNSIGNED_BYTE, pixmap);
  glFlush();
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


class ConvexQuadrilateral
{
	public: 
	Point2f A,B,C,D;
	/* Points A,B,C,D go in clockwise direction starting with A in bottom left corner and C in top right corner
	 * isLeft(A,B,C) Function IS sensitive to the ordering of A and B
	 * A = line point 1; B = line point 2; C = point to check against.
	 * Vector of AB, with A as tail, A --------> B = (B-A)
	 * Hence, isLeft(B,A,C) will yield true, if C is either above or to the left of the vector AB
	 */
	ConvexQuadrilateral()
	{
		A = Point2f(10,10);
		C = Point2f(100,100);
		B = Point2f(10,100);
		D = Point2f(100,10);
	};
	ConvexQuadrilateral(Point2f a, Point2f b, Point2f c, Point2f d)
	{
		A = a,B=b,C=c,D=d;
	};
	
	bool isInsideQuadrilateral(Point2f p)
	{
		//if( !insideLine(A,B,p) && !insideLine(B,C,p) && insideLine(D,C,p) && insideLine(A,D,p))
		if( !isLeft(A,B,p) && !isLeft(B,C,p) && !isLeft(C,D,p) && !isLeft(D,A,p))
		{

		//	printf("value of f: \n");
			return true;
		}
		return false;
	}
	bool isLeft(Point2f a, Point2f b, Point2f c){
		     return ((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x)) > 0;
	}
	bool insideLine(Point2f A, Point2f B, Point2f p)
	{
		double f = p.x*(B.y-A.y)+p.y*(A.x-B.x)+B.x*A.y-A.x*B.y;
		if( f <= 0 )
		{
		//	printf("value of f: %f \n",f);
			return true;
		}
		else
			return false;
		return false;
	}
};
class Sphere
{
	public:
	double radius;
	Point3f center,color,specularPoint, specularColor;
	Sphere()
	{
		radius= 0.0;
		center=Point3f(0,0,0);
		color=Point3f(0,0,0);
	}
	Sphere(Point3f c, double rad,Point3f col)
	{
		radius= rad;
		center=c;
		color=col;
	}
	Sphere(Point3f c, double rad)
	{
		radius= rad;
		center=c;
		color=Point3f(1,1,1);
	}
	void setSpecularPoint(Point3f p,Point3f q= Point3f(1,1,1))
	{
		specularPoint = p;
		specularColor = q;
	}
	
	Point3f getNormal(Point3f point)
	{
		Point3f n=point-center;
		n.Normalize();
		return n;
	}
	bool isPointInside(Point3f Pe)
	{
		double c=(center-Pe)%(center-Pe)-(radius*radius);
		
		if(c>0)		//Point is inside the sphere 
			return false;
		else 	return true;
	}
	Point3f getColor(){	return color;}
	
};
class Star
{
	/*
	 * A,B,C,D,E are all named in couter clockwise direction in the same order with A in bottom left corner of star 
	 * 
	 */
	public: 
	Point2f A,B,C,D, E;
	Star()
	{
		A = Point2f(10,10);
		C = Point2f(100,100);
		B = Point2f(100,10);
		D = Point2f(50,130);
		E = Point2f(10,100);
	};
	Star(Point2f a, Point2f b, Point2f c, Point2f d,Point2f e)
	{
		A = a,B=b,C=c,D=d, E=e;
	};

	// returns 1 if point c is to the left of vector AB: A ---> B
	bool isLeft(Point2f a, Point2f b, Point2f c){
		     return ((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x)) > 0;
	}
	
	bool isInsideStar(Point2f p)
	{
		// if( !insideLine(A,C,p) && insideLine(A,D,p) && !insideLine(B,D,p)&& !insideLine(E,B,p)&& insideLine(E,C,p))

		if( isLeft(A,C,p) + !isLeft(A,D,p) + isLeft(B,D,p) + !isLeft(B,E,p) + !isLeft(E,C,p) >= 4)
		{

		//	printf("value of f: \n");
			return true;
		}
		return false;

	}
	bool insideLine(Point2f A, Point2f B, Point2f p)
	{
		double f = p.x*(B.y-A.y)+p.y*(A.x-B.x)+B.x*A.y-A.x*B.y;
		if( f <= 0 )
		{
		//	printf("value of f: %f \n",f);
			return true;
		}
		else
			return false;
		return false;
	}
};

class Blob
{
	public:
		vector<Point2f> centers;
		double R , L;
		Blob()
		{
			centers.pb(Point2f(200,200));
			centers.pb(Point2f(250,250));
			centers.pb(Point2f(300,300));
			centers.pb(Point2f(250,300));
			centers.pb(Point2f(300,250));
			centers.pb(Point2f(300,230));
			centers.pb(Point2f(330,230));
			centers.pb(Point2f(210,210));
			R = 50, L = 0.15;
		}
		bool isInsideBlob(Point2f p)
		{
			return F(p)>0;
		}
		double F(Point2f p)
		{
			double dist = 0, r=0;
			for(int i = 0;i<centers.size();i++)
			{
				r = sqrt((centers[i]-p)%(centers[i]-p));
				if(r < R)
					dist += (1-3*(r/R)*(r/R) +3*pow((r/R),4) - pow(r/R,6)); 
				else
					dist+=0;
			}
			return dist - L;
		}
};

bool functionSine(Point2f p,double scaleXPi = 100, double translateY = 100, double scaleY = 50)
{
	return (p.y > translateY + scaleY*sin(p.x/scaleXPi));
}

double distanceBetweenPoints(Point3f A,Point3f B)
{	return sqrt((A-B)%(A-B));}

double distanceBetweenPoints(Point2f A,Point2f B)
{	return sqrt((A-B)%(A-B));}

void readPPM(const string filename)
{ 
	  int ch, bit, comment;
	    FILE *fp;
	    fp=fopen(filename.c_str(),"r");				
	    //fp=fopen("red.ppm","r");					
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
		    pixmap = new unsigned char[width * height * 3];         
		    int y, x, pixel;
		    unsigned char red, green, blue;
		    for(y = height-1; y >= 0; y--) 
		    {
			    for(x = 0; x < width; x++) 
			    {
				    fscanf(fp, "%c%c%c", &red, &green, &blue);
				    pixel = (y * width + x) * 3; 
				    pixmap[pixel] = red;
				    pixel++;
				    pixmap[pixel] = green;
				    pixel++;
				    pixmap[pixel] = blue;
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
		    pixmap = new unsigned char[width * height * 1];         
		    int y, x, pixel;
		    unsigned char red, green, blue;
		    for(y = height-1; y >= 0; y--) 
		    {
			    for(x = 0; x < width; x++) 
			    {
				    fscanf(fp, "%c", &red);
				    pixel = (y * width + x); 
				    pixmap[pixel] = red;
				   
			    }
		    }

	    }
	    fclose(fp);

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
	maxcolor = 255;
	pixmap = new unsigned char[width * height * 3]; 
	if(argc == 2)
	{
		readPPM(input);
	
	}
	else if(argc == 3)
	{
		string output = argv[2];
		readPPM(input);
		writePPM(output);
	
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
	glutCreateWindow("pr03");
	init();
	glutReshapeFunc(windowResize);
	glutDisplayFunc(windowDisplay);
	glutMouseFunc(processMouse);
	glutMainLoop();

	
	return 0; //This line never gets reached. 
}

