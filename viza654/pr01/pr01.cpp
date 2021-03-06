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

using namespace std;

// =============================================================================
// These variables will store the input ppm image's width, height, and color
// =============================================================================
int width, height;
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
	
	cout<< "input arg is : " <<input<<endl;
	//initialize the global variables
	width = 640;
	height = 480;
	pixmap = new unsigned char[width * height * 3]; 
	if(input == "red")
		fillColor(255,0,0);
	else if(input == "green")
		fillColor(0,255,0);
	else  if(input == "blue")
		fillColor(0,0,255);
	else  if(input == "all")
	{
		fillColor(0,0,255,0,0,width/2,height/2);
		fillColor(255,255,0,(width/2) +1,0,width,height/2);
		fillColor(255,0,0,0,(height/2) +1,width/2,height);
		fillColor(0,255,0,(width/2)+1,(height/2) +1,width,height);

	}
	else  if(input == "circle")
	{
		double radius = 135,centerY= height/2,centerX= width/2;
		int r = 0,g = 0,b=0;
		for(int y = 0; y < height ; y++) {
			for(int x = 0; x< width; x++) {
				int i = (y * width + x) * 3; 
				// cout<<"x: "<<x<<"y: "<<y<<"i: "<<i<<endl;
				if( (x-centerX)*(x- centerX) + (y-centerY)*(y-centerY) - radius*radius <=0 )
				{
					r = 255; b = 0;g=255;
				}
				else
				{
					r = 0; b = 255;g=0;
				}
				pixmap[i++] = r;
				pixmap[i++] = g; 
				pixmap[i] = b;     
			}
		}

	}
	else  if(input == "random")
	{
		double diff = 135,centerY= height/2,centerX= width/2, radius = 30;
		int r = 0,g = 0,b=0;
		for(int y = 0; y < height ; y++) {
			for(int x = 0; x< width; x++) {
				int i = (y * width + x) * 3; 
				// cout<<"x: "<<x<<"y: "<<y<<"i: "<<i<<endl;
				double t = 0;
				double distToCenter = sqrt((x-centerX)*(x- centerX) + (y-centerY)*(y-centerY)); 
				if(distToCenter > 7*radius )
				{
					r = 0; b = 0;g=0;	// black (0,0,0)
				}
				else if(distToCenter < 7*radius && distToCenter > 6*radius )
				{
					t = (distToCenter - 6*radius)/radius;
					r = 0;
					g = 0;
					b = 255*(1-t)+0*t;	// (0,0,1)
				}
				else if(distToCenter < 6*radius && distToCenter > 5*radius )
				{
					t = (distToCenter - 5*radius)/radius; // (0,1,0)
					r = 0; 
					g = 255*(1-t)+0*t;
					b = 0*(1-t)+255*t;
				}
				else if(distToCenter < 5*radius && distToCenter > 4*radius )
				{
					t = (distToCenter - 4*radius)/radius;	// (0,1,1)
					r  = 0; 
					g  = 255;
					b  = 255*(1-t)+0*t;
				}
				else if(distToCenter < 4*radius && distToCenter > 3*radius )
				{
					t = (distToCenter - 3*radius)/radius;	// (1,0,0)
					r = 255*(1-t);
					g = 255*t;
					b = 255*t;
				}
				else if(distToCenter < 3*radius && distToCenter > 2*radius )
				{
					t = (distToCenter - 2*radius)/radius;	// (1,0,1)
					r = 255;
					g = 0;
					b = 255*(1-t);
				}
				else if(distToCenter < 2*radius && distToCenter > 1*radius )
				{
					t = (distToCenter - 1*radius)/radius;	// (1,1,0)
					r = 255;
					g = 255*(1-t);
					b = 255*t;
				}
				else
				{
					r = 255; b = 255;g=255;
				}
				pixmap[i++] = r;
				pixmap[i++] = g; 
				pixmap[i] = b;     
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
	glutCreateWindow("pr01");
	init();
	glutReshapeFunc(windowResize);
	glutDisplayFunc(windowDisplay);
	glutMouseFunc(processMouse);
	glutMainLoop();

	return 0; //This line never gets reached. 
}

