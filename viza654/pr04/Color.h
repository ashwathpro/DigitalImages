#ifndef _COLOR_H
#define _COLOR_H


class Color {
	double red, green, blue, alpha;

	public:
		Color();
		Color(double, double, double, double);
		Color(double, double, double);
		//Method fucntions
		double getColorRed() { return red; }
		double getColorGreen() { return green; }
		double getColorBlue() { return blue; }
		double getColoralpha() { return alpha; }
		
		double setColorRed(double redValue) {red = redValue;} 
		double setColorGreen(double greenValue) {green = greenValue;} 
		double setColorBlue(double blueValue) {blue = blueValue;} 
		double setColoralpha(double alphaValue) {alpha = alphaValue;} 
		
		Color operator+(Color &c)
		{
		    return Color(red + c.red, green + c.green, blue + c.blue);
		}
		
};

Color::Color(){
	red = 0.5;
	green = 0.5;
	blue = 0.5;
}


Color::Color(double r, double g, double b, double s){
	red = r;
	green = g;
	blue = b;
	alpha = s;
}


Color::Color(double r, double g, double b){
	red = r;
	green = g;
	blue = b;
}

#endif

