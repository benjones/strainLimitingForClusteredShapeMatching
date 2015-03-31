#include "color_spaces.h"

HSLColor RGBColor::to_hsl()  {
	double c_max = std::max(r, std::max(g,b));
	double c_min = std::min(r, std::min(g,b));
	double delta = c_max-c_min;

	double l = 0.5*(c_max+c_min);
	double s = l < 0.5 ? delta/(c_max+c_min) : delta/(2.0-c_max-c_min);

	double pi_over_3 = acos(-1)/3.0;
	double h = 0;
	if(r > g && r > b)
		h = pi_over_3*((g-b)/delta);
	else if(g > b)
		h = pi_over_3*((b-r)/delta+2);
	else
		h = pi_over_3*((r-g)/delta+4);
	h = h < 0 ? h+2*acos(-1) : h;
	return HSLColor(h,s,l);
}

RGBColor HSLColor::to_rgb()  {
	double t_1 = l < 0.5 ? l*(1.0+s) : l+s-l*s;
	double t_2 = 2.0*l - t_1;
	if(t_2 < 0 || t_2 > 1)
		std::cout << "bad t_2!" << std::endl;

	double norm_h = h/(2.0*acos(-1));
	double t_r = norm_h+(1.0/3.0);
	t_r = t_r > 1 ? t_r-1 : t_r;
	double t_g = norm_h;
	double t_b = norm_h-(1.0/3.0);
	t_b = t_b < 0 ? t_b+1 : t_b;

	RGBColor rgb;
	if(t_r < 1.0/6.0)
		rgb.r = t_2 + (t_1-t_2)*6.0*t_r;
	else if(t_r < 1.0/2.0)
		rgb.r = t_1;
	else if(t_r < 2.0/3.0)
		rgb.r = t_2 + (t_1-t_2)*6.0*(2.0/3.0-t_r);
	else
		rgb.r = t_2;

	if(t_g < 1.0/6.0)
		rgb.g = t_2 + (t_1-t_2)*6.0*t_g;
	else if(t_g < 1.0/2.0)
		rgb.g = t_1;
	else if(t_g < 2.0/3.0)
		rgb.g = t_2 + (t_1-t_2)*6.0*(2.0/3.0-t_g);
	else
		rgb.g = t_2;

	if(t_b < 1.0/6.0)
		rgb.b = t_2 + (t_1-t_2)*6.0*t_b;
	else if(t_b < 1.0/2.0)
		rgb.b = t_1;
	else if(t_b < 2.0/3.0)
		rgb.b = t_2 + (t_1-t_2)*6.0*(2.0/3.0-t_b);
	else
		rgb.b = t_2;

	return rgb;
}


