#ifndef COLORSPACES_H
#define COLORSPACES_H

#include <algorithm>
#include <iostream>
#include <math.h>
#include <vector>

class RGBColor;
class HSLColor;

class RGBColor  {
	public:
		RGBColor() : r(0), g(0), b(0)  {}
		RGBColor(double _r, double _g, double _b) : r(_r), g(_g), b(_b)  {}

		HSLColor to_hsl();

		double r,g,b;

		friend std::ostream& operator <<(std::ostream &out, const RGBColor& _vec)  {
			out << "<" << _vec.r << ", " << _vec.g << ", " << _vec.b << ">";
			return out;
		}

	private:
};

class HSLColor  {
	public:
		HSLColor() : h(0), s(0), l(0)  {}
		HSLColor(double _h, double _s, double _l) : h(_h), s(_s), l(_l)  {}

		RGBColor to_rgb();

		friend std::ostream& operator <<(std::ostream &out, const HSLColor& _vec)  {
			out << "<" << _vec.h << ", " << _vec.s << ", " << _vec.l << ">";
			return out;
		}

		double h,s,l;
};

#endif
