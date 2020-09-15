#pragma once

#include "Point.h"


class Sphere {

private :
	float f_Radius;
	Point pt_Center;

public:
	Sphere(float rad, Point c) {
		f_Radius = rad;
		pt_Center = c;
	}

	float GetRadius() { return f_Radius; }
	Vector3 GetCenter() { return pt_Center.GetPos(); }



};
