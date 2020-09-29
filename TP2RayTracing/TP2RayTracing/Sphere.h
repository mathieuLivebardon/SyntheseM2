#pragma once

#include "Point.h"


class Sphere {

private :
	float f_Radius;
	Point pt_Center;
	bool b_Mirror;

public:
	Sphere(float rad, Point c, bool m = false) {
		b_Mirror = m;
		f_Radius = rad;
		pt_Center = c;
	}

	float GetRadius() { return f_Radius; }
	Vector3 GetCenter() { return pt_Center.GetPos(); }
	bool Mirror() { return b_Mirror; }



};
