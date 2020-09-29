#pragma once

#include "Point.h"


class Sphere {

private :
	float f_Radius;
	Point pt_Center;
	bool b_Mirror;
	Vector3 vec3_Albedo;


public:
	Sphere(float rad, Point c, Vector3 albedo = Vector3(1, 1, 1), bool m = false) {
		b_Mirror = m;
		f_Radius = rad;
		pt_Center = c;
		vec3_Albedo = albedo;
	}

	float GetRadius() { return f_Radius; }
	Vector3 GetCenter() { return pt_Center.GetPos(); }
	bool Mirror() { return b_Mirror; }
	Vector3 GetAlbedo() { return vec3_Albedo; }



};
