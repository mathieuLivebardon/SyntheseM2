#pragma once
#include "vector3.h"
#include "Sphere.h"
class Box3
{
	Vector3 vec3_Albedo;
	Point pt_Center;
public :
	Vector3 bounds[2];

	
	Box3(const Vector3& vmin, const Vector3& vmax)
	{
		bounds[0] = vmin;
		bounds[1] = vmax;
	}

	Box3(Sphere s)
	{
		bounds[0] = Vector3(s.GetCenter().x - s.GetRadius(), s.GetCenter().y - s.GetRadius(), s.GetCenter().z - s.GetRadius());
		bounds[1] = Vector3(s.GetCenter().x + s.GetRadius(), s.GetCenter().y + s.GetRadius(), s.GetCenter().z + s.GetRadius());
		pt_Center = s.GetCenter();
		vec3_Albedo = s.GetAlbedo();
	}


	Vector3 GetAlbedo() { return vec3_Albedo; }


};

