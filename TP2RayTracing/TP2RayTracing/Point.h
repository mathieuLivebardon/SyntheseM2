#pragma once
#include "vector3.h"

class Point {
private :
	Vector3 vec3_Position;

public :
	Point()
	{
		vec3_Position = Vector3(0, 0, 0);
	}


	Point(float x, float y, float z)
	{
		vec3_Position = Vector3(x, y, z);
	}

	Vector3 GetPos()
	{
		return vec3_Position;
	}

};
