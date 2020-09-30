#pragma once

#include "vector3.h"
class Color{

private:
	Vector3 vec3_Color;

public:
	Color()
	{
		vec3_Color = Vector3(0, 0, 0);
	}


	Color(float x, float y, float z)
	{
		vec3_Color = Vector3(x, y, z);
	}

	Color(Vector3 vec3)
	{
		vec3_Color = vec3;
	}


	Vector3 GetColorRGB()
	{
		return vec3_Color;
	}
};