#pragma once

#pragma once
#include "vector3.h"

class Direction {
private:
	Vector3 vec3_Dir;

public:
	Direction()
	{
		vec3_Dir = Vector3(0, 0, 0);
	}

	Direction(Vector3 vec)
	{
		vec3_Dir = vec;
	}

	Direction(float x, float y, float z)
	{
		vec3_Dir = Vector3(x, y, z);
	}

	Vector3 GetPos()
	{
		return vec3_Dir;
	}

};
