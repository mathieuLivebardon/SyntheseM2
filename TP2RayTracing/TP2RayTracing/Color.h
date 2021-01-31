#pragma once

#include "vector3.h"
class Color {
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

	float r()
	{
		return vec3_Color.x;
	}
	float g()
	{
		return vec3_Color.y;
	}
	float b()
	{
		return vec3_Color.z;
	}

	Color operator+(const Color col) {
		Color thisCol;
		thisCol.vec3_Color.x = this->vec3_Color.x + col.vec3_Color.x;
		thisCol.vec3_Color.y = this->vec3_Color.y + col.vec3_Color.y;
		thisCol.vec3_Color.z = this->vec3_Color.z + col.vec3_Color.z;
		return thisCol;
	}
};