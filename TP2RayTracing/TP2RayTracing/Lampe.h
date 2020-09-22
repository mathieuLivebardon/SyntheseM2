#pragma once
#include "vector3.h"


class Lampe
{
	Point pt_pos;
	Vector3 vec3_quant;

public:
	Lampe(Point Lpos, Vector3 Lq)
	{
		pt_pos = Lpos;
		vec3_quant = Lq;
	}


	Vector3 GetPos() { return pt_pos.GetPos(); }
	Vector3 GetQuant() { return vec3_quant; }


private:

};



