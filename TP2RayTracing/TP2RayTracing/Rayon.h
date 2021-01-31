#pragma once
#include "Point.h"
#include "Direction.h"

class Rayon
{
private:
	Point pt_Start;
	Direction dir_Direction;

public:
	Rayon(Point p, Direction d)
	{
		pt_Start = p;
		dir_Direction = d;
	};

	Vector3 GetOrigin() {
		return pt_Start.GetPos();
	}
	Vector3 GetDirection() {
		return dir_Direction.GetPos();
	}

	void SetDirection(Direction d) { dir_Direction = d; }
};
