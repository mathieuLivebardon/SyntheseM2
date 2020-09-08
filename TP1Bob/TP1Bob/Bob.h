#pragma once
#include "Vector3.h"
class Bob
{
private:
	int m;
	float g;
	float Cd;
	Vector3 V0;
	Vector3 P0;
	int Dt;
public :
	Bob(int _m, float _g, float _Cd, Vector3 _V0, Vector3 _P0, int _Dt)
	{
		m = _m;
		g = _g;
		Cd = _Cd;
		V0 = _V0;
		P0 = _P0;
		Dt = _Dt;
	}

	int GetM() { return m; }
	float GetG() { return g; }
	float GetCd() { return Cd; }
	Vector3 GetV0() { return V0; }
	Vector3 GetP0() { return P0; }
	int GetDt() { return Dt; }

	void SetBob(int _m, float _g, float _Cd, Vector3 _V0, Vector3 _P0, int _Dt) {
		m = _m;
		g = _g;
		Cd = _Cd;
		V0 = _V0;
		P0 = _P0;
		Dt = _Dt;
	}
	


};