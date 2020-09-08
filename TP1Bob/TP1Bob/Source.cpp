#include "Vector3.h"
#include "Bob.h"
#include <iostream>
#include <fstream>
#include<list>
using namespace std;




Vector3 calcA(float m, float g,float Cd,Vector3 Z,Vector3 V,Vector3 P);
Vector3 calcV(Vector3 V, Vector3 A);
Vector3 calcP(Vector3 P, Vector3 V);
int main() {


	string const nomFichier("D:/M2/Synthese/TP1Bob/TP1Bob/Vent.csv");
	ofstream monFlux(nomFichier.c_str());
	monFlux << "T posX posY posZ VecX VecY VecZ AX AY AZ"<<endl;

	float m = 80;
	float g = 9.8f;
	float Cd = 0.25f;
	Vector3 Z = Vector3(0, 0, 1);
	Vector3 V0 = Vector3(50, 0, 0);
	Vector3 P0 = Vector3(1, 0, 4000);
	int Dt= 1;
	list<Vector3> V;
	list<Vector3> P;
	list<Vector3> A;

	V.push_back(V0);
	P.push_back(P0);
	A.push_back(calcA(m, g, Cd, Z, V.back(), P.back()));
	Vector3 vent = Vector3(3, 0, 0);
	cout << "T : 0\n";
	V.back().print("V");
	P.back().print("P");
	A.back().print("A");
	cout << "\n";
	monFlux << "0 " << to_string(P.back().GetX())<<" " << to_string(P.back().GetY()) << " " << to_string(P.back().GetZ()) << " " << to_string(V.back().GetX()) << " " << to_string(V.back().GetY()) << " " << to_string(V.back().GetZ()) << " " << to_string(A.back().GetX()) << " " << to_string(A.back().GetY()) << " " << to_string(A.back().GetZ())<<"\n";

	for (size_t i = 1; i < 100; i++)
	{
		cout << "T : "<<i<<"\n";
		Vector3 Vt = V.back();
		Vector3 Pt = P.back();
		Vector3 At = A.back();
		V.push_back(calcV(Vt, At));
		P.push_back(calcP(Pt, Vt));
		A.push_back(calcA(m, g, Cd, Z, V.back(), P.back()));
		V.back().print("V");
		P.back().print("P");
		A.back().print("A");
		cout << "\n";
		monFlux <<to_string(i)<< " " << to_string(P.back().GetX()) << " " << to_string(P.back().GetY()) << " " << to_string(P.back().GetZ()) << " " << to_string(V.back().GetX()) << " " << to_string(V.back().GetY()) << " " << to_string(V.back().GetZ()) << " " << to_string(A.back().GetX()) << " " << to_string(A.back().GetY()) << " " << to_string(A.back().GetZ()) << "\n";
	}



	


	system("pause");

	return 0;
}


Vector3 calcA(float m, float g, float Cd, Vector3 Z, Vector3 V, Vector3 P)
{
	Vector3 A;
	Z = Z * g * -m;
	A = (Z - V * Cd * V.norme())*(1/m);
	return A;
}

Vector3 calcV(Vector3 V, Vector3 A)
{
	return V + A;
}

Vector3 calcP(Vector3 P, Vector3 V)
{
	Vector3 vent(3, 0, 0);
	Vector3 test = P + V + vent;
	return P + V +vent;
}