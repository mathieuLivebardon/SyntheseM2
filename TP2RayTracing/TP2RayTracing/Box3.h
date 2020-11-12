#pragma once
#include "vector3.h"
#include "Sphere.h"
#include <vector>
using namespace std;
class Box3
{
	Vector3 vec3_Albedo;
	Point pt_Center;


public:
	Vector3 bounds[2];
	vector<Sphere> lst_spheres;

	Box3(const Vector3& vmin, const Vector3& vmax, vector<Sphere> spheres, Vector3 albedo)
	{
		cout << "Creation de la boite || bMIN : " << vmin << " bMax : " << vmax << endl;
		bounds[0] = vmin;
		bounds[1] = vmax;
		SetSpheres(spheres);
		vec3_Albedo = albedo;
	}

	Box3(Sphere s)
	{
		bounds[0] = Vector3(s.GetCenter().x - s.GetRadius(), s.GetCenter().y - s.GetRadius(), s.GetCenter().z - s.GetRadius());
		bounds[1] = Vector3(s.GetCenter().x + s.GetRadius(), s.GetCenter().y + s.GetRadius(), s.GetCenter().z + s.GetRadius());
		pt_Center = s.GetCenter();
		vec3_Albedo = s.GetAlbedo();
	}

	Vector3 GetAlbedo() { return vec3_Albedo; }

	void SetSpheres(vector<Sphere> spheres)
	{
		for (int i = 0; i < spheres.size(); i++)
		{
			Vector3 Scenter = spheres[i].GetCenter();
			cout << Scenter << " : " << i;
			if (Scenter.x >= bounds[0].x && Scenter.y >= bounds[0].y && Scenter.z >= bounds[0].z &&
				Scenter.x <= bounds[1].x && Scenter.y <= bounds[1].y && Scenter.z <= bounds[1].z)
			{
				lst_spheres.push_back(spheres[i]);
				cout << " IN BOX ";
			}
			cout << endl;
		}
	}
};
