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

	Box3(const Vector3& vmin, const Vector3& vmax, vector<Sphere> spheres, Vector3 albedo, bool createString = false)
	{
		
		bounds[0] = vmin;
		bounds[1] = vmax;
		SetSpheres(spheres);
		vec3_Albedo = albedo;
		if(createString)
		{ 
			cout << "Creation de la boite || bMIN : " << vmin << " bMax : " << vmax << " Nombre de spheres : "<<lst_spheres.size() << endl;
		}
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
			
			//cout << Scenter << " : " << i;
			if (TestSphereAABB(spheres[i]))
			{
				lst_spheres.push_back(spheres[i]);
				//cout << " IN BOX ";
			}
			//cout << endl;
		}
	}

	// Returns true if sphere s intersects AABB b, false otherwise
	bool TestSphereAABB(Sphere s)
	{
		float sqDist = 0.0f;
		Point p = s.GetCenter();

		for (int i = 0; i < 3; i++) {
			// for each axis count any excess distance outside box extents
			float v = p.GetPos()[i];
			if (v < bounds[0][i]) sqDist += (bounds[0][i] - v) * (bounds[0][i] - v);
			if (v > bounds[1][i]) sqDist += (v - bounds[1][i]) * (v - bounds[1][i]);
		}

		return sqDist <= s.GetRadius() * s.GetRadius();
	}


};
