#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "lodepng.h"
#include "vector3.h"
#include <cstdint>
#include <optional>
#include "Rayon.h"
#include "Sphere.h"
#include "Color.h"
#include "Lampe.h"
#include "Box3.h"
#include <map>
#include <utility>

#define _USE_MATH_DEFINES
#include <math.h>
#include <random>

using namespace std;

default_random_engine generator;
uniform_real_distribution<float> distribution(-0.00f, 0.00f);

// Comparator function to sort pairs
// according to second value
bool cmp(pair<int, float>& a,
	pair<int, float>& b)
{
	return a.second < b.second;
}

// Function to sort the map according
// to value in a (key-value) pairs
void sort(map<int, float>& M)
{
	// Declare vector of pairs
	vector<pair<int, float> > A;

	// Copy key-value pair from Map
	// to vector of pairs
	for (auto& it : M) {
		A.push_back(it);
	}
	// Sort using comparator function
	sort(A.begin(), A.end(), cmp);
}

/*float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}*/

float RandomFloat(float a, float b) {
	uniform_real_distribution<float> distributionFloat(a, b);

	return distributionFloat(generator);
}



Vector3 RandAlbedo()
{
	return Vector3(RandomFloat(0, 1), RandomFloat(0, 1), RandomFloat(0, 1));
}

void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

optional<float> raySphereIntersect(Rayon r, Sphere s) {
	float a = r.GetDirection().dot(r.GetDirection());
	Vector3 s0_r0 = r.GetOrigin() - s.GetCenter();
	float b = 2.0 * r.GetDirection().dot(s0_r0);
	float c = s0_r0.dot(s0_r0) - (s.GetRadius() * s.GetRadius());
	double delta = pow(b, 2) - 4.0 * a * c;
	if (delta >= 0.0) {
		float resultat = (-b - sqrt(delta)) / (2.0 * a);
		if (resultat >= 0)
		{
			return resultat;
		}
	}
	return nullopt;
}

optional<float> rayAABIntersect(Rayon r, Box3 b)
{
	float t0x = (b.bounds[0].x - r.GetOrigin().x) / r.GetDirection().x;
	float t1x = (b.bounds[1].x - r.GetOrigin().x) / r.GetDirection().x;
	float t0y = (b.bounds[0].y - r.GetOrigin().y) / r.GetDirection().y;
	float t1y = (b.bounds[1].y - r.GetOrigin().y) / r.GetDirection().y;
	float t0z = (b.bounds[0].z - r.GetOrigin().z) / r.GetDirection().z;
	float t1z = (b.bounds[1].z - r.GetOrigin().z) / r.GetDirection().z;

	float tmin = max(max(min(t0x, t1x), min(t0y, t1y)), min(t0z, t1z));
	float tmax = min(min(max(t0x, t1x), max(t0y, t1y)), max(t0z, t1z));

	if (tmax > 0 && tmin < tmax)
	{
		return tmin;
	}

	return nullopt;
}

void color(std::vector<double>& img, Vector3 pixel, float w, float r, float g, float b, float a)
{
	img[4 * w * pixel.y + 4 * pixel.x + 0] += r;

	img[4 * w * pixel.y + 4 * pixel.x + 1] += g;

	img[4 * w * pixel.y + 4 * pixel.x + 2] += b;

	img[4 * w * pixel.y + 4 * pixel.x + 3] = a;
}

void color2(std::vector<double>& img, Vector3 pixel, float w, Color c, float a)
{
	img[4 * w * pixel.y + 4 * pixel.x + 0] += c.r();

	img[4 * w * pixel.y + 4 * pixel.x + 1] += c.g();

	img[4 * w * pixel.y + 4 * pixel.x + 2] += c.b();

	img[4 * w * pixel.y + 4 * pixel.x + 3] = a;
}

Color computeColor(Lampe lampe, Sphere S, Vector3 A, Point X)
{
	Vector3 col;

	Vector3 N = X.GetPos() - S.GetCenter();
	N = N.normalize();

	Vector3 L = X.GetPos() - lampe.GetPos();
	float d = L.length();
	L = L.normalize();

	float NdotXL = abs(N.dot(L));

	col = (lampe.GetQuant() * A * NdotXL) / (pow(d, 2) * M_PI);

	return Color(col.x, col.y, col.z);
}

void searchCloserBox(vector<Box3> boxes, Rayon r, int& iboxes, optional<float>& min_dst)
{
	for (int i = 0; i < boxes.size(); i++)
	{
		auto dst = rayAABIntersect(r, boxes[i]);
		if (dst)
		{
			if (!min_dst || min_dst.value() > dst.value())
			{
				min_dst = dst;
				iboxes = i;
			}
		}
	}
}

void searchCloserObject(vector<Sphere> spheres, Rayon r, int& ispheres, optional<float>& min_dst)
{
	for (int i = 0; i < spheres.size(); i++)
	{
		auto dst = raySphereIntersect(r, spheres[i]);
		if (dst)
		{
			if (!min_dst || min_dst.value() > dst.value())
			{
				min_dst = dst;
				ispheres = i;
			}
		}
	}
}

map<int, float> searchBox(vector<Box3> boxes, Rayon r)
{
	map<int, float> mapBiDst;
	for (int bi = 0; bi < boxes.size(); bi++)
	{
		auto dst = rayAABIntersect(r, boxes[bi]);
		if (dst)
		{
			mapBiDst.insert(pair<int, float>(bi, dst.value()));
		}
	}

	sort(mapBiDst);

	return mapBiDst;
}

optional<pair<Sphere, float>> SearchSphereDich(Rayon r, vector<Box3> boxes)
{
	map<int, float> RIntersectB = searchBox(boxes, r);
	/*for (auto it = RIntersectB.begin(); it != RIntersectB.end(); ++it)
	{
		cout <<"Box"<< it->first << " => " << it->second << " bMin : "<< boxes[it->first].bounds[0] << " bMax : " << boxes[it->first].bounds[1] << '\n';
	}*/
	optional<pair<Sphere, float>> sdst = nullopt;
	auto it = RIntersectB.begin();
	while (it != RIntersectB.end())
	{
		vector<Sphere> spheres = boxes[it->first].lst_spheres;
		optional<float>min_dstSphereinBox = nullopt;
		int sphereindex = 0;
		searchCloserObject(spheres, r, sphereindex, min_dstSphereinBox);
		if (min_dstSphereinBox)
		{
			sdst = pair<Sphere, float>(spheres[sphereindex], min_dstSphereinBox.value());
			return sdst;
		}
		++it;
	}
	return nullopt;
}

void lancerRayonBox(Rayon r, vector<Lampe> lampes, vector<Box3> boxes, vector<double>& image, unsigned width, Vector3 pixel)
{
	for (int l = 0; l < lampes.size(); l++)
	{
		Point X;
		Vector3 vec_dir;
		Direction dir;
		int boxindex = 0;
		optional<float>lampe = nullopt;
		optional<float>min_dst = nullopt;
		searchCloserBox(boxes, r, boxindex, min_dst);

		if (min_dst)
		{
			float epsilone = -0.05;
			X = Point(r.GetOrigin().x + min_dst.value() * r.GetDirection().x, r.GetOrigin().y + min_dst.value() * r.GetDirection().y, r.GetOrigin().z + min_dst.value() * r.GetDirection().z);
			X = Point(X.GetPos() + r.GetDirection() * epsilone);

			vec_dir = (lampes[l].GetPos() - X.GetPos()).normalize();
			dir = Direction(vec_dir.x, vec_dir.y, vec_dir.z);
			Vector3 L = X.GetPos() - lampes[l].GetPos();
			L = Vector3(L.x, L.y, L.z);
			float d = L.length();
			int i;
			searchCloserBox(boxes, Rayon(X, dir), i, lampe);
			if (lampe)
			{
				if (lampe.value() > d)
				{
					color(image, pixel, width, boxes[boxindex].GetAlbedo().x * 255, boxes[boxindex].GetAlbedo().y * 255, boxes[boxindex].GetAlbedo().z * 255, 100);
				}
				else
				{
					color(image, pixel, width, 0, 0, 0, 255);
				}
			}
			else
			{
				color(image, pixel, width, boxes[boxindex].GetAlbedo().x * 255, boxes[boxindex].GetAlbedo().y * 255, boxes[boxindex].GetAlbedo().z * 255, 100);
			}
		}
		else
		{
			color(image, pixel, width, 0, 0, 0, 255);
		}
	}
}

void lancerRayonSphere(Rayon r, vector<Lampe> lampes, vector<Sphere> spheres, vector<double>& image, unsigned width, Vector3 pixel)
{
	r.SetDirection(Direction(Vector3(r.GetDirection().x + distribution(generator), r.GetDirection().y + distribution(generator), r.GetDirection().z + distribution(generator))));
	for (int l = 0; l < lampes.size(); l++)
	{
		Point X;
		Vector3 vec_dir;
		Direction dir;
		int sphereindex = 0;
		optional<float>lampe = nullopt;
		optional<float>min_dst = nullopt;

		searchCloserObject(spheres, r, sphereindex, min_dst);

		if (min_dst)
		{
			float epsilone = -0.05;

			X = Point(r.GetOrigin().x + min_dst.value() * r.GetDirection().x, r.GetOrigin().y + min_dst.value() * r.GetDirection().y, r.GetOrigin().z + min_dst.value() * r.GetDirection().z);
			X = Point(X.GetPos() + r.GetDirection() * epsilone);
			if (!spheres[sphereindex].Mirror())
			{
				vec_dir = (lampes[l].GetPos() - X.GetPos()).normalize();
				dir = Direction(vec_dir.x, vec_dir.y, vec_dir.z);
				Vector3 L = X.GetPos() - lampes[l].GetPos();
				L = Vector3(L.x, L.y, L.z);
				float d = L.length();
				int sphereCloserLampe = 0;
				searchCloserObject(spheres, Rayon(X, dir), sphereCloserLampe, lampe);
				if (lampe)
				{
					if (lampe.value() > d)
					{
						Color c = computeColor(lampes[l], spheres[sphereindex], spheres[sphereindex].GetAlbedo(), X);
						color(image, pixel, width, c.GetColorRGB().x, c.GetColorRGB().y, c.GetColorRGB().z, 255);
					}
					else
					{
						color(image, pixel, width, 0, 0, 0, 255);
					}
				}
				else
				{
					Color c = computeColor(lampes[l], spheres[sphereindex], spheres[sphereindex].GetAlbedo(), X);
					color(image, pixel, width, c.GetColorRGB().x, c.GetColorRGB().y, c.GetColorRGB().z, 255);
				}
			}
			else
			{
				Vector3 I = r.GetDirection();
				Vector3 N = X.GetPos() - spheres[sphereindex].GetCenter();
				N = N.normalize();
				Direction R((-I.dot(N) * N * 2 + I).normalize());
				Rayon r2(X, R);
				lancerRayonSphere(r2, lampes, spheres, image, width, pixel);
			}
		}
		else
		{
			color(image, pixel, width, 0, 0, 0, 255);
		}
	}
}

Color lancerRayonDich(Rayon r, vector<Lampe>lampes, vector<Box3> boxes)
{
	auto SpherePick = SearchSphereDich(r, boxes);
	Color c = Color();
	if (SpherePick) {
		float dst = SpherePick.value().second;
		Point X;
		Vector3 vec_dir;
		Direction dir;
		float epsilone = -0.05;
		X = Point(r.GetOrigin().x + dst * r.GetDirection().x, r.GetOrigin().y + dst * r.GetDirection().y, r.GetOrigin().z + dst * r.GetDirection().z);
		X = Point(X.GetPos() + r.GetDirection() * epsilone);

		for (int li = 0; li < lampes.size(); li++)
		{
			vec_dir = (lampes[li].GetPos() - X.GetPos()).normalize();
			dir = Direction(vec_dir.x, vec_dir.y, vec_dir.z);
			Vector3 L = X.GetPos() - lampes[li].GetPos();
			L = Vector3(L.x, L.y, L.z);
			float d = L.length();
			int sphereindexCloserLampe = 0;
			auto SphereCloserLampe = SearchSphereDich(Rayon(X, dir), boxes);
			if (SphereCloserLampe)
			{
				if (SphereCloserLampe.value().second > d)
				{
					c = c + computeColor(lampes[li], SpherePick.value().first, SpherePick.value().first.GetAlbedo(), X);
				}
			}
			else
			{
				c = c + computeColor(lampes[li], SpherePick.value().first, SpherePick.value().first.GetAlbedo(), X);
			}
		}
	}
	return c;
}
void uniformeImg(vector<double>& image, vector<unsigned char>& imageOut)
{
	double max = 0;
	for (size_t i = 0; i < image.size(); i += 4)
	{
		if (image[i] > max)
		{
			max = image[i];
		}
		if (image[i + 1] > max)
		{
			max = image[i + 1];
		}
		if (image[i + 2] > max)
		{
			max = image[i + 2];
		}
	}
	max /= 4;
	for (size_t i = 0; i < image.size(); i += 4)
	{
		imageOut[i] = clamp((int)(image[i] / max * 255), 0, 255);
		imageOut[i + 1] = clamp((int)(image[i + 1] / max * 255), 0, 255);
		imageOut[i + 2] = clamp((int)(image[i + 2] / max * 255), 0, 255);
		imageOut[i + 3] = 255;
	}
}
void createObject(vector<Sphere>& spheres, vector<Box3>& boxes, int nSphere)
{
	for (int i = 0; i < nSphere; i++)
	{
		Sphere s(Sphere(RandomFloat(10, 100), Point(RandomFloat(0, 1100), RandomFloat(0, 1100), RandomFloat(700, 5000)), Vector3(RandomFloat(0, 1), RandomFloat(0, 1), RandomFloat(0, 1))));
		spheres.push_back(s);
		boxes.push_back(Box3(s));
	}
}

void createSpheres(vector<Sphere>& spheres, int nSphere)
{
	for (int i = 0; i < nSphere; i++)
	{
		Sphere s(Sphere(RandomFloat(10, 20), Point(RandomFloat(-1000, 2000), RandomFloat(-1000, 2000), RandomFloat(700, 5000)), Vector3(RandomFloat(0, 1), RandomFloat(0, 1), RandomFloat(0, 1))));
		spheres.push_back(s);
	}
}

void CreateStructGridRec(Box3 b, vector<Box3>& boxes, int nSphere)
{
	if (b.lst_spheres.size() <= nSphere)
	{
		boxes.push_back(b);
	}
	else
	{
		Vector3 bMin = b.bounds[0];
		Vector3 bMax = b.bounds[1];
		float dX = bMax.x - bMin.x;
		float dY = bMax.y - bMin.y;
		float dZ = bMax.z - bMin.z;

		Vector3 boundsC1;
		Vector3 boundsC2;
		if (dX >= dY && dX >= dZ)
		{
			boundsC1 = Vector3(bMax.x - dX / 2, bMax.y, bMax.z);
			boundsC2 = Vector3(bMin.x + dX / 2, bMin.y, bMin.z);;
		}
		else if
			(dY >= dX && dY >= dZ)
		{
			boundsC1 = Vector3(bMax.x, bMax.y - dY / 2, bMax.z);
			boundsC2 = Vector3(bMin.x, bMin.y + dY / 2, bMin.z);;
		}
		else
		{
			boundsC1 = Vector3(bMax.x, bMax.y, bMax.z - dZ / 2);
			boundsC2 = Vector3(bMin.x, bMin.y, bMin.z + dZ / 2);;
		}
		Box3 Child1(bMin, boundsC1, b.lst_spheres, RandAlbedo());
		Box3 Child2(boundsC2, bMax, b.lst_spheres, RandAlbedo());

		CreateStructGridRec(Child1, boxes,nSphere);
		CreateStructGridRec(Child2, boxes,nSphere);
	}
}

void GetMaxMin(Vector3& Min, Vector3& Max, vector<Sphere> spheres)
{
	Vector3 min = Vector3(spheres[0].GetCenter().x - spheres[0].GetRadius(), spheres[0].GetCenter().y - spheres[0].GetRadius(), spheres[0].GetCenter().z - spheres[0].GetRadius());
	Vector3 max = Vector3(spheres[0].GetCenter().x + spheres[0].GetRadius(), spheres[0].GetCenter().y + spheres[0].GetRadius(), spheres[0].GetCenter().z + spheres[0].GetRadius());
	for (int si = 1; si < spheres.size(); si++)
	{
		min.x > spheres[si].GetCenter().x - spheres[si].GetRadius() ? min.x = spheres[si].GetCenter().x - spheres[si].GetRadius() : min.x;
		min.y > spheres[si].GetCenter().y - spheres[si].GetRadius() ? min.y = spheres[si].GetCenter().y - spheres[si].GetRadius() : min.y;
		min.z > spheres[si].GetCenter().z - spheres[si].GetRadius() ? min.z = spheres[si].GetCenter().z - spheres[si].GetRadius() : min.z;		
		
		max.x < spheres[si].GetCenter().x + spheres[si].GetRadius() ? max.x = spheres[si].GetCenter().x + spheres[si].GetRadius() : max.x;
		max.y < spheres[si].GetCenter().y + spheres[si].GetRadius() ? max.y = spheres[si].GetCenter().y + spheres[si].GetRadius() : max.y;
		max.z < spheres[si].GetCenter().z + spheres[si].GetRadius() ? max.z = spheres[si].GetCenter().z + spheres[si].GetRadius() : max.z;
	}
	float offset = 10.0f;
	/*Min = Vector3(min.x-offset,min.y-offset,min.z-offset);
	Max = Vector3(max.x+offset,max.y+offset,max.z+offset);*/
	Min = min-offset;
	Max = max+offset;
}

int main(int argc, char* argv[])
{
	//NOTE: this sample will overwrite the file or test.png without warning!
	const char* filename = argc > 1 ? argv[1] : "RayTracing.png";
	vector<Sphere> spheres;

	vector<Lampe> lampes;
	lampes.push_back(Lampe(Point(100, 150, 200), Vector3(600000000, 600000000, 600000000)));
	lampes.push_back(Lampe(Point(700, 0, 200), Vector3(600000000, 0, 0)));
	lampes.push_back(Lampe(Point(500, 500, 0), Vector3(6000000000, 6000000000, 6000000000)));

	vector<Box3> boxes;
	createSpheres(spheres, 1000);
	//spheres.push_back(Sphere(20, Point(100, 50, 400)));
	spheres.push_back(Sphere(100, Point(50, 500, 600)));
	Vector3 leftBot = Vector3();
	Vector3 rightTop = Vector3();
	GetMaxMin(leftBot, rightTop, spheres);
	Box3 bEnglobante(leftBot, rightTop,spheres, RandAlbedo());
	CreateStructGridRec(bEnglobante, boxes ,20);

	cout <<"Nombre de boites creees : " << boxes.size() << endl;
	//generate some image
	unsigned width = 1000, height = 1000;
	vector<double> image((width * height * 4), 0.0);
	vector<unsigned char> imageOut((width * height * 4), 0);
	Point Camera(500, 500, -1000);
	int nbRayonAntialiasing = 1;
	int nbRayonSmoothShadow = 1;

	//SearchSphereDich(Rayon(Point(500, 500, 0), Direction(0, 0, 1)),boxes);
	//SearchSphereDich(Rayon(Point(600, 600, 0), Direction(0, 0, 1)),boxes);
	//cout << lancerRayonDich(Rayon(Point(10, 50, 0), Direction(0, 0, 1)), lampes, boxes).GetColorRGB() << endl;
	//cout << lancerRayonDich(Rayon(Point(311, 500, 0), Direction(0, 0, 1)), lampes, boxes).GetColorRGB()<<endl;
	for (unsigned y = 0; y < height; y++)
	{
		for (unsigned x = 0; x < width; x++)
		{
			//cout << "pixel : " << x << " ; " << y << endl;
			for (size_t i = 0; i < nbRayonAntialiasing; i++)
			{
				Direction d = (Vector3(x, y, 0) - Camera.GetPos()).normalize();
				Rayon r(Point((float)x, (float)y, 0), d);
				/*Color c = lancerRayonDich(r, lampes, boxes);
				color2(image, Vector3(x, y, 0), width, c, 255);*/
				//SearchSphereDich(Rayon(Point((float)x, (float)y, 0), d), boxes);
				lancerRayonSphere(Rayon(Point((float)x, (float)y, 0), d), lampes, spheres, image, width, Vector3(x, y, 0));
				//lancerRayonBox(Rayon(Point((float)x, (float)y, 0), d), lampes, boxes, image, width, Vector3(x, y, 0));
			}
		}
	}

	uniformeImg(image, imageOut);
	encodeOneStep(filename, imageOut, width, height);

	return 0;
}