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
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>

using namespace std;
default_random_engine generator;
uniform_real_distribution<float> distribution(-0.00f, 0.00f);

float RandomFloat(float a, float b) {
    float random = ((float)rand()) / (float)RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
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
    
    img[4 *w * pixel.y + 4 * pixel.x + 3] = a;
}

Color computeColor(Lampe lampe, Sphere S, Vector3 A,Point X)
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

void searchCloserObject2(vector<Box3> boxes,Rayon r, int &iboxes, optional<float> &min_dst)
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


void lancerRayon2(Rayon r, vector<Lampe> lampes, vector<Box3> boxes, vector<double>& image, unsigned width, Vector3 pixel)
{
    for (int l = 0; l < lampes.size(); l++)
    {
        Point X;
        Vector3 vec_dir;
        Direction dir;
        int boxindex = 0;
        optional<float>lampe = nullopt;
        optional<float>min_dst = nullopt;
        searchCloserObject2(boxes, r, boxindex, min_dst);

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
            searchCloserObject2(boxes, Rayon(X, dir), i, lampe);
            if (lampe)
            {
                if (lampe.value() > d)
                {

                    color(image, pixel, width, boxes[boxindex].GetAlbedo().x *255, boxes[boxindex].GetAlbedo().y *255, boxes[boxindex].GetAlbedo().z*255, 255);
                }
                else
                {
                    color(image, pixel, width, 0, 0, 0, 255);
                }
            }
            else
            {
                color(image, pixel, width, boxes[boxindex].GetAlbedo().x * 255, boxes[boxindex].GetAlbedo().y * 255, boxes[boxindex].GetAlbedo().z * 255, 255);

            }

        }
        else
        {
            color(image, pixel, width, 0, 0, 0, 255);
        }
    }
}





void lancerRayon(Rayon r, vector<Lampe> lampes, vector<Sphere> spheres, vector<double>& image, unsigned width,Vector3 pixel)
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
                    lancerRayon(r2, lampes, spheres, image, width, pixel);
            }
        }
        else
        {
            color(image, pixel, width, 0, 0, 0, 255);
        }
    }
}

void uniformeImg(vector<double> &image, vector<unsigned char> &imageOut)
{
    double max = 0;
    for (size_t i = 0; i < image.size(); i += 4)
    {
        if (image[i] > max)
        {
            max = image[i];
        }
        if (image[i+1] > max)
        {
            max = image[i+1];
        }
        if (image[i+2] > max)
        {
            max = image[i+2];
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
void createObject(vector<Sphere> &spheres, vector<Box3>& boxes, int nSphere)
{
    for (int i = 0; i < nSphere; i++)
    {
        Sphere s(Sphere(RandomFloat(10, 100), Point(RandomFloat(0, 1100), RandomFloat(0, 1100), RandomFloat(700, 1000)), Vector3(RandomFloat(0, 1), RandomFloat(0, 1), RandomFloat(0, 1))));
        spheres.push_back(s);
        boxes.push_back(Box3(s));
    }
}



int main(int argc, char* argv[])
{
    //NOTE: this sample will overwrite the file or test.png without warning!
    const char* filename = argc > 1 ? argv[1] : "RayTracing.png";
    vector<Sphere> spheres;   
    /*spheres.push_back(Sphere(100.0f, Point(100, 700, 400),Vector3(1,0,0)));
    spheres.push_back(Sphere(200.0f, Point(200, 120, 700),Vector3(1,0,1)));
    spheres.push_back(Sphere(150.0f, Point(600, 702, 400),Vector3(0, 0, 1)));
    spheres.push_back(Sphere(100000.0f, Point(500, 101000, 500),Vector3(0, 1, 1)));
    spheres.push_back(Sphere(100000.0f, Point(101000, 500, 500),Vector3(0, 1, 0)));*/
    //spheres.push_back(Sphere(500.0f, Point(0, 0, 700),Vector3(0,0,0), true));



    vector<Lampe> lampes;
    lampes.push_back(Lampe(Point(350, 400, 0), Vector3(600000000, 600000000, 600000000)));
    lampes.push_back(Lampe(Point(700, 0, 200), Vector3(600000000, 0, 0)));
    lampes.push_back(Lampe(Point(750, 110, 300), Vector3(0, 0, 600000000)));


    vector<Box3> boxes;
    /*boxes.push_back(spheres[0]);
    boxes.push_back(spheres[1]);
    boxes.push_back(spheres[2]);*/
    
    createObject(spheres,boxes, 100);
    //generate some image
    unsigned width = 1000, height = 1000;
    vector<double> image((width * height * 4), 0.0);
    vector<unsigned char> imageOut((width * height * 4), 0);
    Point Camera(500, 500, -1000);
    int nbRayonAntialiasing = 1;
    int nbRayonSmoothShadow = 1;

    for (unsigned y = 0; y < height; y++)
    {
        for (unsigned x = 0; x < width; x++) 
        {
            //cout << "pixel : " << x << " ; " << y << endl;
            for (size_t i = 0; i < nbRayonAntialiasing; i++)
            {
                Direction d = (Vector3(x, y, 0) - Camera.GetPos()).normalize();
                //lancerRayon(Rayon(Point((float)x, (float)y, 0), d), lampes, spheres, image, width, Vector3(x, y, 0));
                lancerRayon2(Rayon(Point((float)x, (float)y, 0), d), lampes, boxes, image, width, Vector3(x, y, 0));
            }         
        }
    }

    uniformeImg(image, imageOut);
    encodeOneStep(filename, imageOut, width, height);
	 
	return 0;
}