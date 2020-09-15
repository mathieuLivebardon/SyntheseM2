#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "lodepng.h"
#include "vector3.h"
#include <optional>
#include "Rayon.h"
#include "Sphere.h"


using namespace std;

void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
    //Encode the image
    unsigned error = lodepng::encode(filename, image, width, height);

    //if there's an error, display it
    if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

optional<float> raySphereIntersect(Rayon r, Sphere s) {
    // - r0: ray origin
    // - rd: normalized ray direction
    // - s0: sphere center
    // - sr: sphere radius
    // - Returns distance from r0 to first intersecion with sphere,
    //   or -1.0 if no intersection.
    float a = r.GetDirection().dot(r.GetDirection());
    Vector3 s0_r0 = r.GetOrigin() - s.GetCenter();
    float b = 2.0 * r.GetDirection().dot(s0_r0);
    float c = s0_r0.dot(s0_r0) - (s.GetRadius() * s.GetRadius());
    float delta = b * b - 4.0 * a * c;
    if ( delta < 0.0) {
        return std::nullopt;
    }
    
    float resultat = (-b - sqrt(delta)) / (2.0 * a);
    if (resultat >= 0)
    {
        return resultat;
    }
    else
    {
        return std::nullopt;
    }
}

void color(std::vector<unsigned char>& img, Vector3 pixel, float w, float r, float g, float b, float a)
{
    img[4 *w * pixel.y + 4 * pixel.x + 0] = r;
    img[4 *w * pixel.y + 4 * pixel.x + 1] = g;
    img[4 *w * pixel.y + 4 * pixel.x + 2] = b;
    img[4 *w * pixel.y + 4 * pixel.x + 3] = a;
}



int main(int argc, char* argv[])
{
    //NOTE: this sample will overwrite the file or test.png without warning!
    const char* filename = argc > 1 ? argv[1] : "test.png";
    vector<Sphere> spheres;   
    spheres.push_back(Sphere(100.0f, Point(200, 300, 500)));
    spheres.push_back(Sphere(10.0f, Point(100, 100, 500)));

    //generate some image
    unsigned width = 512, height = 512;
    std::vector<unsigned char> image;
    image.resize(width * height * 4);
    for (unsigned y = 0; y < height; y++)
        for (unsigned x = 0; x < width; x++) {
            
            Rayon r(Point((float)x, (float)y, 0),Direction(0,0,1));

            auto dst = raySphereIntersect(r,spheres[0]);

            float dstMin = spheres[0].GetCenter().z - spheres[0].GetRadius();
            float dstMax = spheres[0].GetCenter().z;
            
            if (dst)
            {
                if(dstMin-dstMax !=0)
                {
                    float coef = (dst.value() - dstMax) / (dstMin - dstMax);
                    color(image, Vector3(x, y, 0), width,  coef* 255, 0*255, coef*255,255);
                }
                else
                {
                color(image, Vector3(x, y, 0), width, 255,255,255, 255);
                }
            }
            else
            {
                color(image, Vector3(x, y, 0), width, 0,0,0,255);
            }
            
        }

    encodeOneStep(filename, image, width, height);
	 
	return 0;
}

