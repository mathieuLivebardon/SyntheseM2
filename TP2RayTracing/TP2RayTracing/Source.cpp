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
        return nullopt;
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
    spheres.push_back(Sphere(110.0f, Point(200, 100, 400)));

    //generate some image
    unsigned width = 512, height = 512;
    std::vector<unsigned char> image;
    image.resize(width * height * 4);
    for (unsigned y = 0; y < height; y++)
        for (unsigned x = 0; x < width; x++) {
            Rayon r(Point((float)x, (float)y, 0), Direction(0, 0, 1));
            optional<float>min_dst = nullopt;
            for (int i = 0; i < spheres.size(); i++)
            {
                auto dst = raySphereIntersect(r, spheres[i]);
                if (dst)
                {
                    if (!min_dst || min_dst.value() > dst.value())
                    {
                        min_dst = dst;
                        float dstMin = spheres[i].GetCenter().z - spheres[i].GetRadius();
                        float dstMax = spheres[i].GetCenter().z;
                        if (dstMin - dstMax != 0)
                        {
                            float coef = (dst.value() - dstMax) / (dstMin - dstMax);
                            color(image, Vector3(x, y, 0), width, coef * 255, 0 * 255, coef * 255, 255);
                        }
                        else
                        {
                            color(image, Vector3(x, y, 0), width, 255, 255, 255, 255);
                        }
                    }
                }
                else if (!min_dst)
                {
                    color(image, Vector3(x, y, 0), width,0, 0, 0, 255);
                }

            }
        }

    encodeOneStep(filename, image, width, height);
	 
	return 0;
}

