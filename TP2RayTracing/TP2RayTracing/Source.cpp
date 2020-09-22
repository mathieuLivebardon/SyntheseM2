#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "lodepng.h"
#include "vector3.h"
#include <optional>
#include "Rayon.h"
#include "Sphere.h"
#include "Color.h"
#include "Lampe.h"
#define _USE_MATH_DEFINES
#include <math.h>


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
    double delta = pow(b,2) - 4.0 * a * c;
    if (delta >=  0.0) {
        float resultat = (-b - sqrt(delta)) / (2.0 * a);
        if (resultat >= 0)
        {
            return resultat;
        }
    }
    return nullopt;
    
}

void color(std::vector<unsigned char>& img, Vector3 pixel, float w, float r, float g, float b, float a)
{
    if (img[4 * w * pixel.y + 4 * pixel.x + 0] + r  > 255) { img[4 * w * pixel.y + 4 * pixel.x + 0] = 255; }
    else{ img[4 * w * pixel.y + 4 * pixel.x + 0] += r; }
    
    if (img[4 * w * pixel.y + 4 * pixel.x + 1] + g > 255) { img[4 * w * pixel.y + 4 * pixel.x + 1] = 255; }
    else{   img[4 * w * pixel.y + 4 * pixel.x + 1] += g; }

    if (img[4 * w * pixel.y + 4 * pixel.x + 2] +b > 255) { img[4 * w * pixel.y + 4 * pixel.x + 2] = 255; }
    else { img[4 * w * pixel.y + 4 * pixel.x + 2] += b; }

    
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



int main(int argc, char* argv[])
{
    //NOTE: this sample will overwrite the file or test.png without warning!
    const char* filename = argc > 1 ? argv[1] : "RayTracing.png";
    vector<Sphere> spheres;   
    spheres.push_back(Sphere(54.0f, Point(100, 300, 200)));
    spheres.push_back(Sphere(50.0f, Point(100, 120, 200)));
    spheres.push_back(Sphere(55.0f, Point(400, 302, 200)));

    vector<Lampe> lampes;
    lampes.push_back(Lampe(Point(252, 230, 200), Vector3(20000000, 20000000, 20000000)));
    lampes.push_back(Lampe(Point(100, 500, 200), Vector3(20000000, 0, 0)));
    /*Rayon r(Point((float)112, (float)138, 0), Direction(0, 0, 1));
    auto dst = raySphereIntersect(r, spheres[1]);

    auto lampe = raySphereIntersect(Rayon(Point(112,138,dst.value()), Direction(0, 1, 0)), spheres[2]);
    cout << dst.value()<<endl;
    cout << lampe.value()<<endl;*/
    //generate some image
    unsigned width = 512, height = 512;
    std::vector<unsigned char> image;
    image.resize(width * height * 4);
    for (unsigned y = 0; y < height; y++)
        for (unsigned x = 0; x < width; x++) {
            Rayon r(Point((float)x, (float)y, 0), Direction(0, 0, 1));                       
            for (int l = 0; l < lampes.size(); l++)
            {
                Point X;
                Vector3 vec_dir;
                Direction dir;
                int sphereindex = 0;
                optional<float>lampe = nullopt;
                optional<float>min_dst = nullopt;
                for (int i = 0; i < spheres.size(); i++)
                {
                    auto dst = raySphereIntersect(r, spheres[i]);
                    if (dst)
                    {
                        if (!min_dst || min_dst.value() > dst.value())
                        {
                            min_dst = dst;                            
                            sphereindex = i;  
                            X = Point((float)x, (float)y, min_dst.value() - 0.02f);
                            vec_dir = (lampes[l].GetPos() - X.GetPos()).normalize();
                            dir = Direction(vec_dir.x, vec_dir.y, vec_dir.z);
                            Vector3 L = X.GetPos() - lampes[l].GetPos();
                            float d = L.length();
                            for (int j = 0; j < spheres.size(); j++)
                            {
                                auto interlampe = raySphereIntersect(Rayon(X, dir), spheres[j]);
                                if (interlampe)
                                {
                                    if(interlampe.value() < d)
                                    {
                                        lampe = interlampe;
                                    }
                                }
                            }
                        }         
                    }
                }         
                if (min_dst)
                {    
                    
                    if (lampe)
                    {
                        color(image, Vector3(x, y, 0), width, 00,0, 0, 255);
                    }
                    else
                    {
                        Color c = computeColor(lampes[l], spheres[sphereindex], Vector3(1, 1, 1), X);
                        //color(image, Vector3(x, y, 0), width,255,255,255, 255);
                        color(image, Vector3(x, y, 0), width, c.GetColorRGB().x, c.GetColorRGB().y, c.GetColorRGB().z, 255);
                    }                   
                }
                else
                {
                    color(image, Vector3(x, y, 0), width, 255, 0, 255, 255);
                }
            }
        }
                
    
    encodeOneStep(filename, image, width, height);
	 
	return 0;
}

