#pragma once
#include <string>
#include <iostream>
#include <math.h>
using namespace std;
class Vector3
{
private:
    float x, y, z;

public:
    Vector3()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector3(float _x, float _y, float _z)
        : x(_x), y(_y), z(_z)
    {

    }

    float GetX() {
        return x;
    }
    float GetY() {
        return y;
    }
    float GetZ() {
        return z;
    }

    void SetXYZ(float _x, float _y, float _z) { x = _x; y = _y; z = _z;}

    Vector3 operator*(const Vector3 Vec2) {
        Vector3 Vec;
        Vec.x = this->x * Vec2.x;
        Vec.y = this->y * Vec2.y;
        Vec.z = this->z * Vec2.z;
        return Vec;
    }
    Vector3 operator*(const int cst) {
        Vector3 Vec;
        Vec.x = this->x * cst;
        Vec.y = this->y * cst;
        Vec.z = this->z * cst;
        return Vec;
    }
    Vector3 operator*(const float cst) {
        Vector3 Vec;
        Vec.x = this->x * cst;
        Vec.y = this->y * cst;
        Vec.z = this->z * cst;
        return Vec;
    }

    Vector3 operator-(const Vector3 Vec2) {
        Vector3 Vec;
        Vec.x = this->x - Vec2.x;
        Vec.y = this->y - Vec2.y;
        Vec.z = this->z - Vec2.z;
        return Vec;
    }

    Vector3 operator+(const Vector3 Vec2) {
        Vector3 Vec;
        Vec.x = this->x + Vec2.x;
        Vec.y = this->y + Vec2.y;
        Vec.z = this->z + Vec2.z;
        return Vec;
    }

    float norme()
    {
        return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }
   

    void print(string vecName)
    {
        cout << vecName << " :\n" << this->x << " ; "<<this->y << " ; "<<this->z<<"\n";

    }


};
