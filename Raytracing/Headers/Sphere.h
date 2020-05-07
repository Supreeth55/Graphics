#pragma once
#include "Shape.h"

class Sphere: public Shape
{

public:
	Vector3f center;
	Vector3f center2;
	float radius_;
	Sphere(Vector3f center, Material* material, float radius);
	Sphere(Vector3f center, Material* material, float radius, Vector3f center2);
	virtual bool Intersect(Ray&ray, Intersection& id);
	virtual float get_area();
	virtual Bbox bbox() { return box_; }
	~Sphere();
};