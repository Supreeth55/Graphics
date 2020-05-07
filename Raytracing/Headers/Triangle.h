#pragma once

#include "Shape.h"
#include "../geom.h"


class Triangle : public Shape
{
private:
	Vector3f v0_;
	Vector3f v1_;
	Vector3f v2_;

	Vector3f n0_;
	Vector3f n1_;
	Vector3f n2_;
	Vector3f normal_;
	Vector3f e1_;
	Vector3f e2_;

	//temp
	float t, v, d, u, o_over_d;
	Vector3f p, s, q;

public:
	Triangle(Material* mat, Vector3f v0, Vector3f v1, Vector3f v2, Vector3f n0, Vector3f n1, Vector3f n2);
	virtual bool Intersect(Ray&ray, Intersection& id);
	virtual Bbox bbox() { return box_; }
	virtual float get_area() { return 0.0f; };
	virtual ~Triangle() {}
};

