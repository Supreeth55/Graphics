#pragma once
#include "Shape.h"

class Cylinder : public Shape
{
public:
	Cylinder(Vector3f base,Vector3f axis , float radius, Material* material);
	virtual float get_area() { return 0.0f; };
	virtual bool Intersect(Ray&ray, Intersection& id);

	virtual Bbox bbox() { return box_; }


	~Cylinder();

private:

	float radius_;
	Material* mat_;
	Vector3f base_;
	Vector3f axis_;
	Slab slab_;
	Quaternionf q_;
};

