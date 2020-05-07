#pragma once
#include "../geom.h"
class Ray
{
public:
	Ray() {};
	Ray(Vector3f start_point, Vector3f direction);
	Vector3f Eval(float intersection_t);
	~Ray();



	Vector3f Q;//Start point
	Vector3f D;//direction
};

