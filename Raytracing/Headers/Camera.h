#pragma once
#include "../geom.h"

class Camera
{
public:
	Camera(Vector3f eye , Quaternionf orientation , float ry, float width, float height);
	~Camera();

	Vector3f eye_;
	Quaternionf orientation_;

	float ry_;
	float rx_;
};

