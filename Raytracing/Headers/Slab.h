#pragma once
#include "../geom.h"

class Slab
{
public:
	Slab():d0_(0.0f),d1_(0.0f),N_(Vector3f()) { }
	Slab(float d0 , float d1 , Vector3f N);
	void Init(float d0, float d1, Vector3f N);
	~Slab();



	float d0_;
	float d1_;
	Vector3f N_;

};

