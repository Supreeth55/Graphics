#include "../Headers/Slab.h"



Slab::Slab(float d0, float d1, Vector3f N){
	d0_ = d0;
	d1_ = d1;
	N_ = N;
}

void Slab::Init(float d0, float d1, Vector3f N) {


	d0_ = d0;
	d1_ = d1;

	if (d0 < d1) {

		d0_ = d1;
		d1_ = d0;
	}

	N_ = N;
}

Slab::~Slab()
{
}
