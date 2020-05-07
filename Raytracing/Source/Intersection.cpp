#include "../Headers/Intersection.h"




Intersection::Intersection()
{
	s = nullptr;
	t = 1000000;
	p = Vector3f();
	n = Vector3f();
}


Intersection::~Intersection()
{

}


//Interval


Interval::Interval(float t0, float t1, Vector3f N) {

	t0_ = t0;
	t1_ = t1;

	if (t0 > t1) {

		t0_ = t1;
		t1_ = t0;
	}


		N0_ = N;
		N1_ = -N;
}

void Interval::SetInterval(float t0, float t1, Vector3f N) {

	t0_ = t0;
	t1_ = t1;

	if (t0 > t1) {

		t0_ = t1;
		t1_ = t0;
	}


	N0_ = N;
	N1_ = -N;
}


void Interval::Intersect(Interval& other) {

}

void Interval::Intersect(Ray& ray, Slab& slab) {


	float n_dot_d = slab.N_.dot(ray.D);
	float n_dot_q = slab.N_.dot(ray.Q);

	float t0 = 0.0f;
	float t1 = 0.0f;

	if (n_dot_d != 0) {

		t0 = -(slab.d0_ + n_dot_q) / n_dot_d;
		t1 = -(slab.d1_ + n_dot_q) / n_dot_d;
		
		SetInterval(t0, t1, slab.N_);
	}
	else {
		float s0 = n_dot_q + slab.d0_;
		float s1 = n_dot_q + slab.d1_;

		if ((s0<=0 && s1>=0 )|| (s1 <= 0 && s0 >= 0)) {

			t0_ = 0.0f;
			t1_ = 1000000;

		}
		else {
			t0_ = 1.0f;
			t1_ = 0.0f;
		}
	}
}


void Interval::Intersect(Ray& ray, Slab& slab, bool &between_box) {

	float n_dot_d = slab.N_.dot(ray.D);
	float n_dot_q = slab.N_.dot(ray.Q);

	float t0 = 0.0f;
	float t1 = 0.0f;

	if (n_dot_d != 0) {

		t0 = -(slab.d0_ + n_dot_q) / n_dot_d;
		t1 = -(slab.d1_ + n_dot_q) / n_dot_d;

		SetInterval(t0, t1, slab.N_);
	}
	else {
		float s0 = n_dot_q + slab.d0_;
		float s1 = n_dot_q + slab.d1_;

		if ((s0 <= 0 && s1 >= 0) || (s1 <= 0 && s0 >= 0)) {

			t0_ = 0.0f;
			t1_ = 1000000;
			between_box = true;
		}
		else {
			t0_ = 1.0f;
			t1_ = 0.0f;
		}
	}
}

void Interval::CompareInterval(Interval& other) {


	//store max of t0
	if (other.t0_ > t0_) {
		t0_ = other.t0_;
		N0_ = other.N0_;
	}

	//store min of t1
	if (other.t1_ < t1_) {
		t1_ = other.t1_;
		N1_ = other.N1_;
	}
}