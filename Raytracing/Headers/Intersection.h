#pragma once
#include "../geom.h"
#include "Ray.h"
#include "Slab.h"


typedef enum {

	E_DIFFUSE,
	E_REFLECTION,
	E_TRANSMISSION,
	E_ALL

}ChoiceType;

class Shape;

class Intersection
{
public:
	Intersection();

	Intersection(Vector3f _p, Vector3f _n, Shape* _s, float _t) : p(_p) ,n(n), s(_s),t(_t){};

	void SetIntersection(Intersection B) {
		t = B.t;
		s = B.s;
		p = B.p;
		n = B.n.normalized();
	}
	//Intersection(Intersection i) : p(i._p), n(n), s(_s), t(_t) {};

	~Intersection();
	Vector3f p;
	Vector3f n;
	Shape* s;
	float t;

};

class Interval {

public:

	float t0_;
	float t1_;
	Vector3f N0_;
	Vector3f N1_;


	Interval() : t0_(0.0f), t1_(INFINITY){}
	Interval(float t0, float t1, Vector3f N);
	void CompareInterval(Interval& other);
	void Intersect(Interval& other);
	void Intersect(Ray& ray, Slab& slab);
	void Intersect(Ray& ray, Slab& slab, bool &between_box);
	void SetInterval(float t0, float t1, Vector3f N);

	inline void Reset() {
		
		t0_ = 0.0f; 
		t1_ = INFINITY;
	}
	inline void Empty() {

		t0_ = 0;
		t1_ = -1;
	}

};

