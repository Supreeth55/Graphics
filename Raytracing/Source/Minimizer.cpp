#include "../Headers/Minimizer.h"
#include "../Headers/Shape.h"
#include "../Headers/Slab.h"
#include "../Headers/Tracer.h"

#define DRAW_BOUNDINGBOX1

// Constructor
Minimizer::Minimizer(Ray& r) : ray(&r) {

	t_ = INFINITY;
	shape_ = nullptr;
	hit_ = false;
}
// Called by BVMinimize to intersect the ray with a Shape. This
// should return the intersection t, but should also track
// the minimum t and it's corresponding intersection info.
// Return INF to indicate no intersection.
float Minimizer::minimumOnObject(Shape* obj) {

#ifdef DRAW_BOUNDINGBOX
	//with bounding box
	
	float tt = minimumOnObjectVolume(obj);

	if (tt < t_) {

		shape_ = obj;
		hit_ = true;
		t_ = tt;
	}
		//return INFINITY;
	return tt;
	

#endif
	Intersection id;

	if (obj->Intersect(*ray, id)) {
		
		if (id.t < t_ && id.t > 0.00001f) {
			id_ = id;
			t_ = id.t;
			id_.s = obj;
			shape_ = obj;
			hit_ = true;
		}
	}

	return  id.t;

}
// Called by BVMinimize to intersect the ray with a Box3d and
// returns the t value. This should be similar to the already
// written box (3 slab) intersection. (The difference begin a
// distance of zero should be returned if the ray starts within the bbox.)
// Return INF to indicate no intersection.
float Minimizer::minimumOnVolume(const Bbox& box)
{

	float max = 0.0f;

	Vector3f L = box.min();  // Box corner   
	Vector3f U = box.diagonal();  // Box corner

	Slab slab1_;
	slab1_.N_ = Vector3f(1, 0, 0);
	slab1_.d0_ = -L.x();
	slab1_.d1_ = -L.x() - U.x();

	if (slab1_.d0_ < slab1_.d1_) {

		float t = slab1_.d0_;
		slab1_.d0_ = slab1_.d1_;
		slab1_.d1_ = t;
	}


	Slab slab2_;
	slab2_.N_ = Vector3f(0, 1, 0);
	slab2_.d0_ = -L.y();
	slab2_.d1_ = -L.y() - U.y();

	if (slab2_.d0_ < slab2_.d1_) {

		float t = slab2_.d0_;
		slab2_.d0_ = slab2_.d1_;
		slab2_.d1_ = t;
	}

	Slab slab3_;
	slab3_.N_ = Vector3f(0, 0, 1);
	slab3_.d0_ = -L.z();
	slab3_.d1_ = -L.z() - U.z();

	if (slab3_.d0_ < slab3_.d1_) {

		float t = slab3_.d0_;
		slab3_.d0_ = slab3_.d1_;
		slab3_.d1_ = t;
	}

	float final_t0 = 0;
	float final_t1 = INFINITY;

	Interval intr1;
	bool between_box = false;
	//Slab1
	intr1.Reset();
	intr1.Intersect(*ray,slab1_);
	if (between_box) {

		return max;
	}
	/**/
	final_t0 = std::max(final_t0, intr1.t0_);
	final_t1 = std::min(final_t1, intr1.t1_);


	//SLab2
	intr1.Reset();
	between_box = false;
	intr1.Intersect(*ray, slab2_);

	if (between_box) {

		return max;
	}


	final_t0 = std::max(final_t0, intr1.t0_);
	final_t1 = std::min(final_t1, intr1.t1_);

	//Slab3
	intr1.Reset();
	between_box = false;
	intr1.Intersect(*ray, slab3_);

	if (between_box) {

		return max;
	}

	final_t0 = std::max(final_t0, intr1.t0_);
	final_t1 = std::min(final_t1, intr1.t1_);

	if (final_t0 < 0 && final_t1 < 0) {

		return INFINITY;
	}

	else if (final_t0 > final_t1) {

		return INFINITY;
	}
	else {

		float t = 0.0f;

		if (final_t0 < 0 || final_t1 < 0) {
			t = std::max(final_t0, final_t1);
		}
		else {
			t = std::min(final_t0, final_t1);
		}
		return t;

	}

return INFINITY;
}


Minimizer& Minimizer::operator=(const Minimizer& m) {

	ray  = m.ray;
	t_ = m.t_;
	shape_ = m.shape_;
	id_ = m.id_;
	hit_ = m.hit_;

	return *this;
}

Minimizer::Minimizer(const Minimizer& m) {

	ray = m.ray;
	t_ = m.t_;
	shape_ = m.shape_;
	id_ = m.id_;
	hit_ = m.hit_;
}
//

float Minimizer::minimumOnObjectVolume(Shape* shp) {
	float max = 0.0f;
	Vector3f L = shp->bbox().min();  // Box corner   
	Vector3f U = shp->bbox().diagonal();  // Box corner
										  //	U +=L;
	Slab slab1_;
	slab1_.N_ = Vector3f(1, 0, 0);
	slab1_.d0_ = -L.x();
	slab1_.d1_ = -L.x() - U.x();

	if (slab1_.d0_ < slab1_.d1_) {

		float t = slab1_.d0_;
		slab1_.d0_ = slab1_.d1_;
		slab1_.d1_ = t;
	}

	Slab slab2_;
	slab2_.N_ = Vector3f(0, 1, 0);
	slab2_.d0_ = -L.y();
	slab2_.d1_ = -L.y() - U.y();

	if (slab2_.d0_ < slab2_.d1_) {

		float t = slab2_.d0_;
		slab2_.d0_ = slab2_.d1_;
		slab2_.d1_ = t;
	}


	Slab slab3_;
	slab3_.N_ = Vector3f(0, 0, 1);
	slab3_.d0_ = -L.z();
	slab3_.d1_ = -L.z() - U.z();

	if (slab3_.d0_ < slab3_.d1_) {

		float t = slab3_.d0_;
		slab3_.d0_ = slab3_.d1_;
		slab3_.d1_ = t;
	}

	float final_t0 = 0;
	float final_t1 = INFINITY;

	Interval intr1;
	bool between_box = false;
	//Slab1
	intr1.Reset();
	intr1.Intersect(*ray,slab1_, between_box);
	if (between_box) {

		return max;
	}
	/**/
	final_t0 = std::max(final_t0, intr1.t0_);
	final_t1 = std::min(final_t1, intr1.t1_);


	//SLab2
	intr1.Reset();
	between_box = false;
	intr1.Intersect(*ray, slab2_, between_box);

	if (between_box) {

		return max;
	}


	final_t0 = std::max(final_t0, intr1.t0_);
	final_t1 = std::min(final_t1, intr1.t1_);

	//Slab3
	intr1.Reset();
	between_box = false;
	intr1.Intersect(*ray,slab3_,between_box);

	if (between_box) {

		return max;
	}

	final_t0 = std::max(final_t0, intr1.t0_);
	final_t1 = std::min(final_t1, intr1.t1_);

	if (final_t0 < 0 && final_t1 < 0) {

		return INFINITY;
	}

	else if (final_t0 > final_t1) {

		return INFINITY;
	}
	else {

		float t = 0.0f;

		if (final_t0 < 0 || final_t1 < 0) {
			t = std::max(final_t0, final_t1);
		}
		else {
			t = std::min(final_t0, final_t1);
		}
		return t;

	}

	return INFINITY;
}