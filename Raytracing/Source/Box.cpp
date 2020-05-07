#include "../Headers/Box.h"

Box::Box(Vector3f corner , Vector3f diagonal, Material* material) :Shape(E_BOX) {

	corner_ = corner;
	diagonal_ = diagonal;
	material_ = material;

	//creating slabs here

	slab1_.Init(-corner.x(), -corner.x() - diagonal.x(), Vector3f(1.0f, 0.0f, 0.0f));
	slab2_.Init(-corner.y(), -corner.y() - diagonal.y(), Vector3f(0.0f, 1.0f, 0.0f));
	slab3_.Init(-corner.z(), -corner.z() - diagonal.z(), Vector3f(0.0f, 0.0f, 1.0f));

	/*Bbox b1 = Bbox(corner, Eigen::Vector3f(corner.x() + diagonal.x(), corner.y() + diagonal.y(), corner.z() + diagonal.z()));
	box_.extend(b1);*/

	Vector3f c1 = corner_;
	Vector3f c2 = corner_ + diagonal_;

	Bbox b1(corner_, Eigen::Vector3f(corner_.x() + diagonal_.x(), corner_.y() + diagonal_.y(), corner_.z() + diagonal_.z()));
	box_.extend(b1);

	//box_ = Bbox(c1, c2);
	//box_ = Bbox(c1, c2);
}


bool Box::Intersect(Ray & ray, Intersection & id)
{

	Interval final_interval;
	//float final_t;

	//slab1
	Interval interval;
	interval.Intersect(ray, slab1_);
	final_interval.CompareInterval(interval);

	//slab2
	Interval interval2;
	interval2.Intersect(ray, slab2_);
	final_interval.CompareInterval(interval2);

	//slab3
	Interval interval1;
	interval1.Intersect(ray, slab3_);
	final_interval.CompareInterval(interval1);

	if (final_interval.t0_ > final_interval.t1_) {

		return false;
	}
	else if (final_interval.t0_ < 0 && final_interval.t1_ < 0) {
		return false;
	}
	/*else if (final_interval.t0_ == 0 && final_interval.t1_ == 0) {

		return false;
	}*/
	else 
	{

		if (final_interval.t0_ < 0 || final_interval.t1_< 0) {
			id.t = std::max(final_interval.t0_, final_interval.t1_);
		}
		else {
			id.t = std::min(final_interval.t0_, final_interval.t1_);
		}
		id.p = ray.Eval(id.t);
		id.s = this;


		if (ray.D.dot(final_interval.N0_) < 0)
			id.n = final_interval.N0_;
		else
			id.n = -final_interval.N0_;


		/*if (id.t == final_interval.t0_) {
			id.n = final_interval.N0_;

		}
		else {

			id.n = -final_interval.N1_;

		}*/
		

	
		// set normal


	/*	if(id.n.dot(ray.D) > 0) {
			id.n = -id.n;
		}
		*/
	//	id.n.normalize();
		return true;

		/*if (final_interval.t0_ * final_interval.t1_ >  0 ) {

			if (final_interval.t0_ < final_interval.t1_) {

				id.t = final_interval.t0_;
				id.p = ray.Eval(id.t);
				id.s = this;
				id.n = final_interval.N0_;
				
			}
			else {

				id.t = final_interval.t1_;
				id.p = ray.Eval(id.t);
				id.s = this;
				id.n = final_interval.N1_;
				
			}
			return true;
		}
		else{

			if (final_interval.t0_ >= 0) {

				id.t = final_interval.t0_;
				id.p = ray.Eval(id.t);
				id.s = this;
				id.n = final_interval.N0_;
				return true;
			}
			else if(final_interval.t1_>=0){

				id.t = final_interval.t1_;
				id.p = ray.Eval(id.t);
				id.s = this;
				id.n = final_interval.N1_;
				return true;
			}
			else {
			//	std::cout << "On of the root is zero returning \n";
				return false;
			}
		}*/
	
	}

	return false;
}

Box::~Box()
{
}
