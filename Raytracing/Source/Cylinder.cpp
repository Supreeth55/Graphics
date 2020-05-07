#include "../Headers/Cylinder.h"

Cylinder::Cylinder(Vector3f base, Vector3f axis,float radius, Material* material) : Shape(E_CYLINDER){

	radius_ = radius;
	base_ = base;
	axis_ = axis;
	material_ = material;

	slab_.d0_ = 0.0f;
	slab_.d1_ = -1 * (axis.norm());//-(sqrtf(axis_.x()*axis_.x() + axis_.y()*axis_.y() + axis_.z()*axis_.z()));
	slab_.N_ = Vector3f(0.0f, 0.0f, 1.0f);

	q_ = Quaternionf::FromTwoVectors(axis_, Vector3f::UnitZ());


	Vector3f upper_base = base + axis_;
	Vector3f r(radius, radius, radius);
	auto box1 = Bbox(r + base, -r + base);

	auto box2 = Bbox(r + upper_base, upper_base - r);


	box_ = Bbox(-r + base, r + upper_base);//
	box_.extend(box1);
	box_.extend(box2);



}

bool Cylinder::Intersect(Ray&ray, Intersection& id) {

	Ray new_ray;// (q_._transformVector(ray.Q - base_), q_._transformVector(ray.D));
	new_ray.Q = q_._transformVector(ray.Q - base_);
	new_ray.D = q_._transformVector(ray.D);

	Interval final_interval;//a
	Interval interval;//b

	//slab intersection
	interval.Intersect(new_ray, slab_);
	final_interval.CompareInterval(interval);


	//with cylinder//c
	float a = new_ray.D.x() * new_ray.D.x() + new_ray.D.y()*new_ray.D.y();
	float b = 2 * (new_ray.D.x() * new_ray.Q.x() + new_ray.D.y()*new_ray.Q.y());
	float c = new_ray.Q.x() * new_ray.Q.x() + new_ray.Q.y()*new_ray.Q.y() - radius_*radius_;

	float det = b*b - 4.0f * a*c;

	if (det < 0) {

		return false;
	}

	float sqrt_det = sqrtf(det);
	float one_over_2a = 1.0f / (2 * a);
	
	float c0 = (-b - sqrt_det) * one_over_2a;
	float c1 = (-b + sqrt_det) *  one_over_2a;

	if (c0 < 0 && c1 < 0) {

		return false;
	}

	if (c0 > final_interval.t0_) {

		final_interval.t0_ = c0;
		Vector3f point = new_ray.Eval(c0);
		final_interval.N0_ = Vector3f(point.x(), point.y(), 0.0f);
	}

	if (c1 < final_interval.t1_) {

		final_interval.t1_ = c1;
		Vector3f point = new_ray.Eval(c1);
		final_interval.N1_ = Vector3f(point.x(), point.y(), 0.0f);
	}

	if (final_interval.t0_ > final_interval.t1_) {

		return false;

	}else if (final_interval.t0_ < 0 && final_interval.t1_<0) {

		return false;

	}
	/*else if (final_interval.t0_ == 0 && final_interval.t1_ == 0) {

		return false;
	}*/
	else {

		/*if (final_interval.t0_ == 0 || final_interval.t1_ == 0) {

			std::cout << "one of the root is zero: "<< final_interval.t0_ <<" , "<< final_interval.t1_<<" \n";
		}*/
		

		if (final_interval.t0_ < 0 || final_interval.t1_ < 0) {
			id.t = std::max(final_interval.t0_, final_interval.t1_);
		}
		else {
			id.t = std::min(final_interval.t0_, final_interval.t1_);
		}



		//if (final_interval.t0_ * final_interval.t1_ >  0) {

		//	if (final_interval.t0_ < final_interval.t1_) {


		//		id.t = final_interval.t0_;
		//		id.p = ray.Eval(id.t);
		//		id.s = this;
		//		id.n = q_.conjugate()._transformVector(final_interval.N0_);

		//	}
		//	else {

		//		id.t = final_interval.t1_;
		//		id.p = ray.Eval(id.t);
		//		id.s = this;
		//		id.n = q_.conjugate()._transformVector(final_interval.N1_);

		//	}
		//	return true;
		//}
		//else {

		////	std::cout << "product < 0 \n";
		//	if (final_interval.t0_ >= 0) {

		//		id.t = final_interval.t0_;
		//		id.p = ray.Eval(id.t);
		//		id.s = this;
		//		id.n = q_.conjugate()._transformVector(final_interval.N0_);
		//		//id.n.normalize();
		//		return true;
		//	}
		//	else if(final_interval.t1_ >= 0) {

		//		id.t = final_interval.t1_;
		//		id.p = ray.Eval(id.t);
		//		id.s = this;
		//		id.n = q_.conjugate()._transformVector(final_interval.N1_);
		//	//	id.n.normalize();
		//		return true;
		//	}
		//	else {
		//		//std::cout << "On of the root is zero returning \n";
		//		return false;
		//	}
		//}

		//working code


		id.p = ray.Eval(id.t);
		id.s = this;
	
		if (id.t == c0 || interval.t0_) {
			id.n = q_.conjugate()._transformVector(final_interval.N0_);
		}else {

			id.n = q_.conjugate()._transformVector(final_interval.N1_);
		}
		//id.n = final_interval.N0_;
		return true;

		//setting the normals

	/*	Eigen::Vector3f normal(0, 0, 0);
		Eigen::Vector3f qAxisVec = q_._transformVector(axis_);
		Eigen::Vector3f qBasePt = q_._transformVector(base_);
		Ray r2 = ray;
		r2.Q = q_._transformVector(r2.Q);
		r2.D = q_._transformVector(r2.D);


		if (r2.Eval(id.t).z() >= (qBasePt + qAxisVec).z())
			normal.z() = 1;
		else if (r2.Eval(id.t).z() <= (qBasePt).z())
			normal.z() = -1;
		else
		{
			Eigen::Vector3f m(new_ray.Eval(id.t));
			normal = Eigen::Vector3f(m.x(), m.y(), 0);
			normal = q_.inverse()._transformVector(normal);
		}

		id.n = normal;
		return true;*/
	}
}


Cylinder::~Cylinder()
{
}
