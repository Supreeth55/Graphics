#include "../Headers/Triangle.h"



Triangle::Triangle(Material* mat, Vector3f v0, Vector3f v1, Vector3f v2, Vector3f n0, Vector3f n1, Vector3f n2) :Shape(E_TRIANGLE), v0_(v0), v1_(v1), v2_(v2), n0_(n0), n1_(n1), n2_(n2) {

	material_ = mat;
	Vector3f c1(std::min({ v0.x(),v1.x(),v2.x() }), std::min({ v0.y(),v1.y(),v2.y() }), std::min({ v0.z(),v1.z(),v2.z() }));
	Vector3f c2(std::max({ v0.x(),v1.x(),v2.x() }), std::max({ v0.y(),v1.y(),v2.y() }), std::max({ v0.z(),v1.z(),v2.z() }));

	box_ = Bbox(c1, c2);

	e1_ = v1_ - v0_;
	e2_ = v2_ - v0_;
	normal_ = e1_.cross(e2_);
}

bool Triangle::Intersect(Ray& ray, Intersection& intersection) {


	Vector3f p = ray.D.cross(e2_);

	float d = p.dot(e1_);
	if (d == 0) {

		return false;
	}

	float o_over_d = 1.0f / d;

	Vector3f s = ray.Q - v0_;
	float u = p.dot(s);
	u = u *o_over_d;

	if (u < 0 || u>1) {

		return false;
	}

	Vector3f q = s.cross(e1_);
	float v = ray.D.dot(q);
	v = v*o_over_d;

	if (v < 0 || (u + v) > 1) {

		return false;
	}

	float t = e2_.dot(q);

	t = t *o_over_d;

	if (t < 0) {

		return false;
	}

	intersection.s = this;
	intersection.t = t;
	intersection.p = ray.Eval(t);
	intersection.n = normal_;// (1 - u - v)*n0_ + u * n1_ + v * n2_;
	//intersection.n.normalize();
	return true;
}