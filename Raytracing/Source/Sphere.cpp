#include "../Headers/Sphere.h"



Sphere::Sphere(Vector3f center, Material* material, float radius) :Shape(E_SPHERE), center(center),radius_(radius) {
	material_ = material;
	Vector3f c1(center + Vector3f(-radius, -radius, -radius));
	Vector3f c2(center + Vector3f(radius, radius, radius));
	box_ = Bbox(c1, c2);
}

Sphere::Sphere(Vector3f center, Material* material, float radius, Vector3f _center2) :Shape(E_SPHERE), center(center), radius_(radius) ,center2(_center2){
	material_ = material;

	Vector3f c1(center + Vector3f(-radius, -radius, -radius));
	Vector3f c2(center + Vector3f(radius, radius, radius));

	Vector3f c3(_center2 + Vector3f(-radius, -radius, -radius));
	Vector3f c4(_center2 + Vector3f(radius, radius, radius));

	// Calculate the bounding box that contains initial and final position
	Vector3f c5(std::min(c1.x(),c3.x()), std::min(c1.y(), c3.y()), std::min(c1.z(), c3.z()) );
	Vector3f c6(std::max(c2.x(), c4.x()), std::max(c2.y(), c4.y()), std::max(c2.z(), c4.z()));
	isMoving = true;
	 
	box_ = Bbox(c5, c6);
}

bool Sphere::Intersect(Ray&ray, Intersection& id) {
	float epsilon = 0.00001f;
	Vector3f q = ray.Q - center;
	float qDOTd = q.dot(ray.D);
	float qDOTq = q.dot(q);
	float fvar = sqrt(qDOTd * qDOTd - qDOTq + radius_ * radius_);

	float tempo;
	if (-qDOTd - fvar > epsilon) 
		tempo = -qDOTd - fvar;
	else if (-qDOTd + fvar > epsilon) 
		tempo = -qDOTd + fvar;
	else return false;
		/*id = Intersection(ray.Eval(id.t), (id.p - center).normalized(), this, (t0 < 0 || t1 < 0) ? std::max(t0, t1) : std::min(t0, t1));*/
		id.p = ray.Eval(tempo);
		id.n = (id.p - center).normalized();
		id.s = this;
		id.t = tempo;
		return true;
}

float Sphere::get_area() {
	return (float)(4.0f * M_PI*(radius_*radius_));
}

Sphere::~Sphere(){

}
