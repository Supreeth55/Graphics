#include "../Headers/Ray.h"


Ray::Ray(Vector3f start_point, Vector3f direction) {

	Q = start_point;
	D = direction.normalized();
}


Vector3f Ray::Eval(float intersection_t)
{
	return Vector3f(Q + D*intersection_t);
}

Ray::~Ray()
{
}
