#include "../Headers/Camera.h"



Camera::Camera(Vector3f eye, Quaternionf orientation, float ry,float width , float height)
{
	eye_ = eye;
	orientation_ = orientation;
	ry_ = ry;
	rx_ = ry* (width / height);
}


Camera::~Camera()
{
}
