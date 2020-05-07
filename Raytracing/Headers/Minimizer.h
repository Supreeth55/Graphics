#pragma once
#include <Eigen/StdVector>
#include <Eigen_unsupported/Eigen/BVH>
#include "Ray.h"
#include "../geom.h"
#include <algorithm>
#include "Intersection.h"

class Shape;
class Minimizer
{
public:
	typedef float Scalar; // KdBVH needs Minimizer::Scalar defined
	Ray* ray;
	float t_;
	Shape* shape_;
	Intersection id_;
	bool hit_;
	// Stuff to track the minimal t and its intersection info
	
	// Constructor
	Minimizer(Ray& r);
	Minimizer& operator=(const Minimizer& m);
	Minimizer(const Minimizer& m);


	// Called by BVMinimize to intersect the ray with a Shape. This
	// should return the intersection t, but should also track
	// the minimum t and it's corresponding intersection info.
	// Return INF to indicate no intersection.
	float minimumOnObject(Shape* obj);

	// Called by BVMinimize to intersect the ray with a Box3d and
	// returns the t value. This should be similar to the already
	// written box (3 slab) intersection. (The difference begin a
	// distance of zero should be returned if the ray starts within the bbox.)
	// Return INF to indicate no intersection.
	float minimumOnVolume(const Bbox& box);

	//for intersections with box
	float minimumOnObjectVolume(Shape* shp);

};