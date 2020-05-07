#pragma 
#include "Shape.h"
#include "Slab.h"

class Box : public Shape
{
private:

	Vector3f corner_;
	Vector3f diagonal_;
	Slab slab1_;
	Slab slab2_;
	Slab slab3_;
	
public:
	Box(Vector3f corner, Vector3f diagonal , Material* material);

	virtual float get_area() { return 0.0f; };
	virtual bool Intersect(Ray&ray, Intersection& id);


	virtual Bbox bbox() { return box_; }
	~Box();
};

