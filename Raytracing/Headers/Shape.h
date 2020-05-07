#pragma once
#include "Ray.h"
#include "Intersection.h"
#include "../geom.h"
#include "../raytrace.h"

typedef enum {

	E_NONE,
	E_BOX,
	E_CYLINDER,
	E_SPHERE,
	E_TRIANGLE,
	E_TOTAL

}ShapeType;

class Material;
class Shape 
{
public:
	Shape(ShapeType type) :type_(type) , is_light_(false), isMoving (false) {}

	virtual float get_area() = 0;
	virtual bool Intersect(Ray& ray , Intersection& id) = 0;
	virtual Bbox bbox() = 0;
	Material* get_material() { return material_; };

	ShapeType get_type() { return type_; };

	inline bool get_is_light() { return is_light_; };
	inline void set_is_light(bool is_light) { is_light_ = is_light; }
	virtual ~Shape();

	ShapeType type_;
	Material* material_;
	bool is_light_;
	bool isMoving ;
	Bbox box_;
};

