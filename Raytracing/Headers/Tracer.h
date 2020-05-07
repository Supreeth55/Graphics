#pragma once
#include <vector>
#include "Shape.h"
#include "Intersection.h"
class Sphere;
class Tracer
{
	Color C;
	Color W;
	Intersection  SampleSphere(  Sphere * sphere);// &center, float& radius);
	Intersection  SampleLight( std::vector<Shape*> &lights);
	Vector3f  SampleBrdf(const ChoiceType& choice,Vector3f& normal,  Vector3f& wo , Shape*);
	Vector3f  SampleLobe(Vector3f normal, float r1, float r2);
	float GeometryFactor(Intersection& A, Intersection& B);
	float PdfLight(Intersection& result);
	float PdfBrdf(const ChoiceType& choice,Vector3f &normal, Vector3f &wo, Vector3f &wi, Shape*);
	float PhongD(Vector3f& normal, Vector3f& m, float& alpha);
	float PhongG(Vector3f&wo, Vector3f& wi, Vector3f& m, Vector3f& N, float& alpha);
	float PhondG1(Vector3f& v, Vector3f&m, Vector3f& N, float& alpha);
	Color EvalScattering(const ChoiceType& choice, Vector3f& normal, Vector3f& wi, Vector3f& wo, Shape*, float&t);
	Color EvalBrdf(const ChoiceType& choice,Vector3f& normal, Vector3f& wi, Vector3f& wo , Shape*, float&t);
	Color Fresnel(const float& d, Shape* shp);

public:
	Tracer() {};
	~Tracer() {};
	Color TraceRay( Ray& ray, std::vector<Shape*> &lights, KdBVH<float, 3, Shape*>& tree_, float deltaTime);
};

