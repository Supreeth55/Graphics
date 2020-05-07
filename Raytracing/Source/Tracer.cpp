#include "../Headers/Tracer.h"
#include "../Headers/Sphere.h"
#include "../Headers/Minimizer.h"
#include <random>
#include "../raytrace.h"
std::mt19937_64 RNGen1, RNGen2;
std::uniform_real_distribution<> myrandom1(0.0f, 1.0f), myrandom2(0.0f, 1.0f);
#define RUSSIAN_ROULETTE 0.8f
#define EXPLICIT
#define PI 3.14159265358979323846
#define EPSILON 0.00001

Color Tracer::TraceRay( Ray& ray, std::vector<Shape*> &lights, KdBVH<float, 3, Shape*>& tree_, float deltaTime) {
	C = Color(0.0f, 0.0f, 0.0f);//cOLOR
	W = Color(1.0f, 1.0f, 1.0f);// WEIGHT PONDERATOR
	float p,q,Wmis;
	bool objectMoved = false; //reference to move 1 object
	Vector3f offset = Vector3f(0, 0, 0);
	Vector3f initialPosition = Vector3f(0, 0, 0);

	//INITIAL RAY
	Minimizer m(ray);
	Intersection *P = BVMinimize(tree_, m) == FLT_MAX ? NULL:   &m.id_;// intersection with first object
	
	if (!P)  return C;// no intersection, black
	if (P->s->is_light_)  return P->s->material_->Kd;// is light KD
	
	/////////////SOLUTION #2 
	//// IF THE OBJECT IS MOVING, JUST  APPLY THE OFFSET OF THE OBJECT ON THE RAY AND RECALCULATE 
	//if (P->s->isMoving) {
	//	ray.Q = ray.Q - (static_cast<Sphere*>(P->s)->center2 - static_cast<Sphere*>(P->s)->center) * deltaTime;
	//	m =  Minimizer(ray);
	//	P = BVMinimize(tree_, m) == FLT_MAX ? NULL : &m.id_;// intersection with first object
	//	if (!P)  return C;// no intersection, black
	//	if (P->s->is_light_)  return P->s->material_->Kd;// is light  KD
	//}
		

	Vector3f wi, wo = -ray.D;
	while (myrandom1(RNGen1) <= RUSSIAN_ROULETTE) {
		Vector3f normal = P->n.normalized();
#ifdef EXPLICIT
		Intersection L = SampleLight( lights );// Randomly choose a light 
		p = PdfLight(L) / GeometryFactor( *P, L);//  Probability to hit light ( angular measure )
		q = PdfBrdf(E_ALL, normal, wo, wi, P->s) * RUSSIAN_ROULETTE;//probability of diffuse + reflection + refraction
		Wmis = p * p / (p * p + q * q);
		wi = (L.p - P->p).normalized();
		Ray shadowRay( (P->p - wo * 0.01f) , wi);//Ray goes explicit to light
		Minimizer min(shadowRay);
		Intersection *I =  BVMinimize(tree_, min) == FLT_MAX ? NULL: &min.id_;//Shadow ray Intersection

		if (p > 0 && I  && I->p == L.p) {
			Color f = EvalScattering(E_ALL, normal, wo, wi, P->s, P->t);
			C += W * (f / p) * Wmis * (Color)L.s->get_material()->Kd;
		}
#endif
		float rand = myrandom1(RNGen1);
		ChoiceType choice = rand < P->s->material_->pd ? E_DIFFUSE :rand < P->s->material_->pd + P->s->material_->pr ? E_REFLECTION :E_TRANSMISSION;
		
		wi = SampleBrdf(choice, normal,  wo, P->s);// calculates rebound, based on diffuse, reflection, refraction
		Ray new_ray((P->p + wi * 0.01f), wi);
		Minimizer min1(new_ray);
		BVMinimize(tree_, min1);
		Intersection Q = min1.id_;
		if (!Q.s )   break;



		
		Color f = EvalScattering(choice, normal, wo, wi, P->s, P->t);
		 p = PdfBrdf(E_ALL, normal, wo, wi, P->s) * RUSSIAN_ROULETTE;
		if (p < EPSILON)   break;

		W *= f / p;

		//IMPLICIT LIGHT CONNECTION
		if (Q.s->material_->isLight()) {
			q = PdfLight(Q) / GeometryFactor(*P, Q);//Probability the implicit light could be chosen explicitly
			Wmis = 1;// p * p / (p * p + q * q);
			C += W * Wmis * (Color)Q.s->get_material()->Kd;
			break;
		}
		P->SetIntersection(Q);
		wo = - wi;
	}
	return  C;
}

Color Tracer::EvalScattering(const ChoiceType& choice, Vector3f& normal, Vector3f& wo, Vector3f& wi, Shape*shape, float&t) {
	Color color;
	switch (choice) {
		case E_DIFFUSE: {
			color = EvalBrdf(E_DIFFUSE,normal,wo,wi,shape,t);
		}break;
		case E_REFLECTION: {
			color = EvalBrdf(E_REFLECTION, normal, wo, wi, shape, t);
		}break;
		case E_TRANSMISSION: {
			color = EvalBrdf(E_TRANSMISSION, normal, wo, wi, shape, t);
		}break;
		default: {//ALL
			color =
				EvalBrdf(E_DIFFUSE, normal, wo, wi, shape, t) +
				EvalBrdf(E_REFLECTION, normal, wo, wi, shape, t) +
				EvalBrdf(E_TRANSMISSION, normal, wo, wi, shape, t);
		}break;
	}
	return fabsf(normal.dot(wi)) * color;
}

Color Tracer::EvalBrdf(const ChoiceType& choice, Vector3f& normal, Vector3f& wo, Vector3f& wi , Shape*shape,  float&t) {
	if (choice == E_DIFFUSE) {
		return (shape->get_material()->Kd / PI);
	}
	else if(choice == E_REFLECTION) {
		Vector3f m = (wo + wi).normalized();
		float woDOTn = fabsf(wo.dot(normal));
		float wiDOTn = fabsf(wi.dot(normal));
		float wiDOTm = wi.dot(m);
		float denom = 1.0f/ (4.0f * woDOTn*wiDOTn);

		return (Fresnel(wiDOTm, shape) * PhongD(normal, m, shape->get_material()->alpha) * PhongG(wo, wi, m, normal, shape->get_material()->alpha) * denom);
	}
	else if(choice == E_TRANSMISSION){
		float ior = shape->get_material()->ior_;
		float alpha = shape->get_material()->alpha;
		Vector3f materialKt = shape->get_material()->Kt;

		float n = 0.0f;
		float ni = 1.0f;
		float no = ior;
		float atenX = 1.0f;
		float atenY = 1.0f;
		float atenZ = 1.0f;

		Vector3f attuenation_;
		if (wo.dot(normal) < 0.0f) {
			ni = ior;
			no = 1.0f;
			atenX = powf(2.71f, t * log(materialKt.x()));
			atenY = powf(2.71f, t * log(materialKt.y()));
			atenZ = powf(2.71f, t * log(materialKt.z()));
			attuenation_ = Vector3f(atenX, atenY, atenZ);
		}
		else {
			attuenation_ = Vector3f(atenX, atenY, atenZ);
		}
		n = ni/no;

		Vector3f m = -(wo*ni + wi*no).normalized();
		float woDOTm = wo.dot(m);
		float woDOTn = fabsf(wo.dot(normal));
		float wiDOTn = fabsf(wi.dot(normal));
		float wiDOTm = wi.dot(m);

		float radicand = 1 - (n*n)*(1 - (woDOTm*woDOTm));
		if (radicand < 0) {
			float denom = 1.0f / (4.0f * woDOTn*wiDOTn);
			return (Color)attuenation_* (Fresnel(wiDOTm, shape) * PhongD(normal, m, alpha) * PhongG(wo, wi, m, normal, alpha) * denom);
		}
		else {
			float denom = 1.0f / ( woDOTn*wiDOTn);
			float abs_wi_dot_m = fabsf(wiDOTm);
			float abs_wo_dot_m = fabsf(woDOTm);
			float no_sqrd = no*no;
			float numo = abs_wi_dot_m*abs_wo_dot_m*no_sqrd;
			float deno = (no*(wiDOTm)) + (ni*(woDOTm));
			float deno_sqrd = deno*deno;
			float res = numo / deno_sqrd;
			return  (Color)attuenation_* ((Color(1,1,1) - Fresnel(wiDOTm, shape)) * PhongD(normal, m, alpha) * PhongG(wo, wi, m, normal, alpha) * denom * res);
		}
	}
	return Color(0, 0, 0);//This never happens is to avoid a warning
}

//Returns the Wi. acording to normal. wo and  surface specig
Vector3f  Tracer::SampleBrdf(const ChoiceType& choice,Vector3f& normal,  Vector3f& wo , Shape* shape) {
	Vector3f wi;
	float alpha = shape->get_material()->alpha;
	if (choice == E_DIFFUSE) {
		wi = SampleLobe(normal  ,  sqrtf(myrandom2(RNGen2))  ,   2 * PI * (myrandom2(RNGen2)));
	}
	else if(choice == E_REFLECTION){
		Vector3f m = SampleLobe(normal  ,  pow(myrandom2(RNGen2), 1.0f / (alpha + 1.0f))  ,  2 * PI * (myrandom2(RNGen2)));
		wi = 2.0f * (wo.dot(m)) * m - wo;
	}
	else if(choice == E_TRANSMISSION){
		Vector3f m = SampleLobe(normal, pow(myrandom2(RNGen2), 1.0f / (alpha + 1.0f)),    2 * PI * (myrandom2(RNGen2)));
		float woDOTn = wo.dot(normal);
		float ior = shape->get_material()->ior_;
		float ni = woDOTn<0? ior:  1.0f;
		float no = woDOTn<0? 1: ior;
		float n =  ni / no;
		float woDOTm = wo.dot(m);
		float radicand = 1.0f - ( n * n ) * (1 - (woDOTm*woDOTm));
		if (radicand < 0) 
			wi = 2.0f * woDOTm * m - wo;
		else 
			wi =  (((n * (woDOTm))-(woDOTn >= 0 ? 1 : -1 ) * sqrtf(radicand)) * m) - n * wo;
	}
	return wi.normalized();
}

Vector3f Tracer::SampleLobe(Vector3f normal, float r1, float r2 ) {
	float s = sqrtf(1 - (r1*r1));
	Quaternionf q = Quaternionf::FromTwoVectors(Vector3f::UnitZ(), normal);
	return q._transformVector(Vector3f(s * cosf(r2), s * sinf(r2), r1) );
} 

float Tracer::PdfBrdf(const ChoiceType& choice, Vector3f &normal, Vector3f &wo, Vector3f &wi , Shape* shape) {
	float alpha = shape->get_material()->alpha;
	float pd = shape->get_material()->pd;
	float pr = shape->get_material()->pr;
	float pt = shape->get_material()->pt;
	     if(choice == E_DIFFUSE) {
		return pd * (fabsf(wi.dot(normal)) / PI);
	}
	else if(choice == E_REFLECTION){
		Vector3f m = (wo+wi).normalized();
		return  pr * PhongD(normal, m,  alpha) * (fabsf(m.dot(normal))) / (4.0f * fabsf(wi.dot(m)));
	}
	else if(choice == E_TRANSMISSION) {
		float ior = shape->get_material()->ior_;
		float woDOTn = wo.dot(normal);
		float ni = woDOTn < 0 ? ior : 1.0f;
		float no = woDOTn < 0 ? 1 : ior;
		float n = ni / no;
		Vector3f m = -(wo*ni + wi*no).normalized();
		float woDOTm = wo.dot(m);
		float radicand = 1 - (n*n)*(1 - (woDOTm*woDOTm));
		
		if (radicand < 0.0f) { //total internal reflection
			m = (wo + wi).normalized(); 
			return  pt * PhongD(normal, m, alpha) * (fabsf(m.dot(normal))) / (4.0f * fabsf(wi.dot(m)));
		}
		else {
			float deno = (no * (wi.dot(m))) + (ni * (woDOTm));
			return (pt * PhongD(normal, m, alpha) * fabsf(m.dot(normal)) * no * no * fabsf(wi.dot(m)) / deno / deno);
		}
	}
	else  {
		return PdfBrdf(E_DIFFUSE, normal, wo, wi, shape) + PdfBrdf(E_REFLECTION, normal, wo, wi, shape) + PdfBrdf(E_TRANSMISSION, normal, wo, wi, shape);
	}
}

float Tracer::GeometryFactor(Intersection& A , Intersection& B) {
	Vector3f d = A.p - B.p;
	float a_normal_dot_d = A.n.dot(d);
	float b_normal_dot_d = B.n.dot(d);
	float d_normal_dot_d = d.dot(d);
	float d_dot_sqr = d_normal_dot_d*d_normal_dot_d;
	return fabsf((a_normal_dot_d*b_normal_dot_d) / d_dot_sqr);
}

Intersection Tracer::SampleLight( std::vector<Shape*> &lights) {
	
	int index = myrandom1(RNGen1) * lights.size();
	return SampleSphere((Sphere*)lights[index]);
	
	/*Sphere *sph = (Sphere*)lights[0];
	Intersection result;
	SampleSphere(result, sph );
	return result;*/
}

Intersection Tracer::SampleSphere( Sphere *sph ){
	Intersection result;
	float z = 2 * myrandom1(RNGen1) - 1.0f;
	float r = sqrtf(1 - z * z);
	float a = 2 * PI * myrandom1(RNGen1);
	result.n = Vector3f(r * cos(a), r * sin(a), z);
	result.p = sph->center + result.n * sph->radius_;
	result.s = sph;
	return result;
}

float Tracer::PdfLight(Intersection& result) {
	return (1.0f / result.s->get_area());
}

Color Tracer::Fresnel(const float& wiDOTm,Shape* shp) {//TODO since  when wiDOTm equal to LDOTH
	return shp->material_->Ks + (Vector3f(1, 1, 1) - shp->material_->Ks) * (pow(1-fabs(wiDOTm) , 5));
}

float Tracer::PhongD(Vector3f& normal , Vector3f& m, float& alpha) {
	return ((m.dot(normal))>0?1:0) * ((alpha + 2.0f) / (2.0f * PI)) * pow((m.dot(normal)), alpha);
}

float Tracer::PhongG(Vector3f&wo, Vector3f& wi, Vector3f& m, Vector3f& normal , float& alpha) {
	return PhondG1(wi, m, normal, alpha) * PhondG1(wo, m, normal, alpha);
}

float Tracer::PhondG1(Vector3f& v, Vector3f&m , Vector3f& normal, float& alpha) {
	float vDOTn = v.dot(normal);
	if (vDOTn >1.0f)   return 1.0f;

	float tanTheta = sqrtf(1.0f - (vDOTn * vDOTn)) / vDOTn;
	if (tanTheta == 0.0f)   return 1.0f;

	float a = sqrtf((alpha / 2.0f) + 1.0f) / fabs(tanTheta);

	if (a < 1.6f) 
		return ((v.dot(m) / vDOTn) > 0 ? 1.0f : 0) * (3.535f * a + 2.181f * a * a) / (1.0f + 2.276f * a + 2.577f * a * a);
	else 
		return ((v.dot(m) / vDOTn) > 0 ? 1.0f : 0) ;
}