int passes = 4096;
int  eachnSteps = 500;// if 0, doesnt draw intermedial images 
bool deepthOfField = false;

#include <vector>
#include <windows.h>
#include <cstdlib>
#include <limits>
#include <crtdbg.h>
#include "geom.h"
#include "raytrace.h"
#include "realtime.h"
#include "stb_image.h"
#include <random>// A good  *thread-safe* random number generator.
#include "Headers\Shape.h"
#include "Headers\Box.h"
#include "Headers\Cylinder.h"
#include "Headers\Sphere.h"
#include "Headers\Camera.h"
#include "Headers\Triangle.h"
#include "Headers\Tracer.h"
#include "Headers\Minimizer.h"
#define STB_IMAGE_IMPLEMENTATION
std::mt19937_64 RNGen;	
std::uniform_real_distribution<> myrandom(0.0, 1.0);// Call myrandom(RNGen) to get a uniformly distributed random number in [0,1].
std::vector<Shape*> list_;
std::vector<Shape*> listInMovement_;
std::vector<Shape*> lights_;
Camera* camera_;
#include "rgbe.h"
void WriteHdrImage(const std::string outName, const int width, const int height, Color* image, int pass_num = 1)
{
	float* data = new float[width*height * 3];	// Turn image from a 2D-bottom-up array of Vector3D to an top-down-array of floats
	float* dp = data;
	for (int y = height - 1; y >= 0; --y) {
		for (int x = 0; x<width; ++x) {
			Color pixel = image[y*width + x] / (float)(pass_num);
			*dp++ = pixel[0];
			*dp++ = pixel[1];
			*dp++ = pixel[2];
		}
	}

	// Write image to file in HDR (a.k.a RADIANCE) format
	rgbe_header_info info;
	char errbuf[100] = { 0 };
	FILE* fp = fopen(outName.c_str(), "wb");
	info.valid = false;
	int r = RGBE_WriteHeader(fp, width, height, &info, errbuf);
	if (r != RGBE_RETURN_SUCCESS)
		printf("error: %s\n", errbuf);

	r = RGBE_WritePixels_RLE(fp, data, width, height, errbuf);
	if (r != RGBE_RETURN_SUCCESS)
		printf("error: %s\n", errbuf);
	fclose(fp);

	delete data;
}

std::string hrsMinsSecs(int _secs) {
	std::string result;
	int forHours = _secs / 3600,
		remainder = _secs % 3600,
		forMinutes = remainder / 60,
		forSeconds = remainder % 60;

	return (forHours > 0 ?  std::to_string(forHours) + "h" : "") +
		(forMinutes > 0 ?  std::to_string(forMinutes) + "m" : "") +
		(forSeconds > 0 ? std::to_string(forSeconds)+ "s"  : "")
		;
}

Bbox bounding_box(Shape* shape) {
	return shape->bbox();
}

Scene::Scene() 
{ 
}

Scene::~Scene()
{
	delete camera_;
	for (auto shape : list_) {
		delete shape;
	}
	list_.clear();
}

void Scene::triangleMesh(MeshData* mesh) 
{ 
	for (auto &indices : mesh->triangles) {
		Vector3f v0 = mesh->vertices[indices[0]].pnt;
		Vector3f v1 = mesh->vertices[indices[1]].pnt;
		Vector3f v2 = mesh->vertices[indices[2]].pnt;
		Vector3f n0 = mesh->vertices[indices[0]].nrm;
		Vector3f n1 = mesh->vertices[indices[1]].nrm;
		Vector3f n2 = mesh->vertices[indices[2]].nrm;

		auto triangle = new Triangle(mesh->mat, v0, v1, v2, n0, n1, n2);
		list_.push_back(triangle);
	}
}

Quaternionf Orientation(unsigned int i, const std::vector<std::string>& strings, const std::vector<float>& f)
{
    Quaternionf q(1,0,0,0); // Unit quaternion
    while (i<strings.size()) {
        std::string c = strings[i++];
        if (c == "x")  
            q *= angleAxis(f[i++]*Radians, Vector3f::UnitX());
        else if (c == "y")  
            q *= angleAxis(f[i++]*Radians, Vector3f::UnitY());
        else if (c == "z")  
            q *= angleAxis(f[i++]*Radians, Vector3f::UnitZ());
        else if (c == "q")  {
            q *= Quaternionf(f[i+0], f[i+1], f[i+2], f[i+3]);
            i+=4; }
        else if (c == "a")  {
            q *= angleAxis(f[i+0]*Radians, Vector3f(f[i+1], f[i+2], f[i+3]).normalized());
            i+=4; } }
    return q;
}

////////////////////////////////////////////////////////////////////////
// Material: encapsulates surface properties
void Material::setTexture(const std::string path)
{
	/*int width, height, n;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* image = stbi_load(path.c_str(), &width, &height, &n, 0);

    // Realtime code below:  This sends the texture in *image to the graphics card.
    // The raytracer will not use this code (nor any features of OpenGL nor the graphics card).
    glGenTextures(1, &texid);
    glBindTexture(GL_TEXTURE_2D, texid);
    glTexImage2D(GL_TEXTURE_2D, 0, n, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 100);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (int)GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (int)GL_LINEAR_MIPMAP_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    stbi_image_free(image);*/
}

void Material::CalculateProbablities() {
	totalK = Kd.norm() + Ks.norm() + Kt.norm();
	pd = Kd.norm() / totalK;
	pr = Ks.norm() / totalK;
	pt = Kt.norm() / totalK;//check if this has to be Kt
}

void Scene::Command(const std::vector<std::string>& strings,
                    const std::vector<float>& f)
{
    if (strings.size() == 0) return;
    std::string c = strings[0];
    
    if (c == "screen") {
        // syntax: screen width height
      //  realtime->setScreen(int(f[1]),int(f[2]));
        width = int(f[1]);
        height = int(f[2]); }

    else if (c == "camera") {
        // syntax: camera x y z   ry   <orientation spec>
        // Eye position (x,y,z),  view orientation (qw qx qy qz),  frustum height ratio ry
     //   realtime->setCamera(Vector3f(f[1],f[2],f[3]), Orientation(5,strings,f), f[4]); 

		camera_ = new Camera(Vector3f(f[1], f[2], f[3]), Orientation(5, strings, f), f[4], (float)width, (float)height);
	}

    else if (c == "ambient") {
        // syntax: ambient r g b
        // Sets the ambient color.  Note: This parameter is temporary.
        // It will be ignored once your raytracer becomes capable of
        // accurately *calculating* the true ambient light.
      //  realtime->setAmbient(Vector3f(f[1], f[2], f[3])); 
	}

    else if (c == "brdf")  {
        // syntax: brdf  r g b   r g b  alpha
        // later:  brdf  r g b   r g b  alpha  r g b ior
        // First rgb is Diffuse reflection, second is specular reflection.
        // third is beer's law transmission followed by index of refraction.
        // Creates a Material instance to be picked up by successive shapes

		if(f.size()==12)
			currentMat = new Material(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7], Vector3f(f[8], f[9], f[10]) , f[11]);
		else 
			currentMat = new Material(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7]);
	}

    else if (c == "light") {
        // syntax: light  r g b   
        // The rgb is the emission of the light
        // Creates a Material instance to be picked up by successive shapes
        currentMat = new Light(Vector3f(f[1], f[2], f[3])); }
   
    else if (c == "sphere") {
        // syntax: sphere x y z   r
        // Creates a Shape instance for a sphere defined by a center and radius
     //   realtime->sphere(Vector3f(f[1], f[2], f[3]), f[4], currentMat); 
		Shape* shape;
		if (f.size() == 8)
			shape = new Sphere(Vector3f(f[1], f[2], f[3]), currentMat, f[4], Vector3f(f[5], f[6], f[7]));
		else 
			shape = new Sphere(Vector3f(f[1], f[2], f[3]),  currentMat , f[4]);
		list_.push_back(shape);

		if (currentMat->isLight()) {
			shape->set_is_light(true);
			lights_.push_back(shape);
		}
	}

    else if (c == "box") {
        // syntax: box bx by bz   dx dy dz
        // Creates a Shape instance for a box defined by a corner point and diagonal vector
    //    realtime->box(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), currentMat); 
		Shape* shape = new Box(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), currentMat);
		list_.push_back(shape);
	}

    else if (c == "cylinder") {
        // syntax: cylinder bx by bz   ax ay az  r
        // Creates a Shape instance for a cylinder defined by a base point, axis vector, and radius
      //  realtime->cylinder(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7], currentMat); 

		Shape* shape = new Cylinder(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7], currentMat);
		list_.push_back(shape);
	}

    else if (c == "capsule") {
        // syntax: cylinder bx by bz   ax ay az  r
        // Creates a Shape instance for a cylinder defined by a base point, axis vector, and radius
    //    realtime->cylinder(Vector3f(f[1], f[2], f[3]), Vector3f(f[4], f[5], f[6]), f[7], currentMat); 
	}

    else if (c == "mesh") {
        // syntax: mesh   filename   tx ty tz   s   <orientation>
        // Creates many Shape instances (one per triangle) by reading
        // model(s) from filename. All triangles are rotated by a
        // quaternion (qw qx qy qz), uniformly scaled by s, and
        // translated by (tx ty tz) .
        Matrix4f modelTr = translate(Vector3f(f[2],f[3],f[4]))
                          *scale(Vector3f(f[5],f[5],f[5]))
                          *toMat4(Orientation(6,strings,f));
        ReadAssimpFile(strings[1], modelTr);  
	}
    else {
        fprintf(stderr, "\n*********************************************\n");
        fprintf(stderr, "* Unknown command: %s\n", c.c_str());
        fprintf(stderr, "*********************************************\n\n");
    }
}

void Scene::TraceImage(Color* image)
{
	Vector3f X = camera_->rx_* camera_->orientation_._transformVector(Vector3f::UnitX());
	Vector3f Y = camera_->ry_* camera_->orientation_._transformVector(Vector3f::UnitY());
	Vector3f Z = -1.0f * camera_->orientation_._transformVector(Vector3f::UnitZ());
	//KdBVH<float, 3, Shape*> Tree(list_.begin(), list_.end());

	clock_t begin = clock();
	printf( "Pass :" );
	int cntr = 0;
	while (cntr++ < passes) {
		printf(  "%d,", cntr);
		
		//MOTION BLUR
		 
		float deltaTime = 1.0f * cntr / (1.0f * passes);
		deltaTime = deltaTime < 0.2f ? 1 :	deltaTime;
		listInMovement_.clear();
		for(auto shape : list_) {
			if (shape->isMoving) {
				Vector3f offset = (static_cast<Sphere*>(shape)->center2 - static_cast<Sphere*>(shape)->center) * deltaTime;
				
				Shape * shape2 = new Sphere(
					static_cast<Sphere*>(shape)->center + offset,
					static_cast<Sphere*>(shape)->material_, 
					static_cast<Sphere*>(shape)->radius_
				);
				listInMovement_.push_back(shape2);
			}
			else
				listInMovement_.push_back(shape);
		}

		//printf("list size %d/n", listInMovement_.size());
		KdBVH<float, 3, Shape*> Tree(listInMovement_.begin(), listInMovement_.end());

		#pragma omp parallel for schedule(dynamic, 1) // Magic: Multi-thread  
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {	
				
				float D = 5.0f, 
					W = 0.5,
					r = W * sqrt( myrandom(RNGen)),
				theta = 2.0f* PI * r * myrandom(RNGen),
					rx =  r * cos(theta),
					ry =  r * sin(theta);
				
				float dx = 2.0f * ((x + (float)(myrandom(RNGen))) / (float)width) - 1.0f;
				float dy = 2.0f * ((y + (float)(myrandom(RNGen))) / (float)height) - 1.0f;
				Vector3f direction = dx*X + dy*Y + Z;// no deepth of field
				Vector3f direction2 = ((D * dx - rx) * X + (D * dy - ry) * Y + D * Z).normalized();
				
				Ray ray(
					deepthOfField ? camera_->eye_ + rx *X + ry * Y :  camera_->eye_,
					deepthOfField ? direction2 : direction
				);
				Tracer tracer;
				image[y*width + x] += tracer.TraceRay(ray, lights_, Tree, deltaTime);
			}// width  
		}//  height  
		if (eachnSteps != 0 &&  cntr % eachnSteps == 0) {
			std::string new_filename = "testscene_p"+ std::to_string(cntr)  +  "_t"  +  hrsMinsSecs(int(clock() - begin) / CLOCKS_PER_SEC) + ".hdr";
			WriteHdrImage(new_filename, width, height, image, cntr);
			printf("Sample  :%.2s n", new_filename);
			}
		
	}//passes
	printf("time :%.2f sec \n\n", double(clock() - begin) / CLOCKS_PER_SEC);
}