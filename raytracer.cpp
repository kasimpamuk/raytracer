#include <iostream>
#include <cmath>

#include "parser.h"
#include "ppm.h"

using namespace parser;

Scene scene;
Vec3i BG_COLOR;
float ShadowRayEpsilon;
int MaxRecursionDepth;
Vec3f AmbientLight;
Vec3f** mesh_normals;
Vec3f* triangle_normals;

typedef struct{
    Vec3f o,d;
} Ray;

Vec3f multS(Vec3f vec, float f){  //Scalar multiplication
    Vec3f res;
    res.x = f*vec.x;
    res.y = f*vec.y;
    res.z = f*vec.z;
    return res;
}

Vec3f add(Vec3f vec1, Vec3f vec2){
    Vec3f res;
    res.x = vec1.x+vec2.x;
    res.y = vec1.y+vec2.y;
    res.z = vec1.z+vec2.z;
    return res;
}

float dotProd(Vec3f vec1, Vec3f vec2){
    return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
}

Vec3f crossProd(Vec3f vec1, Vec3f vec2){
    Vec3f res;
    res.x = vec1.y*vec2.z - vec1.z*vec2.y;
    res.y = vec1.z*vec2.x - vec1.x*vec2.z;
    res.z = vec1.x*vec2.y - vec1.y*vec2.x;
    return res;
}

Vec3f normalize(Vec3f vec){
    return multS(vec,1/sqrt(dotProd(vec,vec)));
}

float determinant(Vec3f col1, Vec3f col2, Vec3f col3)
{
    float res =  col1.x*(col2.y*col3.z - col2.z*col3.y)
               + col1.y*(col2.z*col3.x - col2.x*col3.z)
               + col1.z*(col2.x*col3.y - col2.y*col3.x);
    return res;
}

void normal_calculator(){
    mesh_normals = new Vec3f*[scene.meshes.size()];
    for (int mesh = 0; mesh < scene.meshes.size(); mesh++) {

        mesh_normals[mesh] = new Vec3f[scene.meshes[mesh].faces.size()];

        for (int face = 0; face < scene.meshes[mesh].faces.size(); face++) {
            Vec3f a = scene.vertex_data[scene.meshes[mesh].faces[face].v0_id-1];
            Vec3f b = scene.vertex_data[scene.meshes[mesh].faces[face].v1_id-1];
            Vec3f c = scene.vertex_data[scene.meshes[mesh].faces[face].v2_id-1];

            mesh_normals[mesh][face] = normalize(crossProd(add(b,multS(a,-1)), add(c,multS(a,-1))));
        }
    }

    triangle_normals = new Vec3f[scene.triangles.size()];
    for (int triangle = 0; triangle < scene.triangles.size(); triangle++) {
        Vec3f a = scene.vertex_data[scene.triangles[triangle].indices.v0_id-1];
        Vec3f b = scene.vertex_data[scene.triangles[triangle].indices.v1_id-1];
        Vec3f c = scene.vertex_data[scene.triangles[triangle].indices.v2_id-1];

        triangle_normals[triangle] = normalize(crossProd(add(b,multS(a,-1)), add(c,multS(a,-1))));
    }

}

float intersect_sphere(Ray ray, Sphere sphere){
    Vec3f m = add(ray.o,multS(scene.vertex_data[sphere.center_vertex_id-1],-1));

    float A = dotProd(ray.d,ray.d);
    float B = 2*dotProd(ray.d,m);
    float C = dotProd(m,m) - sphere.radius*sphere.radius;

    float delta = B*B - 4*A*C;

    if (delta < 0)
        return -1.0F;
    else
        return (-B - sqrt(delta))/(2*A);
}

float intersect_face(Ray ray, Face face){
    Vec3f v1 = scene.vertex_data[face.v0_id-1];
    Vec3f v2 = scene.vertex_data[face.v1_id-1];
    Vec3f v3 = scene.vertex_data[face.v2_id-1];

    Vec3f column1 = add(v1,multS(v2,-1));
    Vec3f column2 = add(v1,multS(v3,-1));
    Vec3f column3 = ray.d;
    Vec3f col_beta = add(v1, multS(ray.o,-1));

    float A_det = determinant(column1, column2, column3);
    float B_det = determinant(col_beta, column2, column3) / A_det;
    float G_det = determinant(column1, col_beta, column3) / A_det;
    float t= determinant(column1, column2, col_beta) / A_det;

    if (B_det + G_det <= 1 && B_det >= 0 && G_det >= 0)
        return t;
    else
        return -1;
}

void ray_intersect(Ray ray, float& min_t, int& obj_no, int& object_type, int& face_no){
    object_type = -1; // -1:none, 1:mesh, 2:triangle, 3:sphere
    obj_no = -1, face_no = -1;
    float temp = 100000000.0F;
    for(int obj = 0; obj < scene.meshes.size(); obj++){
        for(int face = 0; face<scene.meshes[obj].faces.size(); face++){
            temp = intersect_face(ray, scene.meshes[obj].faces[face]);
            if(temp < min_t && temp>0){
                min_t = temp;
                object_type = 1;
                obj_no = obj;
                face_no = face;
            }
        }
    }

    for(int obj = 0; obj < scene.triangles.size(); obj++){
        temp = intersect_face(ray, scene.triangles[obj].indices);
        if(temp < min_t && temp>0){
            min_t = temp;
            object_type = 2;
            obj_no = obj;
        }
    }

    for(int obj = 0; obj < scene.spheres.size(); obj++){
        temp = intersect_sphere(ray, scene.spheres[obj]);
        if(temp < min_t && temp>0){
            min_t = temp;
            object_type = 3;
            obj_no = obj;
        }
    }
}   

Vec3i getPixelColor(Ray ray, int recursion){
    float min_t = 100000000.0F;
    int object_type, obj_no, face_no;
   
    ray_intersect(ray, min_t, obj_no, object_type, face_no);
    
    Vec3i color = {0,0,0};

    if(obj_no == -1){
        if(recursion) return color;

        color = BG_COLOR;
    }
    else if(min_t < 100000000.0F && obj_no >= 0){
        Material material;
        Vec3f normal;
        Vec3f point = add(ray.o, multS(ray.d,min_t));
        Vec3f o_vec = normalize(add(ray.o, multS(point,-1)));

        for(int light = 0; light < scene.point_lights.size(); light++){
            if(object_type == 1){
                material = scene.materials[scene.meshes[obj_no].material_id-1];
                normal = mesh_normals[obj_no][face_no];
            }
            if(object_type == 2){
                material = scene.materials[scene.triangles[obj_no].material_id-1];
                normal = triangle_normals[obj_no];
            }
            if(object_type == 3){
                material = scene.materials[scene.spheres[obj_no].material_id-1];
                normal = normalize((add(point,multS(scene.vertex_data[scene.spheres[obj_no].center_vertex_id-1],-1))));
            }

            Vec3f l_vec = normalize(add(scene.point_lights[light].position,multS(point,-1)));
            Vec3f intensity = scene.point_lights[light].intensity;
            Vec3f half_vector = normalize(add(o_vec,l_vec));

            float cos_theta = std::max(0.0F, dotProd(normal, l_vec));
            float cos_alpha = std::max(0.0F, dotProd(normal, half_vector)); 

            Vec3f tmp = scene.point_lights[light].position;
            double ir = 1/(pow(tmp.x-point.x,2)+pow(tmp.y-point.y,2)+pow(tmp.z-point.z,2));

            color.x +=  material.diffuse.x*cos_theta*ir*intensity.x;
            color.y +=  material.diffuse.y*cos_theta*ir*intensity.y;
            color.z +=  material.diffuse.z*cos_theta*ir*intensity.z;

            if(cos_theta >= 0.0F){
                color.x += material.specular.x * pow(cos_alpha, material.phong_exponent) *ir* intensity.x;
                color.y += material.specular.y * pow(cos_alpha, material.phong_exponent) *ir* intensity.y;
                color.z += material.specular.z * pow(cos_alpha, material.phong_exponent) * ir*intensity.z;
            }
        }
            

        if (material.is_mirror && recursion < scene.max_recursion_depth) {
            Ray mirror_ray;
            mirror_ray.d = add(multS(normal,(2 * dotProd(normal, o_vec))),  o_vec);
            mirror_ray.o = add(point,multS(normal,scene.shadow_ray_epsilon));
            Vec3i mirror_color = getPixelColor(mirror_ray, recursion+1);

            color.x += mirror_color.x*material.mirror.x;
            color.y += mirror_color.y*material.mirror.y;
            color.z += mirror_color.z*material.mirror.z;
        }

        color.x += scene.ambient_light.x*material.ambient.x;
        color.y += scene.ambient_light.y*material.ambient.y;
        color.z += scene.ambient_light.z*material.ambient.z;
        
    }
    return color;
}


int main(int argc, char* argv[])
{
    scene.loadFromXml(argv[1]);
    BG_COLOR = scene.background_color;
    ShadowRayEpsilon = scene.shadow_ray_epsilon;
    MaxRecursionDepth = scene.max_recursion_depth;
    AmbientLight = scene.ambient_light;

    normal_calculator();

    for(int i=0; i<scene.cameras.size(); i++){
        Camera camera = scene.cameras[i];

        Vec3f eye = camera.position;
        Vec3f gaze = camera.gaze;
        Vec3f v_camera = camera.up;
        Vec3f u_camera = crossProd(v_camera, multS(gaze,-1.0F));
        float plane_left = camera.near_plane.x;
        float plane_rigth = camera.near_plane.y;
        float plane_bottom = camera.near_plane.z;
        float plane_top = camera.near_plane.w;
        float plane_distance = camera.near_distance;
        int width = camera.image_width;
        int height = camera.image_height;

        unsigned char*  image = new unsigned char [width * height * 3];

        Vec3f middle = add(eye, multS(gaze, plane_distance));
        Vec3f q_plane = add(middle, add(multS(v_camera,plane_top),multS(u_camera,plane_left))); 

        int k=0;
        for(int y = 0; y<height; y++){
            for(int x=0; x<width; x++){

                float s_u = ((float) x + 0.5) * (plane_rigth - plane_left) / ((float) width);
                float s_v = ((float) y + 0.5) * (plane_top - plane_bottom) / ((float) height);

                Vec3f s_plane = add(q_plane, add(multS(u_camera, s_u), multS(v_camera, -s_v)));
                Vec3f ray_d = add(multS(eye, -1), s_plane);
                Ray ray;
                ray.d = ray_d;
                ray.o = eye;
                Vec3i pixel_color = getPixelColor(ray,0);
                Vec3f pixel = add(eye, ray_d);
                image[k++] = (unsigned char) std::min(pixel_color.x, 255);
                image[k++] = (unsigned char) std::min(pixel_color.y, 255);
                image[k++] = (unsigned char) std::min(pixel_color.z, 255);
            }   
            //std::cout << height*width*3 - k << "\n";
            
            
            write_ppm(camera.image_name.c_str(), image, width, height);
        }

    }
}
