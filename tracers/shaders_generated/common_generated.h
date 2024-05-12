/////////////////////////////////////////////////////////////////////
/////////////  Required  Shader Features ////////////////////////////
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
/////////////////// include files ///////////////////////////////////
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
/////////////////// declarations in class ///////////////////////////
/////////////////////////////////////////////////////////////////////
#ifndef uint32_t
#define uint32_t uint
#endif
#define FLT_MAX 1e37f
#define FLT_MIN -1e37f
#define FLT_EPSILON 1e-6f
#define DEG_TO_RAD  0.017453293f
#define unmasked
#define half  float16_t
#define half2 f16vec2
#define half3 f16vec3
#define half4 f16vec4
#define SH_WIDTH 9
bool  isfinite(float x)            { return !isinf(x); }
float copysign(float mag, float s) { return abs(mag)*sign(s); }

struct complex
{
  float re, im;
};

complex make_complex(float re, float im) { 
  complex res;
  res.re = re;
  res.im = im;
  return res;
}

complex to_complex(float re)              { return make_complex(re, 0.0f);}
complex complex_add(complex a, complex b) { return make_complex(a.re + b.re, a.im + b.im); }
complex complex_sub(complex a, complex b) { return make_complex(a.re - b.re, a.im - b.im); }
complex complex_mul(complex a, complex b) { return make_complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
complex complex_div(complex a, complex b) {
  const float scale = 1 / (b.re * b.re + b.im * b.im);
  return make_complex(scale * (a.re * b.re + a.im * b.im), scale * (a.im * b.re - a.re * b.im));
}

complex real_add_complex(float value, complex z) { return complex_add(to_complex(value),z); }
complex real_sub_complex(float value, complex z) { return complex_sub(to_complex(value),z); }
complex real_mul_complex(float value, complex z) { return complex_mul(to_complex(value),z); }
complex real_div_complex(float value, complex z) { return complex_div(to_complex(value),z); }

complex complex_add_real(complex z, float value) { return complex_add(z, to_complex(value)); }
complex complex_sub_real(complex z, float value) { return complex_sub(z, to_complex(value)); }
complex complex_mul_real(complex z, float value) { return complex_mul(z, to_complex(value)); }
complex complex_div_real(complex z, float value) { return complex_div(z, to_complex(value)); }

float real(complex z) { return z.re;}
float imag(complex z) { return z.im; }
float complex_norm(complex z) { return z.re * z.re + z.im * z.im; }
float complex_abs(complex z) { return sqrt(complex_norm(z)); }
complex complex_sqrt(complex z) 
{
  float n = complex_abs(z);
  float t1 = sqrt(0.5f * (n + abs(z.re)));
  float t2 = 0.5f * z.im / t1;
  if (n == 0.0f)
    return to_complex(0.0f);
  if (z.re >= 0.0f)
    return make_complex(t1, t2);
  else
    return make_complex(abs(t2), copysign(t1, z.im));
}

struct Cell {
  float density;
  float sh_r[SH_WIDTH];
  float sh_g[SH_WIDTH];
  float sh_b[SH_WIDTH];
};
struct BoundingBox {
  vec3 min;
  vec3 max;
};

struct CellData {
    vec3 color;
    float density;
};

#ifndef SKIP_UBO_INCLUDE
#include "include/RayMarcherExample_generated_ubo.h"
#endif

/////////////////////////////////////////////////////////////////////
/////////////////// local functions /////////////////////////////////
/////////////////////////////////////////////////////////////////////

mat4 translate4x4(vec3 delta)
{
  return mat4(vec4(1.0, 0.0, 0.0, 0.0),
              vec4(0.0, 1.0, 0.0, 0.0),
              vec4(0.0, 0.0, 1.0, 0.0),
              vec4(delta, 1.0));
}

mat4 rotate4x4X(float phi)
{
  return mat4(vec4(1.0f, 0.0f,  0.0f,           0.0f),
              vec4(0.0f, +cos(phi),  +sin(phi), 0.0f),
              vec4(0.0f, -sin(phi),  +cos(phi), 0.0f),
              vec4(0.0f, 0.0f,       0.0f,      1.0f));
}

mat4 rotate4x4Y(float phi)
{
  return mat4(vec4(+cos(phi), 0.0f, -sin(phi), 0.0f),
              vec4(0.0f,      1.0f, 0.0f,      0.0f),
              vec4(+sin(phi), 0.0f, +cos(phi), 0.0f),
              vec4(0.0f,      0.0f, 0.0f,      1.0f));
}

mat4 rotate4x4Z(float phi)
{
  return mat4(vec4(+cos(phi), sin(phi), 0.0f, 0.0f),
              vec4(-sin(phi), cos(phi), 0.0f, 0.0f),
              vec4(0.0f,      0.0f,     1.0f, 0.0f),
              vec4(0.0f,      0.0f,     0.0f, 1.0f));
}

mat4 inverse4x4(mat4 m) { return inverse(m); }
vec3 mul4x3(mat4 m, vec3 v) { return (m*vec4(v, 1.0f)).xyz; }
vec3 mul3x3(mat4 m, vec3 v) { return (m*vec4(v, 0.0f)).xyz; }

mat3 make_float3x3(vec3 a, vec3 b, vec3 c) { // different way than mat3(a,b,c)
  return mat3(a.x, b.x, c.x,
              a.y, b.y, c.y,
              a.z, b.z, c.z);
}

void sh_eval_2(in vec3 d, inout float sh_coeffs[SH_WIDTH]) {
    float x = d.x, y = d.y, z = d.z, z2 = z * z;
    float c0, c1, s0, s1, tmp_a, tmp_b, tmp_c;

    sh_coeffs[0] = 0.28209479177387814;
    sh_coeffs[2] = z * 0.488602511902919923;
    sh_coeffs[6] = z2 * 0.94617469575756008 + -0.315391565252520045;
    c0 = x;
    s0 = y;

    tmp_a = -0.488602511902919978;
    sh_coeffs[3] = tmp_a * c0;
    sh_coeffs[1] = tmp_a * s0;
    tmp_b = z * -1.09254843059207896;
    sh_coeffs[7] = tmp_b * c0;
    sh_coeffs[5] = tmp_b * s0;
    c1 = x * c0 - y * s0;
    s1 = x * s0 + y * c0;

    tmp_c = 0.546274215296039478;
    sh_coeffs[8] = tmp_c * c1;
    sh_coeffs[4] = tmp_c * s1;
}

Cell mult(Cell cell, float number) {
  cell.density *= number;

  for (uint i = 0; i < SH_WIDTH; i++) {
      cell.sh_r[i] *= number;
      cell.sh_g[i] *= number;
      cell.sh_b[i] *= number;
  }

  return cell;
}

Cell add(Cell rh, Cell lh) {
  rh.density += lh.density;

  for (uint i = 0; i < SH_WIDTH; i++) {
      rh.sh_r[i] += lh.sh_r[i];
      rh.sh_g[i] += lh.sh_g[i];
      rh.sh_b[i] += lh.sh_b[i];
  }

  return rh;
}

int indexOf(vec3 pos, int gridSize) {
    return int(pos.x + gridSize * pos.y + gridSize * gridSize * pos.z);
}

float eval_sh(inout float sh[SH_WIDTH], vec3 rayDir) {
  float sh_coeffs[SH_WIDTH];
  sh_eval_2(rayDir, sh_coeffs);

  float sum = 0.0f;
  for (int i = 0; i < SH_WIDTH; i++)
    sum += sh[i] * sh_coeffs[i];

  return sum;
}

Cell trilerp(inout Cell values[8], vec3 pos) {
    vec3 n = vec3(1.0f - pos.x,1.0f - pos.y,1.0f - pos.z);

    Cell first = mult(add(mult(values[0], n.z), mult(values[1], pos.z)), n.y);
    Cell second = mult(add(mult(values[2], n.z), mult(values[3], pos.z)), pos.y);
    Cell third = mult(add(mult(values[4], n.z), mult(values[5], pos.z)), n.y);
    Cell fourth = mult(add(mult(values[6], n.z), mult(values[7], pos.z)), pos.y);

    return add(mult(add(first, second), n.x), mult(add(third, fourth), pos.x));
}

CellData eval_trilinear(vec3 pos, vec3 rd, inout Cell gridData[10], int gridSize) {
    const float EPS = 0.01f;

    CellData result;
    result.color = vec3(0.0);
    result.density = 0;

    if (pos.x < EPS || pos.y < EPS || pos.z < EPS
        || pos.x > 1.0 - EPS || pos.y > 1.0 - EPS || pos.z > 1.0 - EPS) {
        return result;
    }

    pos *= float(gridSize);

    vec3 floor_pos = floor(vec3(pos.x,pos.y,pos.z));

    if (floor_pos.x + 1 >= float(gridSize) || floor_pos.y + 1 >= float(gridSize) || floor_pos.z + 1 >= float(gridSize)) {
        return result;
    }

    vec3 indices[8];
    Cell values[8];

    indices[0] = vec3(floor_pos.x,floor_pos.y,floor_pos.z);
    indices[1] = vec3(floor_pos.x,floor_pos.y,floor_pos.z + 1);
    indices[2] = vec3(floor_pos.x,floor_pos.y + 1,floor_pos.z);
    indices[3] = vec3(floor_pos.x,floor_pos.y + 1,floor_pos.z + 1);
    indices[4] = vec3(floor_pos.x + 1,floor_pos.y,floor_pos.z);
    indices[5] = vec3(floor_pos.x + 1,floor_pos.y,floor_pos.z + 1);
    indices[6] = vec3(floor_pos.x + 1,floor_pos.y + 1,floor_pos.z);
    indices[7] = vec3(floor_pos.x + 1,floor_pos.y + 1,floor_pos.z + 1);

    for (int i = 0; i < 8; i++) {
        //std::cout << "\n" << indices[i].z << "\n";
        values[i] = gridData[indexOf(indices[i], gridSize)];
    }

    vec3 coeffs = vec3(fract(pos.x),fract(pos.y),fract(pos.z));
    Cell cell = trilerp(values, coeffs);

    result.density = cell.density;
    result.color = vec3(eval_sh(cell.sh_r, rd),eval_sh(cell.sh_g, rd),eval_sh(cell.sh_b, rd));

    return result;
}

uint RealColorToUint32(vec4 real_color) {
    int red = max(0, min(255, int((real_color[0] * 255.0f))));
    int green = max(0, min(255, int((real_color[1] * 255.0f))));
    int blue = max(0, min(255, int((real_color[2] * 255.0f))));
    int alpha = max(0, min(255, int((real_color[3] * 255.0f))));

    return red | (green << 8) | (blue << 16) | (alpha << 24);
}

void transform_ray3f(mat4 a_mWorldViewInv, inout vec3 ray_pos, inout vec3 ray_dir) {
  vec4 rayPosTransformed = a_mWorldViewInv*vec4(ray_pos, 1.0f);
  vec4 rayDirTransformed = a_mWorldViewInv*vec4(ray_dir, 0.0f);
  
  (ray_pos) = rayPosTransformed.xyz;
  (ray_dir) = (normalize(rayDirTransformed)).xyz;
}

vec2 RayBoxIntersection(vec3 ray_pos, vec3 ray_dir, vec3 boxMin, vec3 boxMax) {
  ray_dir.x = 1.0f/ray_dir.x; // may precompute if intersect many boxes
  ray_dir.y = 1.0f/ray_dir.y; // may precompute if intersect many boxes
  ray_dir.z = 1.0f/ray_dir.z; // may precompute if intersect many boxes

  float lo = ray_dir.x*(boxMin.x - ray_pos.x);
  float hi = ray_dir.x*(boxMax.x - ray_pos.x);
  
  float tmin = min(lo, hi);
  float tmax = max(lo, hi);

  float lo1 = ray_dir.y*(boxMin.y - ray_pos.y);
  float hi1 = ray_dir.y*(boxMax.y - ray_pos.y);

  tmin = max(tmin, min(lo1, hi1));
  tmax = min(tmax, max(lo1, hi1));

  float lo2 = ray_dir.z*(boxMin.z - ray_pos.z);
  float hi2 = ray_dir.z*(boxMax.z - ray_pos.z);

  tmin = max(tmin, min(lo2, hi2));
  tmax = min(tmax, max(lo2, hi2));
  
  return vec2(tmin,tmax);
}

vec3 EyeRayDir(float x, float y, mat4 a_mViewProjInv) {
  vec4 pos = vec4(2.0f*x - 1.0f,2.0f*y - 1.0f,0.0f,1.0f);
  pos = a_mViewProjInv * pos;
  pos /= pos.w;
  return normalize(pos.xyz);
}

#define KGEN_FLAG_RETURN            1
#define KGEN_FLAG_BREAK             2
#define KGEN_FLAG_DONT_SET_EXIT     4
#define KGEN_FLAG_SET_EXIT_NEGATIVE 8
#define KGEN_REDUCTION_LAST_STEP    16
#define CFLOAT_GUARDIAN 
#define MAXFLOAT FLT_MAX

