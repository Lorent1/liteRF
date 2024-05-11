#include <vector>
#include <chrono>
#include <string>

#include "tracer.h"

// From Mitsuba 3
void sh_eval_2(const float3 &d, float* sh_coeffs) {
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

float eval_sh(float* sh, float3 rayDir) {
  float sh_coeffs[SH_WIDTH];
  sh_eval_2(rayDir, sh_coeffs);

  float sum = 0.0f;
  for (int i = 0; i < SH_WIDTH; i++)
    sum += sh[i] * sh_coeffs[i];

  return sum;
}

// Cell trilerp(Cell* values, float3 pos ) {
//     float3 n = { 1.0f - pos.x, 1.0f - pos.y, 1.0f - pos.z };

//     return ((values[0b000] * n.z + values[0b001] * pos.z) * n.y
//         + (values[0b010] * n.z + values[0b011] * pos.z) * pos.y) * n.x
//         + ((values[0b100] * n.z + values[0b101] * pos.z) * n.y
//             + (values[0b110] * n.z + values[0b111] * pos.z) * pos.y) * pos.x;
// }

Cell mult(Cell cell, float number){
  cell.density *= number;

  for (uint i = 0; i < SH_WIDTH; i++) {
      cell.sh_r[i] *= number;
      cell.sh_g[i] *= number;
      cell.sh_b[i] *= number;
  }

  return cell;
}

Cell add(Cell rh, Cell lh){
  rh.density += lh.density;

  for (uint i = 0; i < SH_WIDTH; i++) {
      rh.sh_r[i] += lh.sh_r[i];
      rh.sh_g[i] += lh.sh_g[i];
      rh.sh_b[i] += lh.sh_b[i];
  }

  return rh;
}

Cell trilerp(Cell* values, float3 pos) {
    float3 n = float3(1.0f - pos.x,1.0f - pos.y,1.0f - pos.z);

    Cell first = mult(add(mult(values[0], n.z), mult(values[1], pos.z)), n.y);
    Cell second = mult(add(mult(values[2], n.z), mult(values[3], pos.z)), pos.y);
    Cell third = mult(add(mult(values[4], n.z), mult(values[5], pos.z)), n.y);
    Cell fourth = mult(add(mult(values[6], n.z), mult(values[7], pos.z)), pos.y);


    // return ((values[0b000] * n.z + values[0b001] * pos.z) * n.y
    //     + (values[0b010] * n.z + values[0b011] * pos.z) * pos.y) * n.x
    //     + ((values[0b100] * n.z + values[0b101] * pos.z) * n.y
    //         + (values[0b110] * n.z + values[0b111] * pos.z) * pos.y) * pos.x;
    return add(mult(add(first, second), n.x), mult(add(third, fourth), pos.x));
}

int indexOf(float3 pos, int gridSize) {
    return pos.x + gridSize * pos.y + gridSize * gridSize * pos.z;
}


CellData eval_trilinear(float3 pos, float3 rd, Cell* gridData, int gridSize) {
    const float EPS = 0.01f;

    CellData result;
    result.color = float3(0.0);
    result.density = 0;

    if (pos.x < EPS || pos.y < EPS || pos.z < EPS
        || pos.x > 1.0 - EPS || pos.y > 1.0 - EPS || pos.z > 1.0 - EPS) {
        return result;
    }

    pos *= gridSize;

    float3 floor_pos = floor({ pos.x, pos.y, pos.z });

    if (floor_pos.x + 1 > gridSize || floor_pos.y + 1 > gridSize || floor_pos.z + 1 > gridSize) {
        return result;
    }

    float3 indices[8];
    Cell values[8];

    indices[0] = { floor_pos.x, floor_pos.y, floor_pos.z };
    indices[1] = { floor_pos.x, floor_pos.y, floor_pos.z + 1 };
    indices[2] = { floor_pos.x, floor_pos.y + 1, floor_pos.z };
    indices[3] = { floor_pos.x, floor_pos.y + 1, floor_pos.z + 1 };
    indices[4] = { floor_pos.x + 1, floor_pos.y, floor_pos.z };
    indices[5] = { floor_pos.x + 1, floor_pos.y, floor_pos.z + 1 };
    indices[6] = { floor_pos.x + 1, floor_pos.y + 1, floor_pos.z };
    indices[7] = { floor_pos.x + 1, floor_pos.y + 1, floor_pos.z + 1 };

    for (int i = 0; i < 8; i++) {
        values[i] = gridData[indexOf(indices[i], gridSize)];
    }

    float3 coeffs = { fract(pos.x), fract(pos.y), fract(pos.z) };
    Cell cell = trilerp(values, coeffs);

    result.density = cell.density;
    result.color = { eval_sh(cell.sh_r, rd), eval_sh(cell.sh_g, rd), eval_sh(cell.sh_b, rd) };

    return result;
}

float2 RayBoxIntersection(float3 ray_pos, float3 ray_dir, float3 boxMin, float3 boxMax){
  ray_dir.x = 1.0f/ray_dir.x; // may precompute if intersect many boxes
  ray_dir.y = 1.0f/ray_dir.y; // may precompute if intersect many boxes
  ray_dir.z = 1.0f/ray_dir.z; // may precompute if intersect many boxes

  float lo = ray_dir.x*(boxMin.x - ray_pos.x);
  float hi = ray_dir.x*(boxMax.x - ray_pos.x);
  
  float tmin = std::min(lo, hi);
  float tmax = std::max(lo, hi);

  float lo1 = ray_dir.y*(boxMin.y - ray_pos.y);
  float hi1 = ray_dir.y*(boxMax.y - ray_pos.y);

  tmin = std::max(tmin, std::min(lo1, hi1));
  tmax = std::min(tmax, std::max(lo1, hi1));

  float lo2 = ray_dir.z*(boxMin.z - ray_pos.z);
  float hi2 = ray_dir.z*(boxMax.z - ray_pos.z);

  tmin = std::max(tmin, std::min(lo2, hi2));
  tmax = std::min(tmax, std::max(lo2, hi2));
  
  return float2(tmin, tmax);
}

static inline float3 EyeRayDir(float x, float y, float4x4 a_mViewProjInv)
{
  float4 pos = float4(2.0f*x - 1.0f, 2.0f*y - 1.0f, 0.0f, 1.0f );
  pos = a_mViewProjInv * pos;
  pos /= pos.w;
  return normalize(to_float3(pos));
}

static inline void transform_ray3f(float4x4 a_mWorldViewInv, float3* ray_pos, float3* ray_dir) 
{
  float4 rayPosTransformed = a_mWorldViewInv*to_float4(*ray_pos, 1.0f);
  float4 rayDirTransformed = a_mWorldViewInv*to_float4(*ray_dir, 0.0f);
  
  (*ray_pos) = to_float3(rayPosTransformed);
  (*ray_dir) = to_float3(normalize(rayDirTransformed));
}

float4 RaymarchSpherical(float3 ro, float3 rd, float tmin, float tmax, float& alpha, Cell* gridData, int gridSize) {
    float stepSize =(tmax - tmin) / 250.0f;

    float3 color = float3(0.0f);
    float densitySum = 0.0f;

    for (int i = 0; i < 250; i++) {
        float3 pos = ro + rd * lerp(tmin, tmax, (float)i / 249.0);
        CellData cell = eval_trilinear(pos, rd, gridData, gridSize);

        cell.density = max(0.0f, cell.density);
        cell.color = clamp(cell.color * 0.32f + 0.54f, float3(0.0), float3(1.0));

        color += cell.color
            * exp(-densitySum)
            * (1.0 - exp(-cell.density * stepSize));

        densitySum += stepSize * cell.density;
    }

    //std::cout << color.x << " " << color.y << " " << color.z << "\n";

    return to_float4(color, alpha);
}

static inline uint32_t RealColorToUint32(float4 real_color) {
    int red = std::max(0, std::min(255, (int)(real_color[0] * 255.0f)));
    int green = std::max(0, std::min(255, (int)(real_color[1] * 255.0f)));
    int blue = std::max(0, std::min(255, (int)(real_color[2] * 255.0f)));
    int alpha = std::max(0, std::min(255, (int)(real_color[3] * 255.0f)));

    return red | (green << 8) | (blue << 16) | (alpha << 24);
}

void RayMarcherExample::kernel2D_RayMarch(uint32_t* out_color, uint32_t width, uint32_t height) {
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float3 rayDir = EyeRayDir((float(x) + 0.5f) / float(width), (float(y) + 0.5f) / float(height), m_worldViewProjInv);
            float3 rayPos = float3(0.0f, 0.26f, 0.0f);

            transform_ray3f(m_worldViewInv, &rayPos, &rayDir);

            float2 tNearAndFar = RayBoxIntersection(rayPos, rayDir, bb.min, bb.max);

            float4 resColor(0.0f);
            float alpha = 1.0f;
            resColor = RaymarchSpherical(rayPos, rayDir, tNearAndFar.x, tNearAndFar.y, alpha, grid.data(), gridSize);

            //std::cout << resColor.x << "\n";
            out_color[y * width + x] = RealColorToUint32(resColor);
        };
    };
}

void RayMarcherExample::RayMarch(uint32_t* out_color, uint32_t width, uint32_t height)
{ 
  auto start = std::chrono::high_resolution_clock::now();
  kernel2D_RayMarch(out_color, width, height);
  rayMarchTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}  

void RayMarcherExample::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  if(std::string(a_funcName) == "RayMarch")
    a_out[0] =  rayMarchTime;
}
