#ifndef RayMarcherExample_UBO_H
#define RayMarcherExample_UBO_H

#ifndef GLSL
#define LAYOUT_STD140
#include "LiteMath.h"
using LiteMath::uint;
typedef LiteMath::float4x4 mat4;
typedef LiteMath::float2   vec2;
typedef LiteMath::float3   vec3;
typedef LiteMath::float4   vec4;
typedef LiteMath::int2     ivec2;
typedef LiteMath::int3     ivec3;
typedef LiteMath::int4     ivec4;
typedef LiteMath::uint2    uvec2;
typedef LiteMath::uint3    uvec3;
typedef LiteMath::uint4    uvec4;
#else
#define M_PI          3.14159265358979323846f
#define M_TWOPI       6.28318530717958647692f
#define INV_PI        0.31830988618379067154f
#define INV_TWOPI     0.15915494309189533577f
#endif

struct RayMarcherExample_Generated_UBO_Data
{
  mat4 m_worldViewInv; 
  mat4 m_worldViewProjInv; 
  BoundingBox bb; 
  #ifndef GLSL 
  uint dummy0; 
  #endif 
  #ifndef GLSL 
  uint dummy1; 
  #endif 
  int gridSize; 
  float rayMarchTime; 
  uint grid_capacity; 
  uint grid_size; 
  uint dummy_last;
};

#endif

