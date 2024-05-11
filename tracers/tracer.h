#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdint>

#include "LiteMath.h"
using namespace LiteMath;

const size_t SH_WIDTH = 9;

struct Cell {
  float density;
  float sh_r[SH_WIDTH];
  float sh_g[SH_WIDTH];
  float sh_b[SH_WIDTH];
};

// static inline Cell operator* (Cell cell, float number) {
//     cell.density *= number;

//     for (size_t i = 0; i < SH_WIDTH; i++) {
//         cell.sh_r[i] *= number;
//         cell.sh_g[i] *= number;
//         cell.sh_b[i] *= number;
//     }

//     return cell;
// };

// static inline Cell operator+ (Cell a, Cell b) {
//     a.density += b.density;

//     for (size_t i = 0; i < SH_WIDTH; i++) {
//         a.sh_r[i] += b.sh_r[i];
//         a.sh_g[i] += b.sh_g[i];
//         a.sh_b[i] += b.sh_b[i];
//     }
    
//     return a;
// };

struct BoundingBox {
  float3 min;
  float3 max;
};

struct CellData {
    float3 color;
    float density;
};

class RayMarcherExample{ // : public IRenderAPI
public:

  RayMarcherExample()
  {
    const float4x4 view = lookAt(float3(0,1.5,-3), float3(0,0,0), float3(0,1,0)); // pos, look_at, up
    const float4x4 proj = perspectiveMatrix(90.0f, 1.0f, 0.1f, 100.0f);
    m_worldViewInv      = inverse4x4(view); 
    m_worldViewProjInv  = inverse4x4(proj); 
  }

  // 28 * 4 gridSize = 128

  void InitGrid(const float _gridSize) {
    gridSize = _gridSize;
    grid.resize(gridSize * gridSize * gridSize);

    for (size_t i = 0; i < gridSize * gridSize * gridSize; i++) {
      grid[i].density = 0.02;
      for (size_t j = 0; j < SH_WIDTH; j++) {
        grid[i].sh_r[j] = 0.1;
        grid[i].sh_g[j] = 0.1;
        grid[i].sh_b[j] = 0.1;
      }
    }
  }

  void SetBoundingBox(const float3 boxMin, const float3 boxMax) {
    bb.min = boxMin;
    bb.max = boxMax;
  }

  void SetWorldViewMProjatrix(const float4x4& a_mat) {m_worldViewProjInv = inverse4x4(a_mat);}
  void SetWorldViewMatrix(const float4x4& a_mat) {m_worldViewInv = inverse4x4(a_mat);}

  virtual void kernel2D_RayMarch(uint32_t* out_color, uint32_t width, uint32_t height);
  virtual void RayMarch(uint32_t* out_color [[size("width*height")]], uint32_t width, uint32_t height);  

  virtual void CommitDeviceData() {}                                       // will be overriden in generated class
  virtual void UpdateMembersPlainData() {}                                 // will be overriden in generated class (optional function)
  //virtual void UpdateMembersVectorData() {}                              // will be overriden in generated class (optional function)
  //virtual void UpdateMembersTexureData() {}                              // will be overriden in generated class (optional function)
  virtual void GetExecutionTime(const char* a_funcName, float a_out[4]);   // will be overriden in generated class
  
  std::vector<Cell> grid;
protected:
  int gridSize;
  BoundingBox bb;

  float4x4 m_worldViewProjInv;
  float4x4 m_worldViewInv;
  float    rayMarchTime;
};
