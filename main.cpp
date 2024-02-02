#include <iostream>
#include <fstream>
#include <vector>
#include <memory>  // for shared pointers
#include <iomanip> // for std::fixed/std::setprecision

#include "example_tracer/example_tracer.h"
#include "Image2d.h"

//#include "vk_context.h"
//std::shared_ptr<TestClass> CreateTestClass_Generated(vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated);

int main(int argc, const char** argv)
{
  #ifndef NDEBUG
  bool enableValidationLayers = true;
  #else
  bool enableValidationLayers = false;
  #endif

  uint WIN_WIDTH  = 512;
  uint WIN_HEIGHT = 512;

  bool onGPU = false;
  std::shared_ptr<RayMarcherExample> pImpl = nullptr;
  //if(onGPU)
  //{
  //  unsigned int a_preferredDeviceId = args.getOptionValue<int>("--gpu_id", 0);
  //  auto ctx = vk_utils::globalContextGet(enableValidationLayers, a_preferredDeviceId);
  //  pImpl = CreateTestClass_Generated(ctx, WIN_WIDTH*WIN_HEIGHT);
  //}
  //else
    pImpl = std::make_shared<RayMarcherExample>();

  pImpl->CommitDeviceData();

  std::vector<uint> pixelData(WIN_WIDTH*WIN_HEIGHT);  

  for(int angleY = 0; angleY < 360; angleY += 10) 
  {
    float4x4 mRot    = rotate4x4Y(float(angleY)*DEG_TO_RAD);
    float4   camPos  = mRot*float4(0,0,-3,0) + float4(0,1.5f,0,1);              // rotate and than translate camera position
    float4x4 viewMat = lookAt(to_float3(camPos), float3(0,0,0), float3(0,1,0)); // pos, look_at, up  
    
    pImpl->SetWorldViewMatrix(viewMat);
    pImpl->UpdateMembersPlainData();                                            // copy all POD members from CPU to GPU in GPU implementation
    pImpl->RayMarch(pixelData.data(), WIN_WIDTH, WIN_HEIGHT);
  
    float timings[4] = {0,0,0,0};
    pImpl->GetExecutionTime("RayMarch", timings);

    std::stringstream strOut;
    if(onGPU)
      strOut << std::fixed << std::setprecision(2) << "out_gpu_" << angleY << ".bmp";
    else
      strOut << std::fixed << std::setprecision(2) << "out_cpu_" << angleY << ".bmp";
    std::string fileName = strOut.str();

    LiteImage::SaveBMP(fileName.c_str(), pixelData.data(), WIN_WIDTH, WIN_HEIGHT);

    std::cout << "angl = " << angleY << ", timeRender = " << timings[0] << " ms, timeCopy = " <<  timings[1] + timings[2] << " ms " << std::endl;
  }
  
  pImpl = nullptr;
  return 0;
}