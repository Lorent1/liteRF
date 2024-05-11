#include <iostream>
#include <fstream>
#include <vector>
#include <memory>  // for shared pointers
#include <iomanip> // for std::fixed/std::setprecision
#include <sstream>

#include "tracers/tracer.h"
#include "Image2d.h"

#ifdef USE_VULKAN
#include "vk_context.h"
std::shared_ptr<RayMarcherExample> CreateRayMarcherExample_Generated(vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated);
#endif


int main(int argc, const char** argv)
{
#ifndef NDEBUG
    bool enableValidationLayers = true;
#else
    bool enableValidationLayers = false;
#endif

    uint WIN_WIDTH = 128;
    uint WIN_HEIGHT = 128;

    std::shared_ptr<RayMarcherExample> pImpl = nullptr;
#ifdef USE_VULKAN
    bool onGPU = true; // TODO: you can read it from command line
    if (onGPU)
    {
        auto ctx = vk_utils::globalContextGet(enableValidationLayers, 0);
        pImpl = CreateRayMarcherExample_Generated(ctx, WIN_WIDTH * WIN_HEIGHT);
    }
    else
#else
    bool onGPU = false;
#endif
    pImpl = std::make_shared<RayMarcherExample>();

    pImpl->CommitDeviceData();

    std::vector<uint> pixelData(WIN_WIDTH * WIN_HEIGHT);
    char* filename = "./model.dat";

    size_t gridSize = 128;
    float factor = 1;
    std::ifstream fin(filename, std::ios::in | std::ios::binary);

    if (!fin) {
        std::cout << "Couldn't open file!";
        exit(1);
    }

    if (filename == "./model2.dat") {
        fin.read((char*)&gridSize, sizeof(int));
        fin.read((char*)&factor, sizeof(float));
    }

    pImpl->InitGrid(gridSize);
    pImpl->SetBoundingBox(float3(0, 0, 0), float3(1, 1, 1));
    pImpl->SetWorldViewMProjatrix(perspectiveMatrix(45, 1, 0.1, 100));

    fin.read((char*)&pImpl->grid[0], pImpl->grid.size() * sizeof(Cell));
    if (filename == "./model2.dat") {
        for (int i = 0; i < pImpl->grid.size(); i++) {
            pImpl->grid[i].density *= factor;
        }
    }
    fin.close();

    const int framesNum = 128;

    for (int k = 0; k < framesNum; k++) {
        float4x4 viewMat = lookAt(float3(0.0, 0.0, 1.3), float3(0.0, 0.0, 0.0), float3(0.0, 1.0, 0.0)) * rotate4x4Y(-float(360.0 / framesNum * k) * DEG_TO_RAD) * translate4x4(float3(-0.5, -0.5, -0.5));
        pImpl->SetWorldViewMatrix(viewMat);

        pImpl->UpdateMembersPlainData();                                            // copy all POD members from CPU to GPU in GPU implementation
        pImpl->RayMarch(pixelData.data(), WIN_WIDTH, WIN_HEIGHT);

        float timings[4] = { 0,0,0,0 };
        pImpl->GetExecutionTime("RayMarch", timings);

        std::stringstream strOut;
        if (onGPU)
            strOut << std::fixed << std::setprecision(2) << "out_gpu_" << k << ".bmp";
        else
            strOut << std::fixed << std::setprecision(2) << "out_cpu_" << k << ".bmp";
        std::string fileName = strOut.str();

        LiteImage::SaveBMP(fileName.c_str(), pixelData.data(), WIN_WIDTH, WIN_HEIGHT);

        std::cout << "img no. = " << k << ", timeRender = " << timings[0] << " ms, timeCopy = " << timings[1] + timings[2] << " ms " << std::endl;
    };

    pImpl = nullptr;
    return 0;
}