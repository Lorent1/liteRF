#include <vector>
#include <memory>
#include <limits>
#include <cassert>
#include <chrono>
#include <array>

#include "vk_copy.h"
#include "vk_context.h"
#include "vk_images.h"

#include "tracer_generated.h"
#include "include/RayMarcherExample_generated_ubo.h"

static uint32_t ComputeReductionSteps(uint32_t whole_size, uint32_t wg_size)
{
  uint32_t steps = 0;
  while (whole_size > 1)
  {
    steps++;
    whole_size = (whole_size + wg_size - 1) / wg_size;
  }
  return steps;
}

constexpr uint32_t KGEN_REDUCTION_LAST_STEP    = 16;


void RayMarcherExample_Generated::UpdatePlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  const size_t maxAllowedSize = std::numeric_limits<uint32_t>::max();
  m_uboData.m_worldViewInv = m_worldViewInv;
  m_uboData.m_worldViewProjInv = m_worldViewProjInv;
  m_uboData.bb = bb;
  m_uboData.gridSize = gridSize;
  m_uboData.rayMarchTime = rayMarchTime;
  a_pCopyEngine->UpdateBuffer(m_classDataBuffer, 0, &m_uboData, sizeof(m_uboData));
}

void RayMarcherExample_Generated::ReadPlainMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
  a_pCopyEngine->ReadBuffer(m_classDataBuffer, 0, &m_uboData, sizeof(m_uboData));
  m_worldViewInv = m_uboData.m_worldViewInv;
  m_worldViewProjInv = m_uboData.m_worldViewProjInv;
  bb = m_uboData.bb;
  gridSize = m_uboData.gridSize;
  rayMarchTime = m_uboData.rayMarchTime;
}

void RayMarcherExample_Generated::UpdateVectorMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{
}

void RayMarcherExample_Generated::UpdateTextureMembers(std::shared_ptr<vk_utils::ICopyEngine> a_pCopyEngine)
{ 
}

void RayMarcherExample_Generated::RayMarchCmd(uint32_t* out_color, Cell* grid, uint32_t width, uint32_t height)
{
  uint32_t blockSizeX = 32;
  uint32_t blockSizeY = 8;
  uint32_t blockSizeZ = 1;
  
  struct KernelArgsPC
  {
    uint32_t m_sizeX;
    uint32_t m_sizeY;
    uint32_t m_sizeZ;
    uint32_t m_tFlags;
  } pcData;
  
  uint32_t sizeX  = uint32_t(height);
  uint32_t sizeY  = uint32_t(width);
  uint32_t sizeZ  = uint32_t(1);
  
  pcData.m_sizeX  = height;
  pcData.m_sizeY  = width;
  pcData.m_sizeZ  = 1;
  pcData.m_tFlags = m_currThreadFlags;
  vkCmdPushConstants(m_currCmdBuffer, RayMarchLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(KernelArgsPC), &pcData);
  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, RayMarchPipeline);
  vkCmdDispatch    (m_currCmdBuffer, (sizeX + blockSizeX - 1) / blockSizeX, (sizeY + blockSizeY - 1) / blockSizeY, (sizeZ + blockSizeZ - 1) / blockSizeZ);
}


void RayMarcherExample_Generated::copyKernelFloatCmd(uint32_t length)
{
  uint32_t blockSizeX = MEMCPY_BLOCK_SIZE;

  vkCmdBindPipeline(m_currCmdBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, copyKernelFloatPipeline);
  vkCmdPushConstants(m_currCmdBuffer, copyKernelFloatLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &length);
  vkCmdDispatch(m_currCmdBuffer, (length + blockSizeX - 1) / blockSizeX, 1, 1);
}

VkBufferMemoryBarrier RayMarcherExample_Generated::BarrierForClearFlags(VkBuffer a_buffer)
{
  VkBufferMemoryBarrier bar = {};
  bar.sType               = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
  bar.pNext               = NULL;
  bar.srcAccessMask       = VK_ACCESS_TRANSFER_WRITE_BIT;
  bar.dstAccessMask       = VK_ACCESS_SHADER_READ_BIT;
  bar.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.buffer              = a_buffer;
  bar.offset              = 0;
  bar.size                = VK_WHOLE_SIZE;
  return bar;
}

VkBufferMemoryBarrier RayMarcherExample_Generated::BarrierForSingleBuffer(VkBuffer a_buffer)
{
  VkBufferMemoryBarrier bar = {};
  bar.sType               = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
  bar.pNext               = NULL;
  bar.srcAccessMask       = VK_ACCESS_SHADER_WRITE_BIT;
  bar.dstAccessMask       = VK_ACCESS_SHADER_READ_BIT;
  bar.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
  bar.buffer              = a_buffer;
  bar.offset              = 0;
  bar.size                = VK_WHOLE_SIZE;
  return bar;
}

void RayMarcherExample_Generated::BarriersForSeveralBuffers(VkBuffer* a_inBuffers, VkBufferMemoryBarrier* a_outBarriers, uint32_t a_buffersNum)
{
  for(uint32_t i=0; i<a_buffersNum;i++)
  {
    a_outBarriers[i].sType               = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
    a_outBarriers[i].pNext               = NULL;
    a_outBarriers[i].srcAccessMask       = VK_ACCESS_SHADER_WRITE_BIT;
    a_outBarriers[i].dstAccessMask       = VK_ACCESS_SHADER_READ_BIT;
    a_outBarriers[i].srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    a_outBarriers[i].dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    a_outBarriers[i].buffer              = a_inBuffers[i];
    a_outBarriers[i].offset              = 0;
    a_outBarriers[i].size                = VK_WHOLE_SIZE;
  }
}

void RayMarcherExample_Generated::RayMarchCmd(VkCommandBuffer a_commandBuffer, uint32_t* out_color, Cell* grid, uint32_t width, uint32_t height)
{
  m_currCmdBuffer = a_commandBuffer;
  VkMemoryBarrier memoryBarrier = { VK_STRUCTURE_TYPE_MEMORY_BARRIER, nullptr, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT }; 
  
  auto start = std::chrono::high_resolution_clock::now();
  vkCmdBindDescriptorSets(a_commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, RayMarchLayout, 0, 1, &m_allGeneratedDS[0], 0, nullptr);
  RayMarchCmd(out_color, grid, width, height);
  vkCmdPipelineBarrier(m_currCmdBuffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memoryBarrier, 0, nullptr, 0, nullptr);
  rayMarchTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;

}



void RayMarcherExample_Generated::RayMarch(uint32_t* out_color, Cell* grid, uint32_t width, uint32_t height)
{
  // (1) get global Vulkan context objects
  //
  VkInstance       instance       = m_ctx.instance;
  VkPhysicalDevice physicalDevice = m_ctx.physicalDevice;
  VkDevice         device         = m_ctx.device;
  VkCommandPool    commandPool    = m_ctx.commandPool; 
  VkQueue          computeQueue   = m_ctx.computeQueue; 
  VkQueue          transferQueue  = m_ctx.transferQueue;
  auto             pCopyHelper    = m_ctx.pCopyHelper;
  auto             pAllocatorSpec = m_ctx.pAllocatorSpecial;

  // (2) create GPU objects
  //
  auto outFlags = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
  if(RayMarch_local.needToClearOutput)
    outFlags |= VK_IMAGE_USAGE_TRANSFER_DST_BIT;
  std::vector<VkBuffer> buffers;
  std::vector<VkImage>  images2;
  std::vector<vk_utils::VulkanImageMem*> images;
  auto beforeCreateObjects = std::chrono::high_resolution_clock::now();
  VkBuffer out_colorGPU = vk_utils::createBuffer(device, width*height*sizeof(uint32_t ), outFlags);
  buffers.push_back(out_colorGPU);
  VkBuffer gridGPU = vk_utils::createBuffer(device, gridSize*gridSize*gridSize*sizeof(Cell ), outFlags);
  buffers.push_back(gridGPU);
  

  VkDeviceMemory buffersMem = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, buffers);
  VkDeviceMemory imagesMem  = VK_NULL_HANDLE; // vk_utils::allocateAndBindWithPadding(device, physicalDevice, std::vector<VkBuffer>(), images);
  
  vk_utils::MemAllocInfo tempMemoryAllocInfo;
  tempMemoryAllocInfo.memUsage = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; // TODO, select depending on device and sample/application (???)
  if(buffers.size() != 0)
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, buffers);
  if(images.size() != 0)
  {
    pAllocatorSpec->Allocate(tempMemoryAllocInfo, images2);
    for(auto imgMem : images)
    {
      VkImageViewCreateInfo imageView{};
      imageView.sType    = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
      imageView.viewType = VK_IMAGE_VIEW_TYPE_2D;
      imageView.image    = imgMem->image;
      imageView.format   = imgMem->format;
      imageView.subresourceRange = {};
      imageView.subresourceRange.aspectMask     = imgMem->aspectMask;
      imageView.subresourceRange.baseMipLevel   = 0;
      imageView.subresourceRange.levelCount     = imgMem->mipLvls;
      imageView.subresourceRange.baseArrayLayer = 0;
      imageView.subresourceRange.layerCount     = 1;
      VK_CHECK_RESULT(vkCreateImageView(device, &imageView, nullptr, &imgMem->view));
    }
  }
  
  auto afterCreateObjects = std::chrono::high_resolution_clock::now();
  m_exTimeRayMarch.msAPIOverhead = std::chrono::duration_cast<std::chrono::microseconds>(afterCreateObjects - beforeCreateObjects).count()/1000.f;
  
  auto afterCopy2 = std::chrono::high_resolution_clock::now(); // just declare it here, replace value later
  
  auto afterInitBuffers = std::chrono::high_resolution_clock::now();
  m_exTimeRayMarch.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(afterInitBuffers - afterCreateObjects).count()/1000.f;
  
  auto beforeSetInOut = std::chrono::high_resolution_clock::now();
  this->SetVulkanInOutFor_RayMarch(out_colorGPU, 0, gridGPU, 0, 0); 

  // (3) copy input data to GPU
  //
  auto beforeCopy = std::chrono::high_resolution_clock::now();
  m_exTimeRayMarch.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(beforeCopy - beforeSetInOut).count()/1000.f;
  auto afterCopy = std::chrono::high_resolution_clock::now();
  m_exTimeRayMarch.msCopyToGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy - beforeCopy).count()/1000.f;
  //  
  m_exTimeRayMarch.msExecuteOnGPU = 0;

  // (4) now execute algorithm on GPU
  //
  {
    VkCommandBuffer commandBuffer = vk_utils::createCommandBuffer(device, commandPool);
    VkCommandBufferBeginInfo beginCommandBufferInfo = {};
    beginCommandBufferInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginCommandBufferInfo.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginCommandBufferInfo);
    RayMarchCmd(commandBuffer, out_color, grid, width, height);      
    vkEndCommandBuffer(commandBuffer);  
    auto start = std::chrono::high_resolution_clock::now();
    vk_utils::executeCommandBufferNow(commandBuffer, computeQueue, device);
    vkFreeCommandBuffers(device, commandPool, 1, &commandBuffer);
    auto stop = std::chrono::high_resolution_clock::now();
    m_exTimeRayMarch.msExecuteOnGPU += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()/1000.f;
  }
  
  // (5) copy output data to CPU
  //
  auto beforeCopy2 = std::chrono::high_resolution_clock::now();
  pCopyHelper->ReadBuffer(out_colorGPU, 0, out_color, width*height*sizeof(uint32_t ));
  pCopyHelper->ReadBuffer(gridGPU, 0, grid, gridSize*gridSize*gridSize*sizeof(Cell ));
  this->ReadPlainMembers(pCopyHelper);
  afterCopy2 = std::chrono::high_resolution_clock::now();
  m_exTimeRayMarch.msCopyFromGPU = std::chrono::duration_cast<std::chrono::microseconds>(afterCopy2 - beforeCopy2).count()/1000.f;

  // (6) free resources 
  //
  vkDestroyBuffer(device, out_colorGPU, nullptr);
  vkDestroyBuffer(device, gridGPU, nullptr);
  if(buffersMem != VK_NULL_HANDLE)
    vkFreeMemory(device, buffersMem, nullptr);
  if(imagesMem != VK_NULL_HANDLE)
    vkFreeMemory(device, imagesMem, nullptr);
  
  m_exTimeRayMarch.msAPIOverhead += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - afterCopy2).count()/1000.f;
}

void RayMarcherExample_Generated::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  vk_utils::ExecTime res = {};
  if(std::string(a_funcName) == "RayMarch" || std::string(a_funcName) == "RayMarchBlock")
    res = m_exTimeRayMarch;
  a_out[0] = res.msExecuteOnGPU;
  a_out[1] = res.msCopyToGPU;
  a_out[2] = res.msCopyFromGPU;
  a_out[3] = res.msAPIOverhead;             
}

