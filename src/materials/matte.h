#pragma once
#include "material.h"

namespace pbrt {

// Sigma��OrenNayarģ����Ҫ�õ�һ������
class MatteMaterial : public Material {
  public:
    MatteMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
                  const std::shared_ptr<Texture<Float>> &sigma,
                  const std::shared_ptr<Texture<Float>> &bumpMap) {}

	void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    std::shared_ptr<Texture<Spectrum>> Kd;
    std::shared_ptr<Texture<Float>> sigma, bumpMap;

};

MatteMaterial *CreateMatteMaterial(const TextureParams &mp);

}  // namespace pbrt
