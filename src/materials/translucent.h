#pragma once
#include "material.h"

namespace pbrt 
{

class TranslucentMaterial : public Material {
  public:
    TranslucentMaterial(std::shared_ptr<Texture<Spectrum>>& kd,
                        std::shared_ptr<Texture<Spectrum>>& ks,
                        std::shared_ptr<Texture<Float>>& rough,
                        std::shared_ptr<Texture<Spectrum>>& refl,
                        std::shared_ptr<Texture<Spectrum>>& trans,
                        std::shared_ptr<Texture<Float>>& bump, bool remap)
        : Kd(kd),
          Ks(ks),
          roughness(rough),
          reflect(refl),
          transmit(trans),
          bumpMap(bump),
          remapRoughness(remap) {}

    void ComputeScatteringFunctions(SurfaceInteraction* si, MemoryArena& arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    std::shared_ptr<Texture<Spectrum>> Kd, Ks;
    std::shared_ptr<Texture<Float>> roughness;
    std::shared_ptr<Texture<Spectrum>> reflect, transmit;
    std::shared_ptr<Texture<Float>> bumpMap;
    bool remapRoughness;
};

TranslucentMaterial* CreateTranslucentMaterial(const TextureParams& mp);

}
