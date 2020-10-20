#pragma once

#include "pbrt.h"
#include "material.h"

namespace pbrt {

class GlassMaterial : public Material {
  public:
    GlassMaterial(const std::shared_ptr<Texture<Spectrum>> &Kr,
                  const std::shared_ptr<Texture<Spectrum>> &Kt,
                  const std::shared_ptr<Texture<Float>> &uRoughness,
                  const std::shared_ptr<Texture<Float>> &vRoughness,
                  const std::shared_ptr<Texture<Float>> &index,
                  const std::shared_ptr<Texture<Float>> &bumpMap,
                  bool remapRoughness)
        : Kr(Kr),
          Kt(Kt),
          uRoughness(uRoughness),
          vRoughness(vRoughness),
          index(index),
          bumpMap(bumpMap),
          remapRoughness(remapRoughness) {}

    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;


  private:
    std::shared_ptr<Texture<Spectrum>> Kr, Kt;
    std::shared_ptr<Texture<Float>> uRoughness, vRoughness;
    std::shared_ptr<Texture<Float>> index;
    std::shared_ptr<Texture<Float>> bumpMap;
    bool remapRoughness;
};

}