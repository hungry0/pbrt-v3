#include "mirror.h"
#include "memory.h"
#include "reflection.h"
#include "interaction.h"
#include "texture.h"
#include "paramset.h"

namespace pbrt
{
    
void MirrorMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                                MemoryArena &arena,
                                                TransportMode mode,
                                                bool allowMultipleLobes) const 
{
    if(bumpMap)
        Bump(bumpMap, si);

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
    Spectrum R = Kr->Evaluate(*si).Clamp();
    if(!R.IsBlack())
        si->bsdf->Add(ARENA_ALLOC(arena, SpecularReflection)(
            R, ARENA_ALLOC(arena, FresnelNoOp)()));
}

MirrorMaterial *CreateMirrorMaterial(const TextureParams& mp) 
{
    std::shared_ptr<Texture<Spectrum>> Kr =
        mp.GetSpectrumTexture("Kr", Spectrum(1.f));

    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");

    return new MirrorMaterial(Kr, bumpMap);
}

}