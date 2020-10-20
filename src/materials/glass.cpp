#include "materials/glass.h"
#include "material.h"
#include "texture.h"
#include "memory.h"
#include "interaction.h"
#include "spectrum.h"
#include "reflection.h"

namespace pbrt {

void GlassMaterial::ComputeScatteringFunctions(SurfaceInteraction* si,
                                               MemoryArena& arena,
                                               TransportMode mode,
                                               bool allowMultipleLobes) const
{
    if (bumpMap) Bump(bumpMap, si);
    Float eta = index->Evaluate(*si);
    Float urough = uRoughness->Evaluate(*si);
    Float vrough = vRoughness->Evaluate(*si);
    Spectrum R = Kr->Evaluate(*si).Clamp();
    Spectrum T = Kt->Evaluate(*si).Clamp();

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, eta);

    if (R.IsBlack() && T.IsBlack()) return;

    bool isSpecular = urough == 0 && vrough == 0;
    if (isSpecular && allowMultipleLobes)
    {
        si->bsdf->Add(
            ARENA_ALLOC(arena, FresnelSpecular)(R, T, 1.f, eta, mode));
    } 
    else 
    {
        if (remapRoughness)
        {
            urough = TrowbridgeReitzDistribution::RoughnessToAlpha(urough);
            vrough = TrowbridgeReitzDistribution::RoughnessToAlpha(vrough);
        }

        MicrofacetDistribution* distrib =
            isSpecular ? nullptr
                : ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(urough, vrough);

        if (!R.IsBlack())
        {
            Fresnel* fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, eta);
            if (isSpecular) 
            {  
                // si->bsdf->Add(
                //    ARENA_ALLOC(arena, SpecularReflection)(R, fresnel));
            } 
            else 
            {
                si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetReflection)(R, distrib, fresnel));
            }
        }

        //if (!T.IsBlack())
        //{
        //    if (isSpecular)
        //        si->bsdf->Add(ARENA_ALLOC(arena, SpecularTransmission)(
        //            T, 1.f, eta, mode));
        //    else
        //        si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetTransmission)(
        //            T, distrib, 1.f, eta, mode));
        //}
    }
}

}


