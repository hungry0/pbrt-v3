#include "translucent.h"
#include "material.h"
#include "memory.h"
#include "interaction.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"

namespace pbrt 
{


void TranslucentMaterial::ComputeScatteringFunctions(
    SurfaceInteraction* si, MemoryArena& arena, TransportMode mode,
    bool allowMultipleLobes) const 
{
    if (bumpMap)
    {
        Bump(bumpMap, si);
    }

    Float eta = 1.5f;
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, eta);

    Spectrum r = reflect->Evaluate(*si).Clamp();
    Spectrum t = transmit->Evaluate(*si).Clamp();
    if (r.IsBlack() && t.IsBlack()) return;

    Spectrum kd = Kd->Evaluate(*si).Clamp();
    if (!kd.IsBlack())
    {
        if (!r.IsBlack()) 
        {
            si->bsdf->Add(ARENA_ALLOC(arena, LambertianReflection)(r * kd));
        }

        if (!t.IsBlack())
        {
            si->bsdf->Add(ARENA_ALLOC(arena, LambertianTransmission)(t * kd));
        }
    }

    Spectrum ks = Ks->Evaluate(*si).Clamp();
    if (!ks.IsBlack() && (!r.IsBlack() || !t.IsBlack()))
    {
        Float rough = roughness->Evaluate(*si);
        if (remapRoughness) 
        {
            rough = TrowbridgeReitzDistribution::RoughnessToAlpha(rough);
        }

        MicrofacetDistribution* distrib =
            ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(rough, rough);

        if (!r.IsBlack())
        {
            Fresnel* fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, eta);
            si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetReflection)(
                r * ks, distrib, fresnel));
        }

        if (!t.IsBlack())
        {
            si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetTransmission)(
                t * ks, distrib, 1.f, eta, mode));
        }
    }

}

TranslucentMaterial* CreateTranslucentMaterial(const TextureParams& mp) 
{
    std::shared_ptr<Texture<Spectrum>> Kd =
        mp.GetSpectrumTexture("Kd", Spectrum(0.5f));
    std::shared_ptr<Texture<Spectrum>> Ks =
        mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    std::shared_ptr<Texture<Spectrum>> reflect =
        mp.GetSpectrumTexture("reflect", Spectrum(0.5f));
    std::shared_ptr<Texture<Spectrum>> transmit =
        mp.GetSpectrumTexture("transmit", Spectrum(0.5f));
    std::shared_ptr<Texture<Float>> roughness =
        mp.GetFloatTexture("roughness", 0.1f);
    std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");

    bool remapRoughness = mp.FindBool("remaproughness", true);
    return new TranslucentMaterial(Kd, Ks, roughness, reflect, transmit,
                                   bumpMap, remapRoughness);
}

}  // namespace pbrt