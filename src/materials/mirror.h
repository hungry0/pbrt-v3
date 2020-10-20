#pragma once
#include "material.h"

namespace pbrt 
{
	class MirrorMaterial : public Material
	{
	public:
		MirrorMaterial(const std::shared_ptr<Texture<Spectrum>>& r,
						const std::shared_ptr<Texture<Float>>& bump)
			:Kr(r),bumpMap(bump){}
	
		void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena, TransportMode mode, bool allowMultipleLobes) const;

	private:
		std::shared_ptr<Texture<Spectrum>> Kr;
		std::shared_ptr<Texture<Float>> bumpMap;
	};

	MirrorMaterial* CreateMirrorMaterial(const TextureParams& mp);
}
