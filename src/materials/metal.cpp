#include "metal.h"

namespace pbrt
{

void MetalMaterial::ComputeScatteringFunctions(SurfaceInteraction* si,
                                               MemoryArena& arena,
                                               TransportMode mode,
                                               bool allowMultipleLobes) const 
{
	if (bumpMap)
	{
        Bump(bumpMap, si);
	}
}

}