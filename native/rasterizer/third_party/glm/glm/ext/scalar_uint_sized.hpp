#pragma once

#include "../detail/setup.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_scalar_uint_sized extension included")
#endif

namespace glm{

	typedef detail::uint8		uint8;

	typedef detail::uint16		uint16;

	typedef detail::uint32		uint32;

	typedef detail::uint64		uint64;

}
