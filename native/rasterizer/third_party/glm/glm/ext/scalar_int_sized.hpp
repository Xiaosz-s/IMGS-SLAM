#pragma once

#include "../detail/setup.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_scalar_int_sized extension included")
#endif

namespace glm{

	typedef detail::int8		int8;

	typedef detail::int16		int16;

	typedef detail::int32		int32;

	typedef detail::int64		int64;

}
