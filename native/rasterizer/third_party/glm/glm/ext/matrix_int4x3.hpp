#pragma once

#include "../mat4x3.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_int4x3 extension included")
#endif

namespace glm
{

	typedef mat<4, 3, int, defaultp>	imat4x3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::imat4x3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::imat4x3>::value);
static_assert(std::is_trivially_copyable<glm::imat4x3>::value);
static_assert(std::is_copy_constructible<glm::imat4x3>::value);
static_assert(glm::imat4x3::length() == 4);
