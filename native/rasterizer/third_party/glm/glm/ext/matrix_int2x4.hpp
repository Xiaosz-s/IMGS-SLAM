#pragma once

#include "../mat2x4.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_int2x4 extension included")
#endif

namespace glm
{

	typedef mat<2, 4, int, defaultp>	imat2x4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::imat2x4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::imat2x4>::value);
static_assert(std::is_trivially_copyable<glm::imat2x4>::value);
static_assert(std::is_copy_constructible<glm::imat2x4>::value);
static_assert(glm::imat2x4::length() == 2);

