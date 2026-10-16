#pragma once

#include "../mat2x3.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_uint2x3 extension included")
#endif

namespace glm
{

	typedef mat<2, 3, uint, defaultp>	umat2x3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::umat2x3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::umat2x3>::value);
static_assert(std::is_trivially_copyable<glm::umat2x3>::value);
static_assert(std::is_copy_constructible<glm::umat2x3>::value);
static_assert(glm::umat2x3::length() == 2);

