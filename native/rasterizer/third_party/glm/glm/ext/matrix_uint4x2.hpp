#pragma once

#include "../mat4x2.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_uint4x2 extension included")
#endif

namespace glm
{

	typedef mat<4, 2, uint, defaultp>	umat4x2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::umat4x2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::umat4x2>::value);
static_assert(std::is_trivially_copyable<glm::umat4x2>::value);
static_assert(std::is_copy_constructible<glm::umat4x2>::value);
static_assert(glm::umat4x2::length() == 4);

