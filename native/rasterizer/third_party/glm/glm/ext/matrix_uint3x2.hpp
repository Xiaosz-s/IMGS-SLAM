#pragma once

#include "../mat3x2.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_uint3x2 extension included")
#endif

namespace glm
{

	typedef mat<3, 2, uint, defaultp>	umat3x2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::umat3x2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::umat3x2>::value);
static_assert(std::is_trivially_copyable<glm::umat3x2>::value);
static_assert(std::is_copy_constructible<glm::umat3x2>::value);
static_assert(glm::umat3x2::length() == 3);

