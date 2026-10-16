#pragma once

#include "../mat3x3.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_uint3x3 extension included")
#endif

namespace glm
{

	typedef mat<3, 3, uint, defaultp>	umat3x3;

	typedef mat<3, 3, uint, defaultp>	umat3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::umat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::umat3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::umat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::umat3>::value);
static_assert(std::is_trivially_copyable<glm::umat3x3>::value);
static_assert(std::is_trivially_copyable<glm::umat3>::value);
static_assert(std::is_copy_constructible<glm::umat3x3>::value);
static_assert(std::is_copy_constructible<glm::umat3>::value);
static_assert(glm::umat3x3::length() == 3);
static_assert(glm::umat3::length() == 3);
