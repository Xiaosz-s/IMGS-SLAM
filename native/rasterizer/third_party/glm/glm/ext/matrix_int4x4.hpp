#pragma once

#include "../mat4x4.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_int4x4 extension included")
#endif

namespace glm
{

	typedef mat<4, 4, int, defaultp>	imat4x4;

	typedef mat<4, 4, int, defaultp>	imat4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::imat4x4>::value);
static_assert(std::is_trivially_default_constructible<glm::imat4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::imat4x4>::value);
static_assert(std::is_trivially_copy_assignable<glm::imat4>::value);
static_assert(std::is_trivially_copyable<glm::imat4x4>::value);
static_assert(std::is_trivially_copyable<glm::imat4>::value);
static_assert(std::is_copy_constructible<glm::imat4x4>::value);
static_assert(std::is_copy_constructible<glm::imat4>::value);
static_assert(glm::imat4x4::length() == 4);
static_assert(glm::imat4::length() == 4);
