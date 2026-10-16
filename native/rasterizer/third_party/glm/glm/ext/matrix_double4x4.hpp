#pragma once
#include "../detail/type_mat4x4.hpp"

namespace glm
{

	typedef mat<4, 4, double, defaultp>		dmat4x4;

	typedef mat<4, 4, double, defaultp>		dmat4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::dmat4x4>::value);
static_assert(std::is_trivially_default_constructible<glm::dmat4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::dmat4x4>::value);
static_assert(std::is_trivially_copy_assignable<glm::dmat4>::value);
static_assert(std::is_trivially_copyable<glm::dmat4x4>::value);
static_assert(std::is_trivially_copyable<glm::dmat4>::value);
static_assert(std::is_copy_constructible<glm::dmat4x4>::value);
static_assert(std::is_copy_constructible<glm::dmat4>::value);
static_assert(glm::dmat4x4::length() == 4);
static_assert(glm::dmat4::length() == 4);

