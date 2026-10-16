#pragma once
#include "../detail/type_vec4.hpp"

namespace glm
{

	typedef vec<4, double, defaultp>		dvec4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::dvec4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::dvec4>::value);
static_assert(std::is_trivially_copyable<glm::dvec4>::value);
static_assert(std::is_copy_constructible<glm::dvec4>::value);
static_assert(glm::dvec4::length() == 4);
