#pragma once
#include "../detail/type_mat4x2.hpp"

namespace glm
{

	typedef mat<4, 2, double, defaultp>		dmat4x2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::dmat4x2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::dmat4x2>::value);
static_assert(std::is_trivially_copyable<glm::dmat4x2>::value);
static_assert(std::is_copy_constructible<glm::dmat4x2>::value);
static_assert(glm::dmat4x2::length() == 4);
