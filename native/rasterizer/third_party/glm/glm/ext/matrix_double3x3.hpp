#pragma once
#include "../detail/type_mat3x3.hpp"

namespace glm
{

	typedef mat<3, 3, double, defaultp>		dmat3x3;

	typedef mat<3, 3, double, defaultp>		dmat3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::dmat3x3>::value);
static_assert(std::is_trivially_default_constructible<glm::dmat3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::dmat3x3>::value);
static_assert(std::is_trivially_copy_assignable<glm::dmat3>::value);
static_assert(std::is_trivially_copyable<glm::dmat3x3>::value);
static_assert(std::is_trivially_copyable<glm::dmat3>::value);
static_assert(std::is_copy_constructible<glm::dmat3x3>::value);
static_assert(std::is_copy_constructible<glm::dmat3>::value);
static_assert(glm::dmat3x3::length() == 3);
static_assert(glm::dmat3::length() == 3);

