#pragma once
#include "../detail/type_mat3x2.hpp"

namespace glm
{

	typedef mat<3, 2, double, lowp>		lowp_dmat3x2;

	typedef mat<3, 2, double, mediump>	mediump_dmat3x2;

	typedef mat<3, 2, double, highp>	highp_dmat3x2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::lowp_dmat3x2>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_dmat3x2>::value);
static_assert(std::is_trivially_default_constructible<glm::highp_dmat3x2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::lowp_dmat3x2>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_dmat3x2>::value);
static_assert(std::is_trivially_copy_assignable<glm::highp_dmat3x2>::value);
static_assert(std::is_trivially_copyable<glm::lowp_dmat3x2>::value);
static_assert(std::is_trivially_copyable<glm::mediump_dmat3x2>::value);
static_assert(std::is_trivially_copyable<glm::highp_dmat3x2>::value);
static_assert(std::is_copy_constructible<glm::lowp_dmat3x2>::value);
static_assert(std::is_copy_constructible<glm::mediump_dmat3x2>::value);
static_assert(std::is_copy_constructible<glm::highp_dmat3x2>::value);
static_assert(glm::lowp_dmat3x2::length() == 3);
static_assert(glm::mediump_dmat3x2::length() == 3);
static_assert(glm::highp_dmat3x2::length() == 3);
