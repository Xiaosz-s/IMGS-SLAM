#pragma once
#include "../detail/type_mat2x2.hpp"

namespace glm
{

	typedef mat<2, 2, double, defaultp>		dmat2x2;

	typedef mat<2, 2, double, defaultp>		dmat2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::dmat2x2>::value);
static_assert(std::is_trivially_default_constructible<glm::dmat2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::dmat2x2>::value);
static_assert(std::is_trivially_copy_assignable<glm::dmat2>::value);
static_assert(std::is_trivially_copyable<glm::dmat2x2>::value);
static_assert(std::is_trivially_copyable<glm::dmat2>::value);
static_assert(std::is_copy_constructible<glm::dmat2x2>::value);
static_assert(std::is_copy_constructible<glm::dmat2>::value);
static_assert(glm::dmat2x2::length() == 2);
static_assert(glm::dmat2::length() == 2);

