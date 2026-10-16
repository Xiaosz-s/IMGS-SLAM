#pragma once
#include "../detail/type_mat2x2.hpp"

namespace glm
{

	typedef mat<2, 2, float, defaultp>		mat2x2;

	typedef mat<2, 2, float, defaultp>		mat2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::mat2x2>::value);
static_assert(std::is_trivially_default_constructible<glm::mat2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::mat2x2>::value);
static_assert(std::is_trivially_copy_assignable<glm::mat2>::value);
static_assert(std::is_trivially_copyable<glm::mat2x2>::value);
static_assert(std::is_trivially_copyable<glm::mat2>::value);
static_assert(std::is_copy_constructible<glm::mat2x2>::value);
static_assert(std::is_copy_constructible<glm::mat2>::value);
static_assert(glm::mat2x2::length() == 2);
static_assert(glm::mat2::length() == 2);
