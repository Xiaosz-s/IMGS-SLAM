#pragma once
#include "../detail/type_mat2x3.hpp"

namespace glm
{

	typedef mat<2, 3, float, defaultp>		mat2x3;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::mat2x3>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::mat2x3>::value);
static_assert(std::is_trivially_copyable<glm::mat2x3>::value);
static_assert(std::is_copy_constructible<glm::mat2x3>::value);
static_assert(glm::mat2x3::length() == 2);

