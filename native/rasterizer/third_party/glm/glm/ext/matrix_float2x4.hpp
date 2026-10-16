#pragma once
#include "../detail/type_mat2x4.hpp"

namespace glm
{

	typedef mat<2, 4, float, defaultp>		mat2x4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::mat2x4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::mat2x4>::value);
static_assert(std::is_trivially_copyable<glm::mat2x4>::value);
static_assert(std::is_copy_constructible<glm::mat2x4>::value);
static_assert(glm::mat2x4::length() == 2);

