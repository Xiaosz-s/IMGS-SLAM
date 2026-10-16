#pragma once
#include "../detail/type_vec4.hpp"

namespace glm
{

	typedef vec<4, float, defaultp>		vec4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::vec4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::vec4>::value);
static_assert(std::is_trivially_copyable<glm::vec4>::value);
static_assert(std::is_copy_constructible<glm::vec4>::value);
static_assert(glm::vec4::length() == 4);

