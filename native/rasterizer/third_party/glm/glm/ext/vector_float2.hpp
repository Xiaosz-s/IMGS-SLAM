#pragma once
#include "../detail/type_vec2.hpp"

namespace glm
{

	typedef vec<2, float, defaultp>	vec2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::vec2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::vec2>::value);
static_assert(std::is_trivially_copyable<glm::vec2>::value);
static_assert(std::is_copy_constructible<glm::vec2>::value);
static_assert(glm::vec2::length() == 2);
