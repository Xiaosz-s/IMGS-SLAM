#pragma once

#include "../detail/type_quat.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_quaternion_double extension included")
#endif

namespace glm
{

	typedef qua<double, defaultp>		dquat;

} 

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::dquat>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::dquat>::value);
static_assert(std::is_trivially_copyable<glm::dquat>::value);
static_assert(std::is_copy_constructible<glm::dquat>::value);
static_assert(glm::dquat::length() == 4);
