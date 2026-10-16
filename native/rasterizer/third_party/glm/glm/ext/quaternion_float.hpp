#pragma once

#include "../detail/type_quat.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_quaternion_float extension included")
#endif

namespace glm
{

	typedef qua<float, defaultp>		quat;

} 

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::quat>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::quat>::value);
static_assert(std::is_trivially_copyable<glm::quat>::value);
static_assert(std::is_copy_constructible<glm::quat>::value);
static_assert(glm::quat::length() == 4);
