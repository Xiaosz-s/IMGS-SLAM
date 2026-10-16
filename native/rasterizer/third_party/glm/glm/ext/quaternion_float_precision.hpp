#pragma once

#include "../detail/type_quat.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_quaternion_float_precision extension included")
#endif

namespace glm
{

	typedef qua<float, lowp>		lowp_quat;

	typedef qua<float, mediump>		mediump_quat;

	typedef qua<float, highp>		highp_quat;

} 

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::lowp_quat>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_quat>::value);
static_assert(std::is_trivially_default_constructible<glm::highp_quat>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::lowp_quat>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_quat>::value);
static_assert(std::is_trivially_copy_assignable<glm::highp_quat>::value);
static_assert(std::is_trivially_copyable<glm::lowp_quat>::value);
static_assert(std::is_trivially_copyable<glm::mediump_quat>::value);
static_assert(std::is_trivially_copyable<glm::highp_quat>::value);
static_assert(std::is_copy_constructible<glm::lowp_quat>::value);
static_assert(std::is_copy_constructible<glm::mediump_quat>::value);
static_assert(std::is_copy_constructible<glm::highp_quat>::value);
static_assert(glm::lowp_quat::length() == 4);
static_assert(glm::mediump_quat::length() == 4);
static_assert(glm::highp_quat::length() == 4);
