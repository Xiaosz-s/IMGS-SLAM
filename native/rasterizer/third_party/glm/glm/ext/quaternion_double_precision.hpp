#pragma once

#include "../detail/type_quat.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_quaternion_double_precision extension included")
#endif

namespace glm
{

	typedef qua<double, lowp>		lowp_dquat;

	typedef qua<double, mediump>	mediump_dquat;

	typedef qua<double, highp>		highp_dquat;

} 

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::lowp_dquat>::value);
static_assert(std::is_trivially_default_constructible<glm::mediump_dquat>::value);
static_assert(std::is_trivially_default_constructible<glm::highp_dquat>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::lowp_dquat>::value);
static_assert(std::is_trivially_copy_assignable<glm::mediump_dquat>::value);
static_assert(std::is_trivially_copy_assignable<glm::highp_dquat>::value);
static_assert(std::is_trivially_copyable<glm::lowp_dquat>::value);
static_assert(std::is_trivially_copyable<glm::mediump_dquat>::value);
static_assert(std::is_trivially_copyable<glm::highp_dquat>::value);
static_assert(std::is_copy_constructible<glm::lowp_dquat>::value);
static_assert(std::is_copy_constructible<glm::mediump_dquat>::value);
static_assert(std::is_copy_constructible<glm::highp_dquat>::value);
static_assert(glm::lowp_dquat::length() == 4);
static_assert(glm::mediump_dquat::length() == 4);
static_assert(glm::highp_dquat::length() == 4);
