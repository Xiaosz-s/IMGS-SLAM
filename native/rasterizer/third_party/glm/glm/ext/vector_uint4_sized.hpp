#pragma once

#include "../ext/vector_uint4.hpp"
#include "../ext/scalar_uint_sized.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_vector_uint4_sized extension included")
#endif

namespace glm
{

	typedef vec<4, uint8, defaultp>		u8vec4;

	typedef vec<4, uint16, defaultp>	u16vec4;

	typedef vec<4, uint32, defaultp>	u32vec4;

	typedef vec<4, uint64, defaultp>	u64vec4;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::u8vec4>::value);
static_assert(std::is_trivially_default_constructible<glm::u16vec4>::value);
static_assert(std::is_trivially_default_constructible<glm::u32vec4>::value);
static_assert(std::is_trivially_default_constructible<glm::u64vec4>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::u8vec4>::value);
static_assert(std::is_trivially_copy_assignable<glm::u16vec4>::value);
static_assert(std::is_trivially_copy_assignable<glm::u32vec4>::value);
static_assert(std::is_trivially_copy_assignable<glm::u64vec4>::value);
static_assert(std::is_trivially_copyable<glm::u8vec4>::value);
static_assert(std::is_trivially_copyable<glm::u16vec4>::value);
static_assert(std::is_trivially_copyable<glm::u32vec4>::value);
static_assert(std::is_trivially_copyable<glm::u64vec4>::value);
static_assert(std::is_copy_constructible<glm::u8vec4>::value);
static_assert(std::is_copy_constructible<glm::u16vec4>::value);
static_assert(std::is_copy_constructible<glm::u32vec4>::value);
static_assert(std::is_copy_constructible<glm::u64vec4>::value);
static_assert(glm::u8vec4::length() == 4);
static_assert(glm::u16vec4::length() == 4);
static_assert(glm::u32vec4::length() == 4);
static_assert(glm::u64vec4::length() == 4);
