#pragma once

#include "../ext/vector_uint2.hpp"
#include "../ext/scalar_uint_sized.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_vector_uint2_sized extension included")
#endif

namespace glm
{

	typedef vec<2, uint8, defaultp>		u8vec2;

	typedef vec<2, uint16, defaultp>	u16vec2;

	typedef vec<2, uint32, defaultp>	u32vec2;

	typedef vec<2, uint64, defaultp>	u64vec2;

}

#if GLM_CONFIG_CTOR_INIT == GLM_DISABLE
static_assert(std::is_trivially_default_constructible<glm::u8vec2>::value);
static_assert(std::is_trivially_default_constructible<glm::u16vec2>::value);
static_assert(std::is_trivially_default_constructible<glm::u32vec2>::value);
static_assert(std::is_trivially_default_constructible<glm::u64vec2>::value);
#endif
static_assert(std::is_trivially_copy_assignable<glm::u8vec2>::value);
static_assert(std::is_trivially_copy_assignable<glm::u16vec2>::value);
static_assert(std::is_trivially_copy_assignable<glm::u32vec2>::value);
static_assert(std::is_trivially_copy_assignable<glm::u64vec2>::value);
static_assert(std::is_trivially_copyable<glm::u8vec2>::value);
static_assert(std::is_trivially_copyable<glm::u16vec2>::value);
static_assert(std::is_trivially_copyable<glm::u32vec2>::value);
static_assert(std::is_trivially_copyable<glm::u64vec2>::value);
static_assert(std::is_copy_constructible<glm::u8vec2>::value);
static_assert(std::is_copy_constructible<glm::u16vec2>::value);
static_assert(std::is_copy_constructible<glm::u32vec2>::value);
static_assert(std::is_copy_constructible<glm::u64vec2>::value);
static_assert(glm::u8vec2::length() == 2);
static_assert(glm::u16vec2::length() == 2);
static_assert(glm::u32vec2::length() == 2);
static_assert(glm::u64vec2::length() == 2);
