#pragma once

#include "../gtc/bitfield.hpp"

#ifndef GLM_ENABLE_EXPERIMENTAL
#	error "GLM: GLM_GTX_bit is an experimental extension and may change in the future. Use #define GLM_ENABLE_EXPERIMENTAL before including it, if you really want to use it."
#elif GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_GTX_bit extension included")
#endif

namespace glm
{

	template<typename genIUType>
	GLM_FUNC_DECL genIUType highestBitValue(genIUType Value);

	template<typename genIUType>
	GLM_FUNC_DECL genIUType lowestBitValue(genIUType Value);

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_DECL vec<L, T, Q> highestBitValue(vec<L, T, Q> const& value);

	template<typename genIUType>
	[[deprecated("Use ceilPowerOfTwo from GTC_round instead")]]
	GLM_FUNC_DECL genIUType powerOfTwoAbove(genIUType Value);

	template<length_t L, typename T, qualifier Q>
	[[deprecated("Use ceilPowerOfTwo from GTC_round instead")]]
	GLM_FUNC_DECL vec<L, T, Q> powerOfTwoAbove(vec<L, T, Q> const& value);

	template<typename genIUType>
	[[deprecated("Use floorPowerOfTwo from GTC_round instead")]]
	GLM_FUNC_DECL genIUType powerOfTwoBelow(genIUType Value);

	template<length_t L, typename T, qualifier Q>
	[[deprecated("Use floorPowerOfTwo from GTC_round instead")]]
	GLM_FUNC_DECL vec<L, T, Q> powerOfTwoBelow(vec<L, T, Q> const& value);

	template<typename genIUType>
	[[deprecated("Use roundPowerOfTwo from GTC_round instead")]]
	GLM_FUNC_DECL genIUType powerOfTwoNearest(genIUType Value);

	template<length_t L, typename T, qualifier Q>
	[[deprecated("Use roundPowerOfTwo from GTC_round instead")]]
	GLM_FUNC_DECL vec<L, T, Q> powerOfTwoNearest(vec<L, T, Q> const& value);

} 

#include "bit.inl"

