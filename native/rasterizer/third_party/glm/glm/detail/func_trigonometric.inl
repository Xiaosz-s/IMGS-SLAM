#include "_vectorize.hpp"
#include <cmath>
#include <limits>

namespace glm
{
	
	template<typename genType>
	GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType radians(genType degrees)
	{
		static_assert(std::numeric_limits<genType>::is_iec559 || GLM_CONFIG_UNRESTRICTED_FLOAT, "'radians' only accept floating-point input");

		return degrees * static_cast<genType>(0.01745329251994329576923690768489);
	}

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<L, T, Q> radians(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(radians, v);
	}

	template<typename genType>
	GLM_FUNC_QUALIFIER GLM_CONSTEXPR genType degrees(genType radians)
	{
		static_assert(std::numeric_limits<genType>::is_iec559 || GLM_CONFIG_UNRESTRICTED_FLOAT, "'degrees' only accept floating-point input");

		return radians * static_cast<genType>(57.295779513082320876798154814105);
	}

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER GLM_CONSTEXPR vec<L, T, Q> degrees(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(degrees, v);
	}

	using ::std::sin;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> sin(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(sin, v);
	}

	using std::cos;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> cos(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(cos, v);
	}

	using std::tan;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> tan(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(tan, v);
	}

	using std::asin;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> asin(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(asin, v);
	}

	using std::acos;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> acos(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(acos, v);
	}

	template<typename genType>
	GLM_FUNC_QUALIFIER genType atan(genType y, genType x)
	{
		static_assert(std::numeric_limits<genType>::is_iec559 || GLM_CONFIG_UNRESTRICTED_FLOAT, "'atan' only accept floating-point input");

		return ::std::atan2(y, x);
	}

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> atan(vec<L, T, Q> const& y, vec<L, T, Q> const& x)
	{
		return detail::functor2<vec, L, T, Q>::call(::std::atan2, y, x);
	}

	using std::atan;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> atan(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(atan, v);
	}

	using std::sinh;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> sinh(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(sinh, v);
	}

	using std::cosh;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> cosh(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(cosh, v);
	}

	using std::tanh;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> tanh(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(tanh, v);
	}

	using std::asinh;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> asinh(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(asinh, v);
	}

	using std::acosh;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> acosh(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(acosh, v);
	}

	using std::atanh;

	template<length_t L, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER vec<L, T, Q> atanh(vec<L, T, Q> const& v)
	{
		return detail::functor1<vec, L, T, T, Q>::call(atanh, v);
	}
}

#if GLM_CONFIG_SIMD == GLM_ENABLE
#	include "func_trigonometric_simd.inl"
#endif

