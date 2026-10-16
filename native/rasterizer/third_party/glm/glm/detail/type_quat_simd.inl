#if GLM_ARCH & GLM_ARCH_SSE2_BIT

namespace glm{
namespace detail
{

	template<qualifier Q>
	struct compute_quat_add<float, Q, true>
	{
		static qua<float, Q> call(qua<float, Q> const& q, qua<float, Q> const& p)
		{
			qua<float, Q> Result;
			Result.data = _mm_add_ps(q.data, p.data);
			return Result;
		}
	};

#	if GLM_ARCH & GLM_ARCH_AVX_BIT
	template<qualifier Q>
	struct compute_quat_add<double, Q, true>
	{
		static qua<double, Q> call(qua<double, Q> const& a, qua<double, Q> const& b)
		{
			qua<double, Q> Result;
			Result.data = _mm256_add_pd(a.data, b.data);
			return Result;
		}
	};
#	endif

	template<qualifier Q>
	struct compute_quat_sub<float, Q, true>
	{
		static qua<float, Q> call(qua<float, Q> const& q, qua<float, Q> const& p)
		{
			qua<float, Q> Result;
			Result.data = _mm_sub_ps(q.data, p.data);
			return Result;
		}
	};

#	if GLM_ARCH & GLM_ARCH_AVX_BIT
	template<qualifier Q>
	struct compute_quat_sub<double, Q, true>
	{
		static qua<double, Q> call(qua<double, Q> const& a, qua<double, Q> const& b)
		{
			qua<double, Q> Result;
			Result.data = _mm256_sub_pd(a.data, b.data);
			return Result;
		}
	};
#	endif

	template<qualifier Q>
	struct compute_quat_mul_scalar<float, Q, true>
	{
		static qua<float, Q> call(qua<float, Q> const& q, float s)
		{
			qua<float, Q> Result;
			Result.data = _mm_mul_ps(q.data, _mm_set_ps1(s));
			return Result;
		}
	};

#	if GLM_ARCH & GLM_ARCH_AVX_BIT
	template<qualifier Q>
	struct compute_quat_mul_scalar<double, Q, true>
	{
		static qua<double, Q> call(qua<double, Q> const& q, double s)
		{
			qua<double, Q> Result;
			Result.data = _mm256_mul_pd(q.data, _mm_set_ps1(s));
			return Result;
		}
	};
#	endif

	template<qualifier Q>
	struct compute_quat_div_scalar<float, Q, true>
	{
		static qua<float, Q> call(qua<float, Q> const& q, float s)
		{
			qua<float, Q> Result;
			Result.data = _mm_div_ps(q.data, _mm_set_ps1(s));
			return Result;
		}
	};

#	if GLM_ARCH & GLM_ARCH_AVX_BIT
	template<qualifier Q>
	struct compute_quat_div_scalar<double, Q, true>
	{
		static qua<double, Q> call(qua<double, Q> const& q, double s)
		{
			qua<double, Q> Result;
			Result.data = _mm256_div_pd(q.data, _mm_set_ps1(s));
			return Result;
		}
	};
#	endif

	template<qualifier Q>
	struct compute_quat_mul_vec4<float, Q, true>
	{
		static vec<4, float, Q> call(qua<float, Q> const& q, vec<4, float, Q> const& v)
		{
#			ifdef GLM_FORCE_QUAT_DATA_WXYZ
				__m128 const q_wwww = _mm_shuffle_ps(q.data, q.data, _MM_SHUFFLE(0, 0, 0, 0));
				__m128 const q_swp0 = _mm_shuffle_ps(q.data, q.data, _MM_SHUFFLE(0, 1, 3, 2));
				__m128 const q_swp1 = _mm_shuffle_ps(q.data, q.data, _MM_SHUFFLE(0, 2, 1, 3));
				__m128 const v_swp0 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(3, 0, 2, 1));
				__m128 const v_swp1 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(3, 1, 0, 2));

				__m128 uv = _mm_sub_ps(_mm_mul_ps(q_swp0, v_swp1), _mm_mul_ps(q_swp1, v_swp0));
				__m128 uv_swp0 = _mm_shuffle_ps(uv, uv, _MM_SHUFFLE(3, 0, 2, 1));
				__m128 uv_swp1 = _mm_shuffle_ps(uv, uv, _MM_SHUFFLE(3, 1, 0, 2));
				__m128 uuv = _mm_sub_ps(_mm_mul_ps(q_swp0, uv_swp1), _mm_mul_ps(q_swp1, uv_swp0));

				__m128 const two = _mm_set1_ps(2.0f);
				uv = _mm_mul_ps(uv, _mm_mul_ps(q_wwww, two));
				uuv = _mm_mul_ps(uuv, two);

				vec<4, float, Q> Result;
				Result.data = _mm_add_ps(v.data, _mm_add_ps(uv, uuv));
				return Result;
#			else
				__m128 const q_wwww = _mm_shuffle_ps(q.data, q.data, _MM_SHUFFLE(3, 3, 3, 3));
				__m128 const q_swp0 = _mm_shuffle_ps(q.data, q.data, _MM_SHUFFLE(3, 0, 2, 1));
				__m128 const q_swp1 = _mm_shuffle_ps(q.data, q.data, _MM_SHUFFLE(3, 1, 0, 2));
				__m128 const v_swp0 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(3, 0, 2, 1));
				__m128 const v_swp1 = _mm_shuffle_ps(v.data, v.data, _MM_SHUFFLE(3, 1, 0, 2));

				__m128 uv      = _mm_sub_ps(_mm_mul_ps(q_swp0, v_swp1), _mm_mul_ps(q_swp1, v_swp0));
				__m128 uv_swp0 = _mm_shuffle_ps(uv, uv, _MM_SHUFFLE(3, 0, 2, 1));
				__m128 uv_swp1 = _mm_shuffle_ps(uv, uv, _MM_SHUFFLE(3, 1, 0, 2));
				__m128 uuv     = _mm_sub_ps(_mm_mul_ps(q_swp0, uv_swp1), _mm_mul_ps(q_swp1, uv_swp0));

				__m128 const two = _mm_set1_ps(2.0f);
				uv  = _mm_mul_ps(uv, _mm_mul_ps(q_wwww, two));
				uuv = _mm_mul_ps(uuv, two);

				vec<4, float, Q> Result;
				Result.data = _mm_add_ps(v.data, _mm_add_ps(uv, uuv));
				return Result;
#			endif
		}
	};
}
}

#endif
