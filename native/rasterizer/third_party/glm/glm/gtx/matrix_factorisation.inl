namespace glm
{
	template <length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER mat<C, R, T, Q> flipud(mat<C, R, T, Q> const& in)
	{
		mat<R, C, T, Q> tin = transpose(in);
		tin = fliplr(tin);
		mat<C, R, T, Q> out = transpose(tin);

		return out;
	}

	template <length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER mat<C, R, T, Q> fliplr(mat<C, R, T, Q> const& in)
	{
		mat<C, R, T, Q> out;
		for (length_t i = 0; i < C; i++)
		{
			out[i] = in[(C - i) - 1];
		}

		return out;
	}

	template <length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER void qr_decompose(mat<C, R, T, Q> const& in, mat<(C < R ? C : R), R, T, Q>& q, mat<C, (C < R ? C : R), T, Q>& r)
	{

		for (length_t i = 0; i < (C < R ? C : R); i++)
		{
			
			q[i] = in[i];

			for (length_t j = 0; j < i; j++)
			{
				q[i] -= dot(q[i], q[j])*q[j];
				r[j][i] = 0;
			}

			q[i] = normalize(q[i]);

			for (length_t j = i; j < C; j++)
			{
				r[j][i] = dot(in[j], q[i]);
			}
		}
	}

	template <length_t C, length_t R, typename T, qualifier Q>
	GLM_FUNC_QUALIFIER void rq_decompose(mat<C, R, T, Q> const& in, mat<(C < R ? C : R), R, T, Q>& r, mat<C, (C < R ? C : R), T, Q>& q)
	{

		mat<R, C, T, Q> tin = transpose(in);
		tin = fliplr(tin);

		mat<R, (C < R ? C : R), T, Q> tr;
		mat<(C < R ? C : R), C, T, Q> tq;
		qr_decompose(tin, tq, tr);

		tr = fliplr(tr);
		r = transpose(tr);
		r = fliplr(r);

		tq = fliplr(tq);
		q = transpose(tq);
	}
} 
