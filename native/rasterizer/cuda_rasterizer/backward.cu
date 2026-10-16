#include "backward.h"
#include "auxiliary.h"
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cmath>
#include "math.h"
namespace cg = cooperative_groups;

__device__ void computeColorFromSH(int idx, int deg, int max_coeffs, const glm::vec3* means, glm::vec3 campos, const float* shs, const bool* clamped, const glm::vec3* dL_dcolor, glm::vec3* dL_dmeans, glm::vec3* dL_dshs, float *dL_dtau)
{
	
	glm::vec3 pos = means[idx];
	glm::vec3 dir_orig = pos - campos;
	glm::vec3 dir = dir_orig / glm::length(dir_orig);

	glm::vec3* sh = ((glm::vec3*)shs) + idx * max_coeffs;

	glm::vec3 dL_dRGB = dL_dcolor[idx];
	dL_dRGB.x *= clamped[3 * idx + 0] ? 0 : 1;
	dL_dRGB.y *= clamped[3 * idx + 1] ? 0 : 1;
	dL_dRGB.z *= clamped[3 * idx + 2] ? 0 : 1;

	glm::vec3 dRGBdx(0, 0, 0);
	glm::vec3 dRGBdy(0, 0, 0);
	glm::vec3 dRGBdz(0, 0, 0);
	float x = dir.x;
	float y = dir.y;
	float z = dir.z;

	glm::vec3* dL_dsh = dL_dshs + idx * max_coeffs;

	float dRGBdsh0 = SH_C0;
	dL_dsh[0] = dRGBdsh0 * dL_dRGB;
	if (deg > 0)
	{
		float dRGBdsh1 = -SH_C1 * y;
		float dRGBdsh2 = SH_C1 * z;
		float dRGBdsh3 = -SH_C1 * x;
		dL_dsh[1] = dRGBdsh1 * dL_dRGB;
		dL_dsh[2] = dRGBdsh2 * dL_dRGB;
		dL_dsh[3] = dRGBdsh3 * dL_dRGB;

		dRGBdx = -SH_C1 * sh[3];
		dRGBdy = -SH_C1 * sh[1];
		dRGBdz = SH_C1 * sh[2];

		if (deg > 1)
		{
			float xx = x * x, yy = y * y, zz = z * z;
			float xy = x * y, yz = y * z, xz = x * z;

			float dRGBdsh4 = SH_C2[0] * xy;
			float dRGBdsh5 = SH_C2[1] * yz;
			float dRGBdsh6 = SH_C2[2] * (2.f * zz - xx - yy);
			float dRGBdsh7 = SH_C2[3] * xz;
			float dRGBdsh8 = SH_C2[4] * (xx - yy);
			dL_dsh[4] = dRGBdsh4 * dL_dRGB;
			dL_dsh[5] = dRGBdsh5 * dL_dRGB;
			dL_dsh[6] = dRGBdsh6 * dL_dRGB;
			dL_dsh[7] = dRGBdsh7 * dL_dRGB;
			dL_dsh[8] = dRGBdsh8 * dL_dRGB;

			dRGBdx += SH_C2[0] * y * sh[4] + SH_C2[2] * 2.f * -x * sh[6] + SH_C2[3] * z * sh[7] + SH_C2[4] * 2.f * x * sh[8];
			dRGBdy += SH_C2[0] * x * sh[4] + SH_C2[1] * z * sh[5] + SH_C2[2] * 2.f * -y * sh[6] + SH_C2[4] * 2.f * -y * sh[8];
			dRGBdz += SH_C2[1] * y * sh[5] + SH_C2[2] * 2.f * 2.f * z * sh[6] + SH_C2[3] * x * sh[7];

			if (deg > 2)
			{
				float dRGBdsh9 = SH_C3[0] * y * (3.f * xx - yy);
				float dRGBdsh10 = SH_C3[1] * xy * z;
				float dRGBdsh11 = SH_C3[2] * y * (4.f * zz - xx - yy);
				float dRGBdsh12 = SH_C3[3] * z * (2.f * zz - 3.f * xx - 3.f * yy);
				float dRGBdsh13 = SH_C3[4] * x * (4.f * zz - xx - yy);
				float dRGBdsh14 = SH_C3[5] * z * (xx - yy);
				float dRGBdsh15 = SH_C3[6] * x * (xx - 3.f * yy);
				dL_dsh[9] = dRGBdsh9 * dL_dRGB;
				dL_dsh[10] = dRGBdsh10 * dL_dRGB;
				dL_dsh[11] = dRGBdsh11 * dL_dRGB;
				dL_dsh[12] = dRGBdsh12 * dL_dRGB;
				dL_dsh[13] = dRGBdsh13 * dL_dRGB;
				dL_dsh[14] = dRGBdsh14 * dL_dRGB;
				dL_dsh[15] = dRGBdsh15 * dL_dRGB;

				dRGBdx += (
					SH_C3[0] * sh[9] * 3.f * 2.f * xy +
					SH_C3[1] * sh[10] * yz +
					SH_C3[2] * sh[11] * -2.f * xy +
					SH_C3[3] * sh[12] * -3.f * 2.f * xz +
					SH_C3[4] * sh[13] * (-3.f * xx + 4.f * zz - yy) +
					SH_C3[5] * sh[14] * 2.f * xz +
					SH_C3[6] * sh[15] * 3.f * (xx - yy));

				dRGBdy += (
					SH_C3[0] * sh[9] * 3.f * (xx - yy) +
					SH_C3[1] * sh[10] * xz +
					SH_C3[2] * sh[11] * (-3.f * yy + 4.f * zz - xx) +
					SH_C3[3] * sh[12] * -3.f * 2.f * yz +
					SH_C3[4] * sh[13] * -2.f * xy +
					SH_C3[5] * sh[14] * -2.f * yz +
					SH_C3[6] * sh[15] * -3.f * 2.f * xy);

				dRGBdz += (
					SH_C3[1] * sh[10] * xy +
					SH_C3[2] * sh[11] * 4.f * 2.f * yz +
					SH_C3[3] * sh[12] * 3.f * (2.f * zz - xx - yy) +
					SH_C3[4] * sh[13] * 4.f * 2.f * xz +
					SH_C3[5] * sh[14] * (xx - yy));
			}
		}
	}

	glm::vec3 dL_ddir(glm::dot(dRGBdx, dL_dRGB), glm::dot(dRGBdy, dL_dRGB), glm::dot(dRGBdz, dL_dRGB));

	float3 dL_dmean = dnormvdv(float3{ dir_orig.x, dir_orig.y, dir_orig.z }, float3{ dL_ddir.x, dL_ddir.y, dL_ddir.z });

	dL_dmeans[idx] += glm::vec3(dL_dmean.x, dL_dmean.y, dL_dmean.z);

	dL_dtau[6 * idx + 0] += -dL_dmean.x;
	dL_dtau[6 * idx + 1] += -dL_dmean.y;
	dL_dtau[6 * idx + 2] += -dL_dmean.z;
}

__global__ void computeCov2DCUDA(int P,
	const float3* means,
	const int* radii,
	const float* cov3Ds,
	const float h_x, float h_y,
	const float tan_fovx, float tan_fovy,
	const float* view_matrix,
	const float* dL_dconics,
	const float2* dL_dray_planes,
	glm::vec3* dL_dmeans,
	float* dL_dcov,
	const float4* __restrict__ conic_opacity,
	float* dL_dopacity,
	float *dL_dtau)
{
	auto idx = cg::this_grid().thread_rank();
	if (idx >= P || !(radii[idx] > 0))
		return;

	const float* cov3D = cov3Ds + 6 * idx;

	float3 mean = means[idx];
	float3 dL_dconic = { dL_dconics[4 * idx], dL_dconics[4 * idx + 1], dL_dconics[4 * idx + 3] };
	const float4 conic = conic_opacity[idx];
	const float combined_opacity = conic.w;

	const float2 dL_dray_plane = dL_dray_planes[idx];

	float3 t = transformPoint4x3(mean, view_matrix);
	
	const float limx = 1.3f * tan_fovx;
	const float limy = 1.3f * tan_fovy;
	float txtz = t.x / t.z;
	float tytz = t.y / t.z;
	t.x = min(limx, max(-limx, txtz)) * t.z;
	t.y = min(limy, max(-limy, tytz)) * t.z;
	
	const float x_grad_mul = txtz < -limx || txtz > limx ? 0 : 1;
	const float y_grad_mul = tytz < -limy || tytz > limy ? 0 : 1;

	txtz = t.x / t.z;
	tytz = t.y / t.z;

	glm::mat3 J = glm::mat3(h_x / t.z, 0.0f, -(h_x * t.x) / (t.z * t.z),
		0.0f, h_y / t.z, -(h_y * t.y) / (t.z * t.z),
		0, 0, 0);

	glm::mat3 W = glm::mat3(
		view_matrix[0], view_matrix[4], view_matrix[8],
		view_matrix[1], view_matrix[5], view_matrix[9],
		view_matrix[2], view_matrix[6], view_matrix[10]);

	glm::mat3 Vrk = glm::mat3(
		cov3D[0], cov3D[1], cov3D[2],
		cov3D[1], cov3D[3], cov3D[4],
		cov3D[2], cov3D[4], cov3D[5]);

	glm::mat3 T = W * J;

	glm::mat3 cov2D = glm::transpose(T) * glm::transpose(Vrk) * T;

	const float det_0 = max(1e-6, cov2D[0][0] * cov2D[1][1] - cov2D[0][1] * cov2D[0][1]);
	const float det_1 = max(1e-6, cov2D[0][0] * cov2D[1][1] - cov2D[0][1] * cov2D[0][1]);
	
	const float coef = sqrt(det_0 / (det_1+1e-6) + 1e-6);

	glm::mat3 Vrk_eigen_vector;
	glm::vec3 Vrk_eigen_value;
	int D = glm_modification::findEigenvaluesSymReal(Vrk,Vrk_eigen_value,Vrk_eigen_vector);

	unsigned int min_id = Vrk_eigen_value[0]>Vrk_eigen_value[1]? (Vrk_eigen_value[1]>Vrk_eigen_value[2]?2:1):(Vrk_eigen_value[0]>Vrk_eigen_value[2]?2:0);

	glm::mat3 Vrk_inv;
	glm::vec3 eigenvector_min;
	bool well_conditioned = Vrk_eigen_value[min_id]>0.00000001;
	if(well_conditioned)
	{
		glm::mat3 diag = glm::mat3(1/Vrk_eigen_value[0],0,0,
									0,1/Vrk_eigen_value[1],0,
									0,0,1/Vrk_eigen_value[2]);
		Vrk_inv = Vrk_eigen_vector * diag * glm::transpose(Vrk_eigen_vector);
	}
	else
	{
		eigenvector_min = Vrk_eigen_vector[min_id];
		Vrk_inv = glm::outerProduct(eigenvector_min,eigenvector_min);
	}

	glm::mat3 cov_cam_inv = glm::transpose(W) * Vrk_inv * W;
	glm::vec3 uvh = {txtz, tytz, 1};
	glm::vec3 uvh_m = cov_cam_inv * uvh;
	glm::vec3 uvh_mn = glm::normalize(uvh_m);

	float u2 = txtz * txtz;
	float v2 = tytz * tytz;
	float uv = txtz * tytz;

	glm::mat3 dL_dVrk;
	glm::vec3 plane;
	float dL_du;
	float dL_dv;
	float dL_dl;
	float l;
	float nl;
	if(isnan(uvh_mn.x)||D==0)
	{
		dL_dVrk = glm::mat3(0,0,0,0,0,0,0,0,0);
		plane = glm::vec3(0,0,0);
		nl = 1;
		l = 1;
		dL_du = 0;
		dL_dv = 0;
		dL_dl = 0;
	}
	else
	{
		float vb = glm::dot(uvh_m, uvh);
		float vbn = glm::dot(uvh_mn, uvh);

		l = sqrt(t.x*t.x+t.y*t.y+t.z*t.z);
		glm::mat3 nJ = glm::mat3(
			1 / t.z, 0.0f, -(t.x) / (t.z * t.z),
			0.0f, 1 / t.z, -(t.y) / (t.z * t.z),
			t.x/l, t.y/l, t.z/l);

		glm::mat3 nJ_inv = glm::mat3(
			v2 + 1,	-uv, 		0,
			-uv,	u2 + 1,		0,
			-txtz,	-tytz,		0
		);

		float clamp_vb = max(vb, 0.0000001f);
		float clamp_vbn = max(vbn, 0.0000001f);
		nl = u2+v2+1;
		float factor_normal = l / nl;
		glm::vec3 uvh_m_vb = uvh_mn/clamp_vbn;
		plane = nJ_inv * uvh_m_vb;

		glm::vec2 ray_plane = {plane[0]*factor_normal, plane[1]*factor_normal};
		glm::vec3 ray_normal_vector = {-plane[0]*factor_normal, -plane[1]*factor_normal, -1};

		glm::vec3 cam_normal_vector = nJ * ray_normal_vector;
		glm::vec3 normal_vector = glm::normalize(cam_normal_vector);
		float lv = glm::length(cam_normal_vector);
		dL_dl = (plane[0] * dL_dray_plane.x + plane[1] * dL_dray_plane.y) / nl;
		
		glm::vec2 dL_dplane = glm::vec2(
			(dL_dray_plane.x * l) / nl,
			(dL_dray_plane.y * l) / nl
		);
		glm::vec3 dL_dplane_append = glm::vec3(dL_dplane.x, dL_dplane.y, 0);

		float dL_dnl = (-dL_dray_plane.x * ray_plane.x - dL_dray_plane.y * ray_plane.y) / nl;

		float tmp = dL_dplane.x * plane.x + dL_dplane.y * plane.y;

		glm::vec3 W_uvh = W * uvh;

		if(well_conditioned){
			dL_dVrk = - glm::outerProduct(Vrk_inv * W_uvh, (Vrk_inv/clamp_vb) * (W_uvh * (-tmp) + W * glm::transpose(nJ_inv) * dL_dplane_append));
		}
		else{
			dL_dVrk = glm::mat3(0,0,0,0,0,0,0,0,0);
			float dL_dvb = -tmp / clamp_vb;
			glm::vec3 nJ_inv_dL_dplane = glm::transpose(nJ_inv) * glm::vec3(dL_dplane.x / clamp_vb, dL_dplane.y / clamp_vb, 0);
			glm::mat3 dL_dVrk_inv = glm::outerProduct(W_uvh, W_uvh * dL_dvb + W * nJ_inv_dL_dplane);
			glm::vec3 dL_dv = (dL_dVrk_inv + glm::transpose(dL_dVrk_inv)) * eigenvector_min;
			for(int j =0;j<3;j++)
			{
				if(j!=min_id)
				{
					float scale = glm::dot(Vrk_eigen_vector[j],dL_dv)/min(Vrk_eigen_value[min_id] - Vrk_eigen_value[j], - 0.0000001f);
					dL_dVrk += glm::outerProduct(Vrk_eigen_vector[j] * scale, eigenvector_min);
				}
			}
		}

		glm::vec3 dL_duvh = 2 * (-tmp) * uvh_m_vb + (cov_cam_inv/clamp_vb) * glm::transpose(nJ_inv) * dL_dplane_append;
		
		glm::mat3 dL_dnJ_inv = glm::outerProduct(dL_dplane_append, uvh_m_vb);
		
		dL_du = dL_dnl * 2 * txtz
				+ dL_duvh.x 
				+ (dL_dnJ_inv[0][1] + dL_dnJ_inv[1][0]) * (-tytz) + 2 * dL_dnJ_inv[1][1] * txtz - dL_dnJ_inv[2][0];
		dL_dv = dL_dnl * 2 * tytz
				+ dL_duvh.y 
				+ (dL_dnJ_inv[0][1] + dL_dnJ_inv[1][0]) * (-txtz) + 2 * dL_dnJ_inv[0][0] * tytz - dL_dnJ_inv[2][1];
	}

	const float opacity = combined_opacity / (coef + 1e-6);
	const float dL_dcoef = dL_dopacity[idx] * opacity;
	const float dL_dsqrtcoef = dL_dcoef * 0.5 * 1. / (coef + 1e-6);
	const float dL_ddet0 = dL_dsqrtcoef / (det_1+1e-6);
	const float dL_ddet1 = dL_dsqrtcoef * det_0 * (-1.f / (det_1 * det_1 + 1e-6));
	
	const float dcoef_da = dL_ddet0 * cov2D[1][1] + dL_ddet1 * (cov2D[1][1]);
	const float dcoef_db = dL_ddet0 * (-2. * cov2D[0][1]) + dL_ddet1 * (-2. * cov2D[0][1]);
	const float dcoef_dc = dL_ddet0 * cov2D[0][0] + dL_ddet1 * (cov2D[0][0]);
	
	float a = cov2D[0][0];
	float b = cov2D[0][1];
	float c = cov2D[1][1];

	float denom = a * c - b * b;
	float dL_da = 0, dL_db = 0, dL_dc = 0;
	float denom2inv = 1.0f / ((denom * denom) + 0.0000001f);

	if (denom2inv != 0)
	{

		dL_da = denom2inv * (-c * c * dL_dconic.x + 2 * b * c * dL_dconic.y + (denom - a * c) * dL_dconic.z);
		dL_dc = denom2inv * (-a * a * dL_dconic.z + 2 * a * b * dL_dconic.y + (denom - a * c) * dL_dconic.x);
		dL_db = denom2inv * 2 * (b * c * dL_dconic.x - (denom + 2 * b * b) * dL_dconic.y + a * b * dL_dconic.z);

		if (det_0 <= 1e-6 || det_1 <= 1e-6){
			dL_dopacity[idx] = 0;
		} else {
			
			dL_da += dcoef_da;
			dL_dc += dcoef_dc;
			dL_db += dcoef_db;

			dL_dopacity[idx] = dL_dopacity[idx] * coef;
		}

		dL_dcov[6 * idx + 0] = (T[0][0] * T[0][0] * dL_da + T[0][0] * T[1][0] * dL_db + T[1][0] * T[1][0] * dL_dc);
		dL_dcov[6 * idx + 3] = (T[0][1] * T[0][1] * dL_da + T[0][1] * T[1][1] * dL_db + T[1][1] * T[1][1] * dL_dc);
		dL_dcov[6 * idx + 5] = (T[0][2] * T[0][2] * dL_da + T[0][2] * T[1][2] * dL_db + T[1][2] * T[1][2] * dL_dc);

		dL_dcov[6 * idx + 1] = 2 * T[0][0] * T[0][1] * dL_da + (T[0][0] * T[1][1] + T[0][1] * T[1][0]) * dL_db + 2 * T[1][0] * T[1][1] * dL_dc;
		dL_dcov[6 * idx + 2] = 2 * T[0][0] * T[0][2] * dL_da + (T[0][0] * T[1][2] + T[0][2] * T[1][0]) * dL_db + 2 * T[1][0] * T[1][2] * dL_dc;
		dL_dcov[6 * idx + 4] = 2 * T[0][2] * T[0][1] * dL_da + (T[0][1] * T[1][2] + T[0][2] * T[1][1]) * dL_db + 2 * T[1][1] * T[1][2] * dL_dc;
	}
	else
	{
		for (int i = 0; i < 6; i++)
			dL_dcov[6 * idx + i] = 0;
	}
	dL_dcov[6 * idx + 0] += dL_dVrk[0][0];
	dL_dcov[6 * idx + 3] += dL_dVrk[1][1];
	dL_dcov[6 * idx + 5] += dL_dVrk[2][2];
	dL_dcov[6 * idx + 1] += dL_dVrk[0][1] + dL_dVrk[1][0];
	dL_dcov[6 * idx + 2] += dL_dVrk[0][2] + dL_dVrk[2][0];
	dL_dcov[6 * idx + 4] += dL_dVrk[1][2] + dL_dVrk[2][1];

	float dL_dT00 = 2 * (T[0][0] * Vrk[0][0] + T[0][1] * Vrk[0][1] + T[0][2] * Vrk[0][2]) * dL_da +
		(T[1][0] * Vrk[0][0] + T[1][1] * Vrk[0][1] + T[1][2] * Vrk[0][2]) * dL_db;
	float dL_dT01 = 2 * (T[0][0] * Vrk[1][0] + T[0][1] * Vrk[1][1] + T[0][2] * Vrk[1][2]) * dL_da +
		(T[1][0] * Vrk[1][0] + T[1][1] * Vrk[1][1] + T[1][2] * Vrk[1][2]) * dL_db;
	float dL_dT02 = 2 * (T[0][0] * Vrk[2][0] + T[0][1] * Vrk[2][1] + T[0][2] * Vrk[2][2]) * dL_da +
		(T[1][0] * Vrk[2][0] + T[1][1] * Vrk[2][1] + T[1][2] * Vrk[2][2]) * dL_db;
	float dL_dT10 = 2 * (T[1][0] * Vrk[0][0] + T[1][1] * Vrk[0][1] + T[1][2] * Vrk[0][2]) * dL_dc +
		(T[0][0] * Vrk[0][0] + T[0][1] * Vrk[0][1] + T[0][2] * Vrk[0][2]) * dL_db;
	float dL_dT11 = 2 * (T[1][0] * Vrk[1][0] + T[1][1] * Vrk[1][1] + T[1][2] * Vrk[1][2]) * dL_dc +
		(T[0][0] * Vrk[1][0] + T[0][1] * Vrk[1][1] + T[0][2] * Vrk[1][2]) * dL_db;
	float dL_dT12 = 2 * (T[1][0] * Vrk[2][0] + T[1][1] * Vrk[2][1] + T[1][2] * Vrk[2][2]) * dL_dc +
		(T[0][0] * Vrk[2][0] + T[0][1] * Vrk[2][1] + T[0][2] * Vrk[2][2]) * dL_db;

	float dL_dJ00 = W[0][0] * dL_dT00 + W[0][1] * dL_dT01 + W[0][2] * dL_dT02;
	float dL_dJ02 = W[2][0] * dL_dT00 + W[2][1] * dL_dT01 + W[2][2] * dL_dT02;
	float dL_dJ11 = W[1][0] * dL_dT10 + W[1][1] * dL_dT11 + W[1][2] * dL_dT12;
	float dL_dJ12 = W[2][0] * dL_dT10 + W[2][1] * dL_dT11 + W[2][2] * dL_dT12;

	float tz = 1.f / t.z;
	float tz2 = tz * tz;
	float tz3 = tz2 * tz;

	float l3 = l * l * l;
	float dL_dtx = x_grad_mul * (-h_x * tz2 * dL_dJ02 + dL_du * tz
								+dL_dl*t.x/l);
	float dL_dty = y_grad_mul * (-h_y * tz2 * dL_dJ12 + dL_dv * tz
								+dL_dl*t.y/l);
	float dL_dtz = -h_x * tz2 * dL_dJ00 - h_y * tz2 * dL_dJ11 + (2 * h_x * t.x) * tz3 * dL_dJ02 + (2 * h_y * t.y) * tz3 * dL_dJ12
					- (dL_du * t.x + dL_dv * t.y) * tz2
					+ dL_dl*t.z/l;

	SE3 T_CW(view_matrix);
	mat33 R = T_CW.R().data();
	mat33 RT = R.transpose();
	float3 t_ = T_CW.t();
	mat33 dpC_drho = mat33::identity();
	mat33 dpC_dtheta = -mat33::skew_symmetric(t);
	float dL_dt[6];
	for (int i = 0; i < 3; i++) {
		float3 c_rho = dpC_drho.cols[i];
		float3 c_theta = dpC_dtheta.cols[i];
		dL_dt[i] = dL_dtx * c_rho.x + dL_dty * c_rho.y + dL_dtz * c_rho.z;
		dL_dt[i + 3] = dL_dtx * c_theta.x + dL_dty * c_theta.y + dL_dtz * c_theta.z;
	}
	for (int i = 0; i < 6; i++) {
		dL_dtau[6 * idx + i] += dL_dt[i];
	}

	float3 dL_dmean = transformVec4x3Transpose({ dL_dtx, dL_dty, dL_dtz }, view_matrix);

	dL_dmeans[idx] = glm::vec3(dL_dmean.x, dL_dmean.y, dL_dmean.z);
	
	float dL_dW00 = J[0][0] * dL_dT00;
	float dL_dW01 = J[0][0] * dL_dT01;
	float dL_dW02 = J[0][0] * dL_dT02;
	float dL_dW10 = J[1][1] * dL_dT10;
	float dL_dW11 = J[1][1] * dL_dT11;
	float dL_dW12 = J[1][1] * dL_dT12;
	float dL_dW20 = J[0][2] * dL_dT00 + J[1][2] * dL_dT10;
	float dL_dW21 = J[0][2] * dL_dT01 + J[1][2] * dL_dT11;
	float dL_dW22 = J[0][2] * dL_dT02 + J[1][2] * dL_dT12;

	float3 c1 = R.cols[0];
	float3 c2 = R.cols[1];
	float3 c3 = R.cols[2];

	float dL_dW_data[9];
	dL_dW_data[0] = dL_dW00;
	dL_dW_data[3] = dL_dW01;
	dL_dW_data[6] = dL_dW02;
	dL_dW_data[1] = dL_dW10;
	dL_dW_data[4] = dL_dW11;
	dL_dW_data[7] = dL_dW12;
	dL_dW_data[2] = dL_dW20;
	dL_dW_data[5] = dL_dW21;
	dL_dW_data[8] = dL_dW22;

	mat33 dL_dW(dL_dW_data);
	float3 dL_dWc1 = dL_dW.cols[0];
	float3 dL_dWc2 = dL_dW.cols[1];
	float3 dL_dWc3 = dL_dW.cols[2];

	mat33 n_W1_x = -mat33::skew_symmetric(c1);
	mat33 n_W2_x = -mat33::skew_symmetric(c2);
	mat33 n_W3_x = -mat33::skew_symmetric(c3);

	float3 dL_dtheta = {};
	dL_dtheta.x = dot(dL_dWc1, n_W1_x.cols[0]) + dot(dL_dWc2, n_W2_x.cols[0]) +
				dot(dL_dWc3, n_W3_x.cols[0]);
	dL_dtheta.y = dot(dL_dWc1, n_W1_x.cols[1]) + dot(dL_dWc2, n_W2_x.cols[1]) +
				dot(dL_dWc3, n_W3_x.cols[1]);
	dL_dtheta.z = dot(dL_dWc1, n_W1_x.cols[2]) + dot(dL_dWc2, n_W2_x.cols[2]) +
				dot(dL_dWc3, n_W3_x.cols[2]);

	dL_dtau[6 * idx + 3] += dL_dtheta.x;
	dL_dtau[6 * idx + 4] += dL_dtheta.y;
	dL_dtau[6 * idx + 5] += dL_dtheta.z;
}

__device__ void computeCov3D(int idx, const glm::vec3 scale, float mod, const glm::vec4 rot, const float* dL_dcov3Ds, glm::vec3* dL_dscales, glm::vec4* dL_drots)
{
	
	glm::vec4 q = rot;
	float r = q.x;
	float x = q.y;
	float y = q.z;
	float z = q.w;

	glm::mat3 R = glm::mat3(
		1.f - 2.f * (y * y + z * z), 2.f * (x * y - r * z), 2.f * (x * z + r * y),
		2.f * (x * y + r * z), 1.f - 2.f * (x * x + z * z), 2.f * (y * z - r * x),
		2.f * (x * z - r * y), 2.f * (y * z + r * x), 1.f - 2.f * (x * x + y * y)
	);

	glm::mat3 S = glm::mat3(1.0f);

	glm::vec3 s = mod * scale;
	S[0][0] = s.x;
	S[1][1] = s.y;
	S[2][2] = s.z;

	glm::mat3 M = S * R;

	const float* dL_dcov3D = dL_dcov3Ds + 6 * idx;

	glm::vec3 dunc(dL_dcov3D[0], dL_dcov3D[3], dL_dcov3D[5]);
	glm::vec3 ounc = 0.5f * glm::vec3(dL_dcov3D[1], dL_dcov3D[2], dL_dcov3D[4]);

	glm::mat3 dL_dSigma = glm::mat3(
		dL_dcov3D[0], 0.5f * dL_dcov3D[1], 0.5f * dL_dcov3D[2],
		0.5f * dL_dcov3D[1], dL_dcov3D[3], 0.5f * dL_dcov3D[4],
		0.5f * dL_dcov3D[2], 0.5f * dL_dcov3D[4], dL_dcov3D[5]
	);

	glm::mat3 dL_dM = 2.0f * M * dL_dSigma;

	glm::mat3 Rt = glm::transpose(R);
	glm::mat3 dL_dMt = glm::transpose(dL_dM);

	glm::vec3* dL_dscale = dL_dscales + idx;
	dL_dscale->x = glm::dot(Rt[0], dL_dMt[0]);
	dL_dscale->y = glm::dot(Rt[1], dL_dMt[1]);
	dL_dscale->z = glm::dot(Rt[2], dL_dMt[2]);

	dL_dMt[0] *= s.x;
	dL_dMt[1] *= s.y;
	dL_dMt[2] *= s.z;

	glm::vec4 dL_dq;
	dL_dq.x = 2 * z * (dL_dMt[0][1] - dL_dMt[1][0]) + 2 * y * (dL_dMt[2][0] - dL_dMt[0][2]) + 2 * x * (dL_dMt[1][2] - dL_dMt[2][1]);
	dL_dq.y = 2 * y * (dL_dMt[1][0] + dL_dMt[0][1]) + 2 * z * (dL_dMt[2][0] + dL_dMt[0][2]) + 2 * r * (dL_dMt[1][2] - dL_dMt[2][1]) - 4 * x * (dL_dMt[2][2] + dL_dMt[1][1]);
	dL_dq.z = 2 * x * (dL_dMt[1][0] + dL_dMt[0][1]) + 2 * r * (dL_dMt[2][0] - dL_dMt[0][2]) + 2 * z * (dL_dMt[1][2] + dL_dMt[2][1]) - 4 * y * (dL_dMt[2][2] + dL_dMt[0][0]);
	dL_dq.w = 2 * r * (dL_dMt[0][1] - dL_dMt[1][0]) + 2 * x * (dL_dMt[2][0] + dL_dMt[0][2]) + 2 * y * (dL_dMt[1][2] + dL_dMt[2][1]) - 4 * z * (dL_dMt[1][1] + dL_dMt[0][0]);

	float4* dL_drot = (float4*)(dL_drots + idx);
	*dL_drot = float4{ dL_dq.x, dL_dq.y, dL_dq.z, dL_dq.w };
}

template<int C>
__global__ void preprocessCUDA(
	int P, int D, int M,
	const float3* means,
	const int* radii,
	const float* shs,
	const bool* clamped,
	const glm::vec3* scales,
	const glm::vec4* rotations,
	const float scale_modifier,
	const float* view,
	const float* proj,
	const float *proj_raw,
	const glm::vec3* campos,
	const float4* dL_dmean2D,
	glm::vec3* dL_dmeans,
	float* dL_dcolor,
	const float* dL_dts,
	float* dL_dcov3D,
	float* dL_dsh,
	glm::vec3* dL_dscale,
	glm::vec4* dL_drot,
	float *dL_dtau)
{
	auto idx = cg::this_grid().thread_rank();
	if (idx >= P || !(radii[idx] > 0))
		return;

	float3 m = means[idx];

	float4 m_hom = transformPoint4x4(m, proj);
	float m_w = 1.0f / (m_hom.w + 0.0000001f);

	glm::vec3 dL_dmean1;
	float mul1 = (proj[0] * m.x + proj[4] * m.y + proj[8] * m.z + proj[12]) * m_w * m_w;
	float mul2 = (proj[1] * m.x + proj[5] * m.y + proj[9] * m.z + proj[13]) * m_w * m_w;
	dL_dmean1.x = (proj[0] * m_w - proj[3] * mul1) * dL_dmean2D[idx].x + (proj[1] * m_w - proj[3] * mul2) * dL_dmean2D[idx].y;
	dL_dmean1.y = (proj[4] * m_w - proj[7] * mul1) * dL_dmean2D[idx].x + (proj[5] * m_w - proj[7] * mul2) * dL_dmean2D[idx].y;
	dL_dmean1.z = (proj[8] * m_w - proj[11] * mul1) * dL_dmean2D[idx].x + (proj[9] * m_w - proj[11] * mul2) * dL_dmean2D[idx].y;

	dL_dmeans[idx] += dL_dmean1;

	float alpha = 1.0f * m_w;
	float beta = -m_hom.x * m_w * m_w;
	float gamma = -m_hom.y * m_w * m_w;

	float a = proj_raw[0];
	float b = proj_raw[5];
	float c = proj_raw[10];
	float d = proj_raw[14];
	float e = proj_raw[11];

	SE3 T_CW(view);
	mat33 R = T_CW.R().data();
	mat33 RT = R.transpose();
	
	float3 p_C = T_CW * m;
	mat33 dp_C_d_rho = mat33::identity();
	mat33 dp_C_d_theta = -mat33::skew_symmetric(p_C);

	float3 d_proj_dp_C1 = make_float3(alpha * a, 0.f, beta * e);
	float3 d_proj_dp_C2 = make_float3(0.f, alpha * b, gamma * e);

	float3 d_proj_dp_C1_d_rho = dp_C_d_rho.transpose() * d_proj_dp_C1; 
	float3 d_proj_dp_C2_d_rho = dp_C_d_rho.transpose() * d_proj_dp_C2;
	float3 d_proj_dp_C1_d_theta = dp_C_d_theta.transpose() * d_proj_dp_C1;
	float3 d_proj_dp_C2_d_theta = dp_C_d_theta.transpose() * d_proj_dp_C2;

	float2 dmean2D_dtau[6];
	dmean2D_dtau[0].x = d_proj_dp_C1_d_rho.x;
	dmean2D_dtau[1].x = d_proj_dp_C1_d_rho.y;
	dmean2D_dtau[2].x = d_proj_dp_C1_d_rho.z;
	dmean2D_dtau[3].x = d_proj_dp_C1_d_theta.x;
	dmean2D_dtau[4].x = d_proj_dp_C1_d_theta.y;
	dmean2D_dtau[5].x = d_proj_dp_C1_d_theta.z;

	dmean2D_dtau[0].y = d_proj_dp_C2_d_rho.x;
	dmean2D_dtau[1].y = d_proj_dp_C2_d_rho.y;
	dmean2D_dtau[2].y = d_proj_dp_C2_d_rho.z;
	dmean2D_dtau[3].y = d_proj_dp_C2_d_theta.x;
	dmean2D_dtau[4].y = d_proj_dp_C2_d_theta.y;
	dmean2D_dtau[5].y = d_proj_dp_C2_d_theta.z;

	float dL_dt2[6];
	for (int i = 0; i < 6; i++) {
		dL_dt2[i] = dL_dmean2D[idx].x * dmean2D_dtau[i].x + dL_dmean2D[idx].y * dmean2D_dtau[i].y;
	}
	for (int i = 0; i < 6; i++) {
		dL_dtau[6 * idx + i] += dL_dt2[i];
	}

	float3 m_view = transformPoint4x3(m, view);
	float t = sqrt(m_view.x*m_view.x+m_view.y*m_view.y+m_view.z*m_view.z);
	float dL_dt = dL_dts[idx];
	float3 dL_dm_view = {m_view.x/t*dL_dt, m_view.y/t*dL_dt, m_view.z/t*dL_dt};
	float3 dL_dmean2 = transformVec4x3Transpose(dL_dm_view, view);
	
	dL_dmeans[idx] += glm::vec3(
		dL_dmean2.x,
		dL_dmean2.y,
		dL_dmean2.z);

	for (int i = 0; i < 3; i++) {
		float3 c_rho = dp_C_d_rho.cols[i];
		float3 c_theta = dp_C_d_theta.cols[i];
		dL_dtau[6 * idx + i] += dL_dm_view.x * c_rho.x;
		dL_dtau[6 * idx + i] += dL_dm_view.y * c_rho.y;
		dL_dtau[6 * idx + i] += dL_dm_view.z * c_rho.z;
		dL_dtau[6 * idx + i + 3] += dL_dm_view.x * c_theta.x;
		dL_dtau[6 * idx + i + 3] += dL_dm_view.y * c_theta.y;
		dL_dtau[6 * idx + i + 3] += dL_dm_view.z * c_theta.z;
	}

	if (shs)
		computeColorFromSH(idx, D, M, (glm::vec3*)means, *campos, shs, clamped, (glm::vec3*)dL_dcolor, (glm::vec3*)dL_dmeans, (glm::vec3*)dL_dsh, dL_dtau);

	if (scales)
		computeCov3D(idx, scales[idx], scale_modifier, rotations[idx], dL_dcov3D, dL_dscale, dL_drot);
}

template <typename T>
__device__ void inline reduce_helper(int lane, int i, T *data) {
  if (lane < i) {
    data[lane] += data[lane + i];
  }
}

template <typename group_t, typename... Lists>
__device__ void render_cuda_reduce_sum(group_t g, Lists... lists) {
  int lane = g.thread_rank();
  g.sync();

  for (int i = g.size() / 2; i > 0; i /= 2) {
    (...,
     reduce_helper(
         lane, i, lists)); 
    g.sync();
  }
}

template <uint32_t C, bool DEPTH = true>
__global__ void __launch_bounds__(BLOCK_X * BLOCK_Y)
renderCUDA(
	const uint2* __restrict__ ranges,
	const uint32_t* __restrict__ point_list,
	int W, int H,
	const float* __restrict__ bg_color,
	const float* __restrict__ view_points,
	const float2* __restrict__ points_xy_image,
	const float4* __restrict__ conic_opacity,
	const float* __restrict__ colors,
	const float* __restrict__ depths,
	const float* __restrict__ ts,
	const float2* __restrict__ ray_planes,
	const float* __restrict__ alphas,
	const float* __restrict__ accum_depth,
	const uint32_t* __restrict__ n_contrib,
	const float* __restrict__ dL_dpixels,
	const float* __restrict__ dL_dpixel_depths,
	const float focal_x, 
	const float focal_y,
	float4* __restrict__ dL_dmean2D,
	float4* __restrict__ dL_dconic2D,
	float* __restrict__ dL_dopacity,
	float* __restrict__ dL_dcolors,
	float* __restrict__ dL_dts,
	float2* __restrict__ dL_dray_planes)
{
	
	auto block = cg::this_thread_block();
	auto tid = block.thread_rank();

	const uint32_t horizontal_blocks = (W + BLOCK_X - 1) / BLOCK_X;
	const uint2 pix_min = { block.group_index().x * BLOCK_X, block.group_index().y * BLOCK_Y };
	const uint2 pix_max = { min(pix_min.x + BLOCK_X, W), min(pix_min.y + BLOCK_Y , H) };
	const uint2 pix = { pix_min.x + block.thread_index().x, pix_min.y + block.thread_index().y };
	const uint32_t pix_id = W * pix.y + pix.x;
	const float2 pixf = { (float)pix.x, (float)pix.y };
	const float2 pixnf = {(pixf.x-W/2.f)/focal_x,(pixf.y-H/2.f)/focal_y};
	const float ln = sqrt(pixnf.x*pixnf.x+pixnf.y*pixnf.y+1);

	const bool inside = pix.x < W&& pix.y < H;
	const uint2 range = ranges[block.group_index().y * horizontal_blocks + block.group_index().x];

	const int rounds = ((range.y - range.x + BLOCK_SIZE - 1) / BLOCK_SIZE);

	bool done = !inside;
	int toDo = range.y - range.x;

	__shared__ int collected_id[BLOCK_SIZE];
	__shared__ float2 collected_xy[BLOCK_SIZE];
	__shared__ float4 collected_conic_opacity[BLOCK_SIZE];
	__shared__ float collected_colors[C * BLOCK_SIZE];
	__shared__ float collected_mean3d[3 * BLOCK_SIZE];
	__shared__ float collected_ts[BLOCK_SIZE];
	__shared__ float2 collected_ray_planes[BLOCK_SIZE];

	__shared__ float2 dL_dmean2D_shared[BLOCK_SIZE];
	__shared__ float4 dL_dconic2D_shared[BLOCK_SIZE];
	__shared__ float dL_dopacity_shared[BLOCK_SIZE];
	__shared__ float3 dL_dcolors_shared[BLOCK_SIZE];
	__shared__ float dL_dts_shared[BLOCK_SIZE];
	__shared__ float2 dL_dray_planes_shared[BLOCK_SIZE];

	const float T_final = inside ? (1 - alphas[pix_id]) : 0;
	const float w_final = inside ? alphas[pix_id] : 0;

	float T = T_final;
	float w = w_final;

	uint32_t contributor = toDo;
	const int last_contributor = inside ? n_contrib[pix_id] : 0;

	float accum_rec[C] = { 0 };
	float dL_dpixel[C];
	float accum_t_rec = 0;
	float dL_dpixel_t;
	float accum_alpha_rec = 0;
	float dL_dalpha= 0;
	
	if (inside) {
		#pragma unroll
		for (int i = 0; i < C; i++) {
			dL_dpixel[i] = dL_dpixels[i * H * W + pix_id];
		}

		float ww = w_final*w_final;
		if constexpr (DEPTH)
		{
			float dL_dpixel_depth_w = dL_dpixel_depths[pix_id];
			float pixel_accum_depth = accum_depth[pix_id];
			dL_dalpha -= dL_dpixel_depth_w*pixel_accum_depth/ww;
			dL_dpixel_t = dL_dpixel_depth_w / w_final/ ln;
		}
	}

	float last_alpha = 0;
	float last_color[C] = { 0 };
	float last_t = 0;
	float last_dL_dw = 0;

	const float ddelx_dx = 0.5 * W;
	const float ddely_dy = 0.5 * H;
	__shared__ int skip_counter;

	for (int i = 0; i < rounds; i++, toDo -= BLOCK_SIZE)
	{

		const int progress = i * BLOCK_SIZE + tid;
		if (range.x + progress < range.y)
		{
			const int coll_id = point_list[range.y - progress - 1];
			collected_id[tid] = coll_id;
			collected_xy[tid] = points_xy_image[coll_id];
			collected_conic_opacity[tid] = conic_opacity[coll_id];
			#pragma unroll
			for (int i = 0; i < C; i++) {
				collected_colors[i * BLOCK_SIZE + tid] = colors[coll_id * C + i];
			}
			if constexpr (DEPTH)
			{
				collected_ray_planes[tid] = ray_planes[coll_id];
				collected_ts[tid] = ts[coll_id];
			}
		}
		
		for (int j = 0; j < min(BLOCK_SIZE, toDo); j++)
		{
			block.sync();
			if (tid == 0) {
				skip_counter = 0;
			}
			block.sync();

			bool skip = done;
			contributor = done ? contributor : contributor - 1;
			skip |= contributor >= last_contributor;

			const float2 xy = collected_xy[j];
			const float2 d = { xy.x - pixf.x, xy.y - pixf.y };
			const float4 con_o = collected_conic_opacity[j];

			float power = -0.5f * (con_o.x * d.x * d.x + con_o.z * d.y * d.y) - con_o.y * d.x * d.y;
			skip |= power > 0.0f;

			const float G = exp(power);
			const float alpha = min(0.99f, con_o.w * G);
			skip |= alpha < 1.0f / 255.0f;

			if (skip) {
				atomicAdd(&skip_counter, 1);
			}
			block.sync();
			if (skip_counter == BLOCK_SIZE) {
				continue;
			}

			T = skip ? T : T / (1.f - alpha);
			const float dchannel_dcolor = alpha * T;
			const float& dpixel_t_dt = dchannel_dcolor;

			float dL_dopa = 0.0f;
			const int global_id = collected_id[j];
			float local_dL_dcolors[3];
			#pragma unroll
			for (int ch = 0; ch < C; ch++)
			{
				const float c = collected_colors[ch * BLOCK_SIZE + j];
				
				accum_rec[ch] = skip ? accum_rec[ch] : last_alpha * last_color[ch] + (1.f - last_alpha) * accum_rec[ch];
				last_color[ch] = skip ? last_color[ch] : c;

				const float dL_dchannel = dL_dpixel[ch];
				dL_dopa += (c - accum_rec[ch]) * dL_dchannel;
				local_dL_dcolors[ch] = skip ? 0.0f : dchannel_dcolor * dL_dchannel;
			}
			dL_dcolors_shared[tid].x = local_dL_dcolors[0];
			dL_dcolors_shared[tid].y = local_dL_dcolors[1];
			dL_dcolors_shared[tid].z = local_dL_dcolors[2];

			float dL_dt;
			float2 ray_plane;

			if constexpr (DEPTH)
			{
				const float t_center = collected_ts[j];
				ray_plane = collected_ray_planes[j];
				float t = t_center + (ray_plane.x * d.x + ray_plane.y * d.y);
				accum_t_rec = last_alpha * last_t + (1.f - last_alpha) * accum_t_rec;
				last_t = t;
				dL_dopa += (t - accum_t_rec) * dL_dpixel_t;
				dL_dt = dpixel_t_dt * dL_dpixel_t;
				
				dL_dts_shared[tid] = skip ? 0.f : dL_dt;
				dL_dray_planes_shared[tid].x = skip ? 0.f : dL_dt * d.x / focal_x;
				dL_dray_planes_shared[tid].y = skip ? 0.f : dL_dt * d.y / focal_y;
			}

			accum_alpha_rec = last_alpha + (1.f - last_alpha) * accum_alpha_rec;
			dL_dopa += (1 - accum_alpha_rec) * dL_dalpha;

			dL_dopa *= T;
			
			last_alpha = skip ? last_alpha : alpha;

			float bg_dot_dpixel = 0;
			#pragma unroll
			for (int i = 0; i < C; i++) {
				bg_dot_dpixel += bg_color[i] * dL_dpixel[i];
			}
			dL_dopa += (-T_final / (1.f - alpha)) * bg_dot_dpixel;

			const float dL_dG = con_o.w * dL_dopa;
			const float gdx = G * d.x;
			const float gdy = G * d.y;
			const float dG_ddelx = -gdx * con_o.x - gdy * con_o.y;
			const float dG_ddely = -gdy * con_o.z - gdx * con_o.y;

			float dL_ddelx = dL_dG * dG_ddelx;
			float dL_ddely = dL_dG * dG_ddely;
			if constexpr (DEPTH)
			{
				dL_ddelx += dL_dt * ray_plane.x;
				dL_ddely += dL_dt * ray_plane.y;
			}

			dL_dmean2D_shared[tid].x = skip ? 0.f : dL_ddelx * ddelx_dx;
			dL_dmean2D_shared[tid].y = skip ? 0.f : dL_ddely * ddely_dy;
            
            atomicAdd(&dL_dmean2D[global_id].z, fabsf(dL_dmean2D_shared[tid].x));
            atomicAdd(&dL_dmean2D[global_id].w, fabsf(dL_dmean2D_shared[tid].y));
			dL_dconic2D_shared[tid].x = skip ? 0.f : -0.5f * gdx * d.x * dL_dG;
			dL_dconic2D_shared[tid].y = skip ? 0.f : -0.5f * gdx * d.y * dL_dG;
			dL_dconic2D_shared[tid].w = skip ? 0.f : -0.5f * gdy * d.y * dL_dG;
			dL_dopacity_shared[tid] = skip ? 0.f : G * dL_dopa;

			render_cuda_reduce_sum(block, 
				dL_dmean2D_shared,
				dL_dconic2D_shared,
				dL_dopacity_shared,
				dL_dcolors_shared, 
				dL_dts_shared,
				dL_dray_planes_shared
			);

			if (tid == 0) {
				float2 dL_dmean2D_acc = dL_dmean2D_shared[0];
				float4 dL_dconic2D_acc = dL_dconic2D_shared[0];
				float dL_dopacity_acc = dL_dopacity_shared[0];
				float3 dL_dcolors_acc = dL_dcolors_shared[0];
				float dL_dts_acc = dL_dts_shared[0];
				float2 dL_dray_planes_acc = dL_dray_planes_shared[0];

				atomicAdd(&dL_dmean2D[global_id].x, dL_dmean2D_acc.x);
				atomicAdd(&dL_dmean2D[global_id].y, dL_dmean2D_acc.y);
				atomicAdd(&dL_dconic2D[global_id].x, dL_dconic2D_acc.x);
				atomicAdd(&dL_dconic2D[global_id].y, dL_dconic2D_acc.y);
				atomicAdd(&dL_dconic2D[global_id].w, dL_dconic2D_acc.w);
				atomicAdd(&dL_dopacity[global_id], dL_dopacity_acc);
				atomicAdd(&dL_dcolors[global_id * C + 0], dL_dcolors_acc.x);
				atomicAdd(&dL_dcolors[global_id * C + 1], dL_dcolors_acc.y);
				atomicAdd(&dL_dcolors[global_id * C + 2], dL_dcolors_acc.z);
				atomicAdd(&dL_dts[global_id], dL_dts_acc);
				atomicAdd(&dL_dray_planes[global_id].x, dL_dray_planes_acc.x);
				atomicAdd(&dL_dray_planes[global_id].y, dL_dray_planes_acc.y);
			}
		}
	}
}

void BACKWARD::preprocess(
	int P, int D, int M,
	const float3* means3D,
	const int* radii,
	const float* shs,
	const bool* clamped,
	const glm::vec3* scales,
	const glm::vec4* rotations,
	const float scale_modifier,
	const float* cov3Ds,
	const float* viewmatrix,
	const float* projmatrix,
	const float* projmatrix_raw,
	const float focal_x, float focal_y,
	const float tan_fovx, float tan_fovy,
	const glm::vec3* campos,
	const float4* dL_dmean2D,
	const float* dL_dconic,
	glm::vec3* dL_dmean3D,
	float* dL_dcolor,
	const float* dL_dts,
	const float2* dL_dray_plane,
	float* dL_dcov3D,
	float* dL_dsh,
	glm::vec3* dL_dscale,
	glm::vec4* dL_drot,
	const float4* conic_opacity,
	float* dL_dopacity,
	float* dL_dtau)
{

	computeCov2DCUDA << <(P + 255) / 256, 256 >> > (
		P,
		means3D,
		radii,
		cov3Ds,
		focal_x,
		focal_y,
		tan_fovx,
		tan_fovy,
		viewmatrix,
		dL_dconic,
		dL_dray_plane,
		(glm::vec3*)dL_dmean3D,
		dL_dcov3D,
		conic_opacity,
		dL_dopacity,
		dL_dtau);

	preprocessCUDA<NUM_CHANNELS> << < (P + 255) / 256, 256 >> > (
		P, D, M,
		(float3*)means3D,
		radii,
		shs,
		clamped,
		(glm::vec3*)scales,
		(glm::vec4*)rotations,
		scale_modifier,
		viewmatrix,
		projmatrix,
		projmatrix_raw,
		campos,
		(float4*)dL_dmean2D,
		(glm::vec3*)dL_dmean3D,
		dL_dcolor,
		dL_dts,
		dL_dcov3D,
		dL_dsh,
		dL_dscale,
		dL_drot,
		dL_dtau);
}

void BACKWARD::render(
	const dim3 grid, const dim3 block,
	const uint2* ranges,
	const uint32_t* point_list,
	int W, int H,
	const float* bg_color,
	const float* view_points,
	const float2* means2D,
	const float4* conic_opacity,
	const float* colors,
	const float* depths,
	const float* ts,
	const float2* ray_planes,
	const float* alphas,
	const float* accum_depth,
	const uint32_t* n_contrib,
	const float* dL_dpixels,
	const float* dL_dpixel_depth,
	const float focal_x, 
	const float focal_y,
	float4* dL_dmean2D,
	float4* dL_dconic2D,
	float* dL_dopacity,
	float* dL_dcolors,
	float* dL_dts,
	float2* dL_dray_planes)
{
    renderCUDA<NUM_CHANNELS, true> <<<grid, block>>> ( \
        ranges, point_list, W, H, bg_color, view_points, means2D, conic_opacity, colors, \
        depths, ts, ray_planes, alphas, \
		accum_depth, \
        n_contrib, dL_dpixels, dL_dpixel_depth, \
        focal_x, focal_y, dL_dmean2D, dL_dconic2D, dL_dopacity, dL_dcolors, \
        dL_dts, dL_dray_planes);
}