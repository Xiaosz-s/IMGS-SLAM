#ifndef EIGEN_REALSVD2X2_H
#define EIGEN_REALSVD2X2_H

namespace Eigen {

namespace internal {

template<typename MatrixType, typename RealScalar, typename Index>
void real_2x2_jacobi_svd(const MatrixType& matrix, Index p, Index q,
                         JacobiRotation<RealScalar> *j_left,
                         JacobiRotation<RealScalar> *j_right)
{
  using std::sqrt;
  using std::abs;
  Matrix<RealScalar,2,2> m;
  m << numext::real(matrix.coeff(p,p)), numext::real(matrix.coeff(p,q)),
       numext::real(matrix.coeff(q,p)), numext::real(matrix.coeff(q,q));
  JacobiRotation<RealScalar> rot1;
  RealScalar t = m.coeff(0,0) + m.coeff(1,1);
  RealScalar d = m.coeff(1,0) - m.coeff(0,1);

  if(abs(d) < (std::numeric_limits<RealScalar>::min)())
  {
    rot1.s() = RealScalar(0);
    rot1.c() = RealScalar(1);
  }
  else
  {

    RealScalar u = t / d;
    RealScalar tmp = sqrt(RealScalar(1) + numext::abs2(u));
    rot1.s() = RealScalar(1) / tmp;
    rot1.c() = u / tmp;
  }
  m.applyOnTheLeft(0,1,rot1);
  j_right->makeJacobi(m,0,1);
  *j_left = rot1 * j_right->transpose();
}

} 

} 

#endif 
