#ifndef EIGEN_GENERALIZEDSELFADJOINTEIGENSOLVER_H
#define EIGEN_GENERALIZEDSELFADJOINTEIGENSOLVER_H

#include "./Tridiagonalization.h"

namespace Eigen { 

template<typename MatrixType_>
class GeneralizedSelfAdjointEigenSolver : public SelfAdjointEigenSolver<MatrixType_>
{
    typedef SelfAdjointEigenSolver<MatrixType_> Base;
  public:

    typedef MatrixType_ MatrixType;

    GeneralizedSelfAdjointEigenSolver() : Base() {}

    explicit GeneralizedSelfAdjointEigenSolver(Index size)
        : Base(size)
    {}

    GeneralizedSelfAdjointEigenSolver(const MatrixType& matA, const MatrixType& matB,
                                      int options = ComputeEigenvectors|Ax_lBx)
      : Base(matA.cols())
    {
      compute(matA, matB, options);
    }

    GeneralizedSelfAdjointEigenSolver& compute(const MatrixType& matA, const MatrixType& matB,
                                               int options = ComputeEigenvectors|Ax_lBx);

  protected:

};

template<typename MatrixType>
GeneralizedSelfAdjointEigenSolver<MatrixType>& GeneralizedSelfAdjointEigenSolver<MatrixType>::
compute(const MatrixType& matA, const MatrixType& matB, int options)
{
  eigen_assert(matA.cols()==matA.rows() && matB.rows()==matA.rows() && matB.cols()==matB.rows());
  eigen_assert((options&~(EigVecMask|GenEigMask))==0
          && (options&EigVecMask)!=EigVecMask
          && ((options&GenEigMask)==0 || (options&GenEigMask)==Ax_lBx
           || (options&GenEigMask)==ABx_lx || (options&GenEigMask)==BAx_lx)
          && "invalid option parameter");

  bool computeEigVecs = ((options&EigVecMask)==0) || ((options&EigVecMask)==ComputeEigenvectors);

  LLT<MatrixType> cholB(matB);

  int type = (options&GenEigMask);
  if(type==0)
    type = Ax_lBx;

  if(type==Ax_lBx)
  {
    
    MatrixType matC = matA.template selfadjointView<Lower>();
    cholB.matrixL().template solveInPlace<OnTheLeft>(matC);
    cholB.matrixU().template solveInPlace<OnTheRight>(matC);

    Base::compute(matC, computeEigVecs ? ComputeEigenvectors : EigenvaluesOnly );

    if(computeEigVecs)
      cholB.matrixU().solveInPlace(Base::m_eivec);
  }
  else if(type==ABx_lx)
  {
    
    MatrixType matC = matA.template selfadjointView<Lower>();
    matC = matC * cholB.matrixL();
    matC = cholB.matrixU() * matC;

    Base::compute(matC, computeEigVecs ? ComputeEigenvectors : EigenvaluesOnly);

    if(computeEigVecs)
      cholB.matrixU().solveInPlace(Base::m_eivec);
  }
  else if(type==BAx_lx)
  {
    
    MatrixType matC = matA.template selfadjointView<Lower>();
    matC = matC * cholB.matrixL();
    matC = cholB.matrixU() * matC;

    Base::compute(matC, computeEigVecs ? ComputeEigenvectors : EigenvaluesOnly);

    if(computeEigVecs)
      Base::m_eivec = cholB.matrixL() * Base::m_eivec;
  }

  return *this;
}

} 

#endif 
