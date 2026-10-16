#ifndef EIGEN_HESSENBERGDECOMPOSITION_H
#define EIGEN_HESSENBERGDECOMPOSITION_H

namespace Eigen { 

namespace internal {
  
template<typename MatrixType> struct HessenbergDecompositionMatrixHReturnType;
template<typename MatrixType>
struct traits<HessenbergDecompositionMatrixHReturnType<MatrixType> >
{
  typedef MatrixType ReturnType;
};

}

template<typename MatrixType_> class HessenbergDecomposition
{
  public:

    typedef MatrixType_ MatrixType;

    enum {
      Size = MatrixType::RowsAtCompileTime,
      SizeMinusOne = Size == Dynamic ? Dynamic : Size - 1,
      Options = MatrixType::Options,
      MaxSize = MatrixType::MaxRowsAtCompileTime,
      MaxSizeMinusOne = MaxSize == Dynamic ? Dynamic : MaxSize - 1
    };

    typedef typename MatrixType::Scalar Scalar;
    typedef Eigen::Index Index; 

    typedef Matrix<Scalar, SizeMinusOne, 1, Options & ~RowMajor, MaxSizeMinusOne, 1> CoeffVectorType;

    typedef HouseholderSequence<MatrixType,typename internal::remove_all<typename CoeffVectorType::ConjugateReturnType>::type> HouseholderSequenceType;
    
    typedef internal::HessenbergDecompositionMatrixHReturnType<MatrixType> MatrixHReturnType;

    explicit HessenbergDecomposition(Index size = Size==Dynamic ? 2 : Size)
      : m_matrix(size,size),
        m_temp(size),
        m_isInitialized(false)
    {
      if(size>1)
        m_hCoeffs.resize(size-1);
    }

    template<typename InputType>
    explicit HessenbergDecomposition(const EigenBase<InputType>& matrix)
      : m_matrix(matrix.derived()),
        m_temp(matrix.rows()),
        m_isInitialized(false)
    {
      if(matrix.rows()<2)
      {
        m_isInitialized = true;
        return;
      }
      m_hCoeffs.resize(matrix.rows()-1,1);
      _compute(m_matrix, m_hCoeffs, m_temp);
      m_isInitialized = true;
    }

    template<typename InputType>
    HessenbergDecomposition& compute(const EigenBase<InputType>& matrix)
    {
      m_matrix = matrix.derived();
      if(matrix.rows()<2)
      {
        m_isInitialized = true;
        return *this;
      }
      m_hCoeffs.resize(matrix.rows()-1,1);
      _compute(m_matrix, m_hCoeffs, m_temp);
      m_isInitialized = true;
      return *this;
    }

    const CoeffVectorType& householderCoefficients() const
    {
      eigen_assert(m_isInitialized && "HessenbergDecomposition is not initialized.");
      return m_hCoeffs;
    }

    const MatrixType& packedMatrix() const
    {
      eigen_assert(m_isInitialized && "HessenbergDecomposition is not initialized.");
      return m_matrix;
    }

    HouseholderSequenceType matrixQ() const
    {
      eigen_assert(m_isInitialized && "HessenbergDecomposition is not initialized.");
      return HouseholderSequenceType(m_matrix, m_hCoeffs.conjugate())
             .setLength(m_matrix.rows() - 1)
             .setShift(1);
    }

    MatrixHReturnType matrixH() const
    {
      eigen_assert(m_isInitialized && "HessenbergDecomposition is not initialized.");
      return MatrixHReturnType(*this);
    }

  private:

    typedef Matrix<Scalar, 1, Size, int(Options) | int(RowMajor), 1, MaxSize> VectorType;
    typedef typename NumTraits<Scalar>::Real RealScalar;
    static void _compute(MatrixType& matA, CoeffVectorType& hCoeffs, VectorType& temp);

  protected:
    MatrixType m_matrix;
    CoeffVectorType m_hCoeffs;
    VectorType m_temp;
    bool m_isInitialized;
};

template<typename MatrixType>
void HessenbergDecomposition<MatrixType>::_compute(MatrixType& matA, CoeffVectorType& hCoeffs, VectorType& temp)
{
  eigen_assert(matA.rows()==matA.cols());
  Index n = matA.rows();
  temp.resize(n);
  for (Index i = 0; i<n-1; ++i)
  {
    
    Index remainingSize = n-i-1;
    RealScalar beta;
    Scalar h;
    matA.col(i).tail(remainingSize).makeHouseholderInPlace(h, beta);
    matA.col(i).coeffRef(i+1) = beta;
    hCoeffs.coeffRef(i) = h;

    matA.bottomRightCorner(remainingSize, remainingSize)
        .applyHouseholderOnTheLeft(matA.col(i).tail(remainingSize-1), h, &temp.coeffRef(0));

    matA.rightCols(remainingSize)
        .applyHouseholderOnTheRight(matA.col(i).tail(remainingSize-1), numext::conj(h), &temp.coeffRef(0));
  }
}

namespace internal {

template<typename MatrixType> struct HessenbergDecompositionMatrixHReturnType
: public ReturnByValue<HessenbergDecompositionMatrixHReturnType<MatrixType> >
{
  public:

    HessenbergDecompositionMatrixHReturnType(const HessenbergDecomposition<MatrixType>& hess) : m_hess(hess) { }

    template <typename ResultType>
    inline void evalTo(ResultType& result) const
    {
      result = m_hess.packedMatrix();
      Index n = result.rows();
      if (n>2)
        result.bottomLeftCorner(n-2, n-2).template triangularView<Lower>().setZero();
    }

    Index rows() const { return m_hess.packedMatrix().rows(); }
    Index cols() const { return m_hess.packedMatrix().cols(); }

  protected:
    const HessenbergDecomposition<MatrixType>& m_hess;
};

} 

} 

#endif 
