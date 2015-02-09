// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Sven Kaulmann

#ifndef DUNESTUFF_PRINTING_HH_INCLUDED
#define DUNESTUFF_PRINTING_HH_INCLUDED

#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <vector>

#include <boost/format.hpp>

#include <dune/common/deprecated.hh>
#include <dune/common/densematrix.hh>
#include <dune/common/densevector.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/string.hh>

namespace Dune {
namespace Stuff {
namespace Common {

//! ensure matlab output is done with highest precision possible, otherwise weird effects are bound to happen
static const auto matlab_output_precision = std::numeric_limits< double >::digits10 + 1;

template< class OutStreamType = std::ostream >
void print(const double& d,
           const std::string name = "double",
           OutStreamType& out = std::cout,
           const std::string prefix = "")
{
  out << prefix << name << " = " << d << "\n";
} // void print(const double& d, ...)

template< class T, class OutStreamType = std::ostream >
void print(const std::vector< T >& vector,
           const std::string name = "vector",
           OutStreamType& out = std::cout,
           const std::string prefix = "")
{
  if (vector.size() == 1)
    print(vector[0], name, out, prefix);
  else {
    out << prefix << name << " = [";
    for (auto i : valueRange(vector.size() - 1)) {
      out << vector[i] << ", ";
    }
    out << vector[vector.size() - 1] << "];" << "\n";
  }
} // void print(const std::vector< double >& ds, ...)

template< class VectorImp,
          class OutStreamType = std::ostream >
void print(const Dune::DenseVector< VectorImp >& vector,
           const std::string name = "DenseVector",
           OutStreamType& out = std::cout,
           const std::string prefix = "")
{
  if (vector.size() == 1)
    print(vector[0], name, out, prefix);
  else {
    out << prefix << name << " = [";
    for (auto i : valueRange(vector.size() - 1)) {
      out << vector[i] << ", ";
    }
    out << vector[vector.size() - 1] << "];" << "\n";
  }
} // void print(const Dune::DenseVector< VectorImp >& vector, ...)

template< class FieldImp, size_t size,
          class OutStreamType = std::ostream >
void print(const std::vector< Dune::FieldVector< FieldImp, size > >& vectors,
           const std::string name = "vector_of_FieldVector",
           OutStreamType& out = std::cout,
           const std::string prefix = "")
{
  if (vectors.size() == 1)
    print(vectors[0], name, out, prefix);
  else {
    for (auto i : valueRange(vectors.size())) {
      print(vectors[i], name + "[" + Dune::Stuff::Common::toString(i) + "]", out, prefix);
    }
  }
} // void print(const std::vector< Dune::FieldVector< FieldImp, size > >& vectors, ...)

/**
 *  \attention  I'm not sure why we need this, print< const std::vector< Dune::DenseVector< VectorImp > >& >
 *              does not seem to work for Dune::FieldVector.
 */
template< class VectorImp,
          class OutStreamType = std::ostream >
void print(const std::vector< Dune::DenseVector< VectorImp > >& vectors,
           const std::string name = "vector_of_DenseVector",
           OutStreamType& out = std::cout,
           const std::string prefix = "")
{
  if (vectors.size() == 1)
    print(vectors[0], name, out, prefix);
  else {
    for (auto i : valueRange(vectors.size())) {
      print(vectors[i], name + "[" + Dune::Stuff::Common::toString(i) + "]", out, prefix);
    }
  }
} // void print(const std::vector< const Dune::DenseVector< VectorImp > >& vectors, ...)

template< class MatrixImp,
          class OutStreamType = std::ostream >
void print(const Dune::DenseMatrix< MatrixImp >& matrix,
           const std::string name = "DenseMatrix",
           OutStreamType& out = std::cout,
           const std::string prefix = "")
{
  if (matrix.rows() == 1 && matrix.cols() == 1)
    print(matrix[0][0], name, out, prefix);
  else {
    out << prefix << name << " = [";
    typedef typename Dune::DenseMatrix< MatrixImp >::const_row_reference RowType;
    const RowType& firstRow = matrix[0];
    for (auto j : valueRange((firstRow.size() - 1))) {
      out << firstRow[j] << ", ";
    }
    out << firstRow[firstRow.size() - 1];
    if (matrix.rows() == 1)
      out << "];" << "\n";
    else
      out << ";" << "\n";
    for (size_t i = 1; i < matrix.rows(); ++i) {
      out << prefix << whitespaceify(name + " = [");
      const RowType& row = matrix[i];
      for (auto j : valueRange(row.size() - 1)) {
        out << row[j] << ", ";
      }
      out << row[row.size() - 1];
      if (i == matrix.rows() - 1)
        out << "];";
      out << "\n";
    }
  }
} // void print(const Dune::DenseMatrix< MatrixImp >& matrix, ...)

/**
 *  \attention  I'm not sure why we need this, print< const std::vector< Dune::DenseMatrix< MatrixImp > >& >
 *              does not seem to work for Dune::FieldMatrix.
 */
template< class FieldImp, size_t rows, size_t cols, class OutStreamType = std::ostream >
void print(const std::vector< Dune::FieldMatrix< FieldImp, rows, cols > >& matrices,
           const std::string name = "vector_of_FieldMatrix",
           OutStreamType& out = std::cout,
           const std::string prefix = "")
{
  if (matrices.size() == 1)
    print(matrices[0], name, out, prefix);
  else {
    for (auto i : valueRange(matrices.size())) {
      print(matrices[i], name + "[" + Dune::Stuff::Common::toString(i) + "]", out, prefix);
    }
  }
} // void print(const std::vector< Dune::FieldMatrix< FieldImp, rows, cols > >& matrices, ...)

template< class MatrixImp,
          class OutStreamType = std::ostream >
void print(const std::vector< Dune::DenseMatrix< MatrixImp > >& matrices,
           const std::string name = "vector_of_DenseMatrix",
           OutStreamType& out = std::cout,
           const std::string prefix = "")
{
  if (matrices.size() == 1)
    print(matrices[0], name, out, prefix);
  else {
    for (auto i : valueRange(matrices.size())) {
      print(matrices[i], name + "[" + Dune::Stuff::Common::toString(i) + "]", out, prefix);
    }
  }
} // void print(const std::vector< Dune::DenseMatrix< MatrixImp > >& matrices, ...)

/**
   *  \brief prints a Dune::FieldVector
   *
   *  or anything compatible in terms of Iterators
   *  \tparam T
   *          should be Dune::FieldVector or compatible
   *  \tparam stream
   *          std::ostream or compatible
   *  \param[in]  arg
   *          vector to be printed
   *  \param[in]  name
   *          name to be printed along
   *  \param[in]  out
   *          where to print
   *  \param[opt] prefix
   *          prefix to be printed before every line
   **/
template< class T, class stream >
void
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
     printFieldVector(T& arg, std::string name, stream& out, std::string prefix = "") {
  out << "\n" << prefix << "printing " << name << " (Dune::FieldVector)" << "\n";
  out << prefix;
  for (auto value : arg)
  {
    out << std::setw(14) << std::setprecision(6) << value;
  }
  out << '\n';
} // printFieldVector

/**
   *  \brief prints a Dune::FieldMatrix
   *
   *  or anything compatible in terms of Iterators
   *  \tparam T
   *          should be Dune::FieldVector or compatible
   *  \tparam stream
   *          std::ostream or compatible
   *  \param[in]  arg
   *          matrix to be printed
   *  \param[in]  name
   *          name to be printed along
   *  \param[in]  out
   *          where to print
   *  \param[opt] prefix
   *          prefix to be printed before every line
   **/
template< class T, class stream >
void
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
     printFieldMatrix(T& arg, std::string name, stream& out, std::string prefix = "") {
  out << "\n" << prefix << "printing " << name << " (Dune::FieldMatrix)";
  size_t row = 1;
  for (auto rIt : arg) {
    out << "\n" << prefix << "  row " << row << ":";
    for (auto vIt : rIt)
      out << std::setw(14) << std::setprecision(6) << vIt;
    row += 1;
  }
} // printFieldMatrix

/** \brief print a SparseRowMatrix (or any interface conforming object) to a given stream in matlab (laodable-) format
   * \ingroup Matlab
   **/
template< class T, class stream >
void
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
     printSparseRowMatrixMatlabStyle( const T& arg, std::string name, stream& out,
                                      const double eps = Config().get("eps", 1e-14) ) {
  name = std::string("fem.") + name;
  out << boost::format("\n%s =sparse( %d, %d );") % name % arg.rows() % arg.cols() << "\n";
  for (auto row : valueRange(arg.rows()))
  {
    for (auto col : valueRange(arg.cols()))
    {
      const typename T::Ttype val = arg(row, col);
      if (std::fabs(val) > eps)
        out << name << "(" << row + 1 << "," << col + 1 << ")=" << std::setprecision(matlab_output_precision) << val
            << ";\n";
    }
  }
} // printSparseRowMatrixMatlabStyle

/** \brief print a ISTLMatrix (or any interface conforming object) to a given stream in matlab (laodable-) format
   * \ingroup Matlab
   **/
template< class MatrixType, class stream >
void
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
     printISTLMatrixMatlabStyle( const MatrixType& arg, std::string name, stream& out,
                                 const double eps = Config().get("eps", 1e-14) ) {
  name = std::string("istl.") + name;
  const auto I = arg.N();
  const auto J = arg.M();
  typedef typename MatrixType::block_type BlockType;
  out << boost::format("\n%s =sparse( %d, %d );") % name % (I* BlockType::rows) % (J* BlockType::cols) << "\n";
  for (auto ii : valueRange(I))
  {
    for (auto jj : valueRange(J))
    {
      if ( arg.exists(ii, jj) )
      {
        const auto& block = arg[ii][jj];
        for (auto i : valueRange(block.N()))
        {
          for (auto j : valueRange(block.M()))
          {
            const auto value = block[i][j];
            if (std::fabs(value) > eps)
            {
              size_t real_row = BlockType::rows * ii + i + 1;
              size_t real_col = BlockType::cols * jj + j + 1;
              out << name << "(" << real_row << "," << real_col << ")="
                  << std::setprecision(matlab_output_precision) << value << ";\n";
            }
          }
        }
      }
    }
  }
} // printISTLMatrixMatlabStyle

/** \brief print a discrete function (or any interface conforming object) to a given stream in matlab (laodable-) format
   * \ingroup Matlab
   **/
template< class T, class stream >
void
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
     printDiscreteFunctionMatlabStyle(const T& arg, const std::string name, stream& out) {
  out << "\n" << name << " = [ " << "\n";
  for (auto val : arg)
  {
    out << std::setprecision(matlab_output_precision) << val;
    out << ";" << "\n";
  }
  out << "];" << "\n";
} // printDiscreteFunctionMatlabStyle

/** \brief print a double vector (or any interface conforming object) to a given stream in matlab (laodable-) format
   * \ingroup Matlab
   **/
template< class T, class stream >
void
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
     printDoubleVectorMatlabStyle(const T* arg, const size_t size, const std::string name, stream& out) {
  out << "\n" << name << " = [ " << "\n";
  for (auto i : valueRange(size))
  {
    out << std::setprecision(matlab_output_precision) << arg[i];
    out << ";" << "\n";
  }
  out << "];" << "\n";
} // printDoubleVectorMatlabStyle

//! simple vector to stream print
template< class Type >
void
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
     printDoubleVec(std::ostream& stream, const Type* vec, const size_t N) {
  stream << "\n [ " << std::setw(5);
  for (auto i : valueRange(N))
    stream << vec[i] << " ";

  stream << " ] " << "\n";
} // printDoubleVec

//! simple discrete function to stream print
template< class DiscFunc >
void
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
     oneLinePrint(std::ostream& stream, const DiscFunc& func) {
  typedef typename DiscFunc::ConstDofIteratorType DofIteratorType;
  DofIteratorType it = func.dbegin();
  stream << "\n" << func.name() << ": \n[ ";
  for ( ; it != func.dend(); ++it)
  {
    // double d = 0.10;// + *it; //stupid hack cause setw/prec ain't working for me
    stream << std::setw(6) << std::setprecision(3) << *it << "  ";
  }
  stream << " ] " << "\n";
} // oneLinePrint

/** \brief localmatrix printing functor for use in Stuff::GridWalk
   * putting this into Stuff::GridWalk::operator() will result in a local matrix being printed for each gird entity\n
   * Example:\n
   * Stuff::GridWalk<GridPartType> gw( gridPart_ );\n
   * Stuff::LocalMatrixPrintFunctor< RmatrixType,FunctorStream> f_R ( Rmatrix, functorStream, "R" );\n
   * gw( f_R );
   * \see Stuff::GridWalk
   * \ingroup GridWalk
   **/
template< class GlobalMatrix>
class
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
      LocalMatrixPrintFunctor
{
public:
  LocalMatrixPrintFunctor(const GlobalMatrix& m, std::ostream& stream, const std::string name)
    : matrix_(m)
      , stream_(stream)
      , name_(name)
  {}

  template< class Entity >
  void operator()(const Entity& en, const Entity& ne, const size_t en_idx, const size_t ne_idx) {
    typename GlobalMatrix::LocalMatrixType localMatrix
      = matrix_.localMatrix(en, ne);
    stream_ << "\nlocal_" << name_ << "_Matrix_" << en_idx << "_" << ne_idx << " = [" << "\n";
    const auto rows = localMatrix.rows();
    const auto cols = localMatrix.columns();
    for (auto i : valueRange(rows))
    {
      for (auto j : valueRange(cols))
      {
        stream_ << std::setw(8) << std::setprecision(2) << localMatrix.get(i, j);
      }
      stream_ << ";" << "\n";
    }
    stream_ << "];" << "\n";
  } // ()

  void preWalk() {
    stream_ << "% printing local matrizes of " << name_ << "\n";
  }

  void postWalk() {
    stream_ << "\n% done printing local matrizes of " << name_ << "\n";
  }

private:
  const GlobalMatrix& matrix_;
  std::ostream& stream_;
  const std::string name_;
};

/** GridWalk functor to print all localfunctions of a given DiscreteFunction
   * \ingroup GridWalk
   **/
template< class DiscreteFunctionType, class QuadratureType >
class
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
      LocalFunctionPrintFunctor
{
public:
  LocalFunctionPrintFunctor(const DiscreteFunctionType& discrete_function, std::ostream& stream)
    : discrete_function_(discrete_function)
      , stream_(stream)
      , name_( discrete_function.name() )
  {}

  template< class Entity >
  void operator()(const Entity& en, const Entity& /*ne*/, const size_t /*en_idx*/, const size_t /*ne_idx */) {
    typename DiscreteFunctionType::LocalFunctionType lf = discrete_function_.localFunction(en);
    QuadratureType quad(en, 2 * discrete_function_.space().order() + 2);
    for (auto qp : valueRange(quad.nop()))
    {
      typename DiscreteFunctionType::RangeType eval(0);
      typename DiscreteFunctionType::DomainType xLocal = quad.point(qp);
      typename DiscreteFunctionType::DomainType xWorld = en.geometry().global(xLocal);
      lf.evaluate(xLocal, eval);
      stream_ << boost::format("xWorld %f \t %s value %f\n") % xWorld % name_ % eval;
    }
  } // ()

  void preWalk() {
    stream_ << "% printing local function values of " << name_ << "\n";
  }

  void postWalk() {
    stream_ << "\n% done printing function values of " << name_ << "\n";
  }

private:
  const DiscreteFunctionType& discrete_function_;
  std::ostream& stream_;
  const std::string name_;
};

/** GridWalk functor to print, w/o transformation, all localfunctions of a given DiscreteFunction
   * \ingroup GridWalk
   **/
template< class DiscreteFunctionType >
class
  DUNE_DEPRECATED_MSG("Will be removed soon, file an issue on https://github.com/wwu-numerik/dune-stuff/issues if you need this (09.02.2015)!")
      LocalFunctionVerbatimPrintFunctor
{
public:
  LocalFunctionVerbatimPrintFunctor(const DiscreteFunctionType& discrete_function, std::ostream& stream)
    : discrete_function_(discrete_function)
      , stream_(stream)
      , name_( discrete_function.name() )
  {}

  template< class Entity >
  void operator()(const Entity& en, const Entity& /*ne*/, const size_t /*en_idx*/, const size_t /*ne_idx */) {
    typename DiscreteFunctionType::LocalFunctionType lf = discrete_function_.localFunction(en);

    for (auto qp : valueRange(lf.numDofs()))
    {
      stream_ << boost::format("%s dof %d value value %f\n") % name_ % qp % lf[qp];
    }
  } // ()

  void preWalk() {
    stream_ << "% printing local function values of " << name_ << "\n";
  }

  void postWalk() {
    stream_ << "\n% done printing function values of " << name_ << "\n";
  }

private:
  const DiscreteFunctionType& discrete_function_;
  std::ostream& stream_;
  const std::string name_;
};

//! useful for visualizing sparsity patterns of matrices
template< class Matrix >
void matrixToGnuplotStream(const Matrix& matrix, std::ostream& stream) {
  unsigned long nz = 0;

  const auto cols = matrix.cols();
  for (auto row : valueRange(matrix.rows()))
  {
    for (auto col : valueRange(cols))
    {
      if ( matrix.find(row, col) )
        stream << row << "\t" << col << "\t" << matrix(row, col) << "\n";
    }
    nz += matrix.numNonZeros(row);
    stream << "#non zeros in row " << row << " " << matrix.numNonZeros(row) << " (of " << matrix.cols() << " cols)\n";
  }
  stream << "#total non zeros " << nz << " of " << matrix.rows() * matrix.cols() << " entries\n";
} // matrixToGnuplotStream

//! proxy to Stuff::matrixToGnuplotStream that redirects its output to a file
template< class Matrix >
void matrixToGnuplotFile(const Matrix& matrix, std::string filename) {
  std::string dir(Config().get( "fem.io.datadir", std::string("data") ) + "/gnuplot/");

  testCreateDirectory(dir);
  std::ofstream file( (dir + filename).c_str() );
  matrixToGnuplotStream(matrix, file);
  file.flush();
  file.close();
} // matrixToGnuplotFile

//! maps 1,2,3 to x,y,z / X,Y,Z
inline std::string dimToAxisName(const size_t dim, const bool capitalize = false) {
  char c = 'x';

  c += dim;
  if (capitalize)
    c -= 32;
  return std::string() += c;
} // dimToAxisName

} // namespace Common
} // namespace Stuff
} // namespace Dune

#endif // ifndef DUNESTUFF_PRINTING_HH_INCLUDED
