// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_FUNCTIONS_AFFINE_HH
#define DUNE_STUFF_FUNCTIONS_AFFINE_HH

#include <memory>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/la/container/pattern.hh>
#include <dune/stuff/functions/constant.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Functions {

/**
 * \brief Simple affine function of the form f(x) = A*x + b
 */
template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class Affine
{
  Affine() { static_assert(AlwaysFalse< EntityImp >::value, "Not available for rangeDimCols > 1!"); }
};


template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim >
class Affine< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
  typedef GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > BaseType;
  typedef Affine< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >                  ThisType;

public:
  typedef typename BaseType::DomainType                                     DomainType;
  typedef typename BaseType::RangeType                                      RangeType;
  typedef typename BaseType::JacobianRangeType                              JacobianRangeType;
  typedef typename Dune::FieldMatrix< RangeFieldImp, rangeDim, domainDim >  MatrixType;
  typedef typename DS::LA::SparsityPatternDefault                           PatternType;

  using typename BaseType::LocalfunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".affine";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["A"] = internal::Get< RangeFieldImp, rangeDim, domainDim >::value_str();
    config["b"] = internal::Get< RangeFieldImp, rangeDim, 1 >::value_str();
    config["sparse"] = "false";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr< ThisType > create(const Common::Configuration config = default_config(),
                                            const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    return Common::make_unique< ThisType >(
          cfg.get("A", default_cfg.get< MatrixType >("A")),
          cfg.get("b", default_cfg.get< RangeType >("b")),
          cfg.get("sparse", default_cfg.get< bool >("sparse")),
          cfg.get("name",  default_cfg.get< std::string >("name")));
  } // ... create(...)

  explicit Affine(const MatrixType matrix, const RangeType vector = RangeType(0), const bool sparse = false, const std::string name_in = static_id())
    : A_(matrix)
    , b_(vector)
    , name_(name_in)
    , b_zero_(check_zero(b_))
    , sparse_(sparse)
    , pattern_(rangeDim)
  {
    if (sparse_)
      calculate_pattern(A_, pattern_);
  }

  Affine(const ThisType& other) = default;

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".affine";
  }

  virtual size_t order() const override final
  {
    return 1;
  }

  virtual void evaluate(const DomainType& x, RangeType& ret) const override final
  {
    if (sparse_) {
      std::fill(ret.begin(), ret.end(), RangeFieldImp(0));
      for (size_t ii = 0; ii < rangeDim; ++ii) {
        const auto& row_pattern = pattern_.inner(ii);
        for (const auto& jj : row_pattern) {
          ret[ii] += A_[ii][jj]*x[jj];
        }
      }
    } else {
      A_.mv(x, ret);
    }
    if (!b_zero_)
      ret += b_;
  }

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& ret) const override final
  {
    ret = A_;
  }

  virtual std::string name() const override final
  {
    return name_;
  }

  static bool check_zero(const RangeType& b)
  {
    for (const auto& entry : b) {
      if (DSC::FloatCmp::ne(entry, RangeFieldImp(0)))
        return false;
    }
    return true;
  }

  static void calculate_pattern(const MatrixType& A, PatternType& pattern)
  {
    for (size_t ii = 0; ii < rangeDim; ++ii) {
      const auto& row = A[ii];
      for (size_t jj = 0; jj < domainDim; ++jj) {
        if (DSC::FloatCmp::ne(row[jj], RangeFieldImp(0)))
          pattern.insert(ii,jj);
      }
    }
  }

private:
  const MatrixType A_;
  const RangeType b_;
  const bool b_zero_;
  const bool sparse_;
  PatternType pattern_;
  const std::string name_;
};


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_AFFINE_HH
