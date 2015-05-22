// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
#define DUNE_STUFF_FUNCTION_CHECKERBOARD_HH

#include <vector>
#include <cmath>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/debug.hh>
#include <dune/stuff/common/fvector.hh>

#include "interfaces.hh"
#include "expression.hh"

namespace Dune {
namespace Stuff {
namespace Functions {


template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class Checkerboard
  : public LocalizableFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
{
  typedef LocalizableFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
    BaseType;
  typedef Checkerboard< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols > ThisType;

  class Localfunction
    : public LocalfunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
  {
    typedef LocalfunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols > BaseType;
  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainType        DomainType;
    typedef typename BaseType::RangeFieldType    RangeFieldType;
    typedef typename BaseType::RangeType         RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const EntityType& ent, const RangeType value)
      : BaseType(ent)
      , value_(value)
    {}

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const override
    {
      return 0;
    }

    virtual void evaluate(const DomainType& UNUSED_UNLESS_DEBUG(xx), RangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      ret = value_;
    }

    virtual void jacobian(const DomainType& UNUSED_UNLESS_DEBUG(xx), JacobianRangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      jacobian_helper(ret, internal::ChooseVariant< rangeDimCols >());
    }

  private:
    template< size_t rC >
    void jacobian_helper(JacobianRangeType& ret, internal::ChooseVariant< rC >) const
    {
      for (auto& col_jacobian: ret) {
        col_jacobian *= RangeFieldType(0);
      }
    }

    void jacobian_helper(JacobianRangeType& ret, internal::ChooseVariant< 1 >) const
    {
      ret *= RangeFieldType(0);
    }
    const RangeType value_;
  }; // class Localfunction

public:
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;

  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;

  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::RangeType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".checkerboard";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["num_elements"] = "[2 2 2]";
    config["values"] = "[1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]";
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
    // calculate number of values and get values
    auto num_elements = cfg.get("num_elements",
                                default_cfg.get< Common::FieldVector< size_t, dimDomain > >("num_elements"), dimDomain);
    size_t num_values = 1;
    for (size_t ii = 0; ii < num_elements.size(); ++ii)
      num_values *= num_elements[ii];
    std::vector< RangeType > values(num_values);
    if (config.has_key("values.0")) { // get every value from its own config entry
      try { // get value directly as RangeType
        for (size_t ii = 0; ii < num_values; ++ii)
          values[ii] = cfg.get< RangeType >("values." + DSC::toString(ii),
                                            dimRange,
                                            dimRangeCols);
      } catch (const Exceptions::conversion_error& e) {
        if (dimRangeCols == 1) { // get every value from its own config entry as the first col of a matrix
          for (size_t ii = 0; ii < num_values; ++ii) {
            const auto values_matrix = cfg.get< Common::FieldMatrix< RangeFieldType, dimRange, 1 > >("values." + DSC::toString(ii),
                                                                                                     dimRange,
                                                                                                     1);
            // this fromString(toString(...)) construct avoids a compilation error if dimRangeCols > 1 and is easier
            // than creating templated helper methods
            values[ii] = DSC::fromString< RangeType >(DSC::toString(values_matrix[ii]));
          }
        } else {
          std::cout << e.what() << std::endl;
        }
      }
    } else {
      // get values as a vector of scalars
      auto values_rf = cfg.get("values", default_cfg.get< std::vector< RangeFieldType > >("values"), num_values);
      for (size_t ii = 0; ii < values_rf.size(); ++ii)
        values[ii] = RangeType(values_rf[ii]);
    }
    // create
    return Common::make_unique< ThisType >(
            cfg.get("lower_left",
                    default_cfg.get< Common::FieldVector< DomainFieldType, dimDomain > >("lower_left"), dimDomain),
            cfg.get("upper_right",
                    default_cfg.get< Common::FieldVector< DomainFieldType, dimDomain > >("upper_right"), dimDomain),
            std::move(num_elements),
            std::move(values),
            cfg.get("name", default_cfg.get< std::string > ("name")));
  } // ... create(...)

  Checkerboard(const Common::FieldVector< DomainFieldType, dimDomain >& lowerLeft,
               const Common::FieldVector< DomainFieldType, dimDomain >& upperRight,
               const Common::FieldVector< size_t, dimDomain >& numElements,
               const std::vector< RangeType >& values,
               const std::string nm = static_id())
    : lowerLeft_(new Common::FieldVector< DomainFieldType, dimDomain >(lowerLeft))
    , upperRight_(new Common::FieldVector< DomainFieldType, dimDomain >(upperRight))
    , numElements_(new Common::FieldVector< size_t, dimDomain >(numElements))
    , values_(new std::vector< RangeType >(values))
    , name_(nm)
  {
    // checks
    size_t totalSubdomains = 1;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      const auto& ll = (*lowerLeft_)[dd];
      const auto& ur = (*upperRight_)[dd];
      const auto& ne = (*numElements_)[dd];
      if (!(ll < ur))
        DUNE_THROW(Dune::RangeError, "lowerLeft has to be elementwise smaller than upperRight!");
      totalSubdomains *= ne;
    }
    if (values_->size() < totalSubdomains)
      DUNE_THROW(Dune::RangeError,
                 "values too small (is " << values_->size() << ", should be " << totalSubdomains << ")");
  } // Checkerboard(...)

  Checkerboard(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".checkerboard";
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityType& entity) const override
  {
    // decide on the subdomain the center of the entity belongs to
    const auto center = entity.geometry().center();
    std::vector< size_t > whichPartition(dimDomain, 0);
    const auto& ll = *lowerLeft_;
    const auto& ur = *upperRight_;
    const auto& ne = *numElements_;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // for points that are on upperRight_[d], this selects one partition too much
      // so we need to cap this
      whichPartition[dd] = std::min(size_t(std::floor(ne[dd]*((center[dd] - ll[dd])/(ur[dd] - ll[dd])))),
                                    ne[dd] - 1);
    }
    size_t subdomain = 0;
    if (dimDomain == 1)
      subdomain = whichPartition[0];
    else if (dimDomain == 2)
      subdomain = whichPartition[0] + whichPartition[1]*ne[0];
    else
      subdomain = whichPartition[0] + whichPartition[1]*ne[0] + whichPartition[2]*ne[1]*ne[0];
    // return the component that belongs to the subdomain
    return std::unique_ptr< Localfunction >(new Localfunction(entity, (*values_)[subdomain]));
  } // ... local_function(...)

private:
  std::shared_ptr< const Common::FieldVector< DomainFieldType, dimDomain > > lowerLeft_;
  std::shared_ptr< const Common::FieldVector< DomainFieldType, dimDomain > > upperRight_;
  std::shared_ptr< const Common::FieldVector< size_t, dimDomain > > numElements_;
  std::shared_ptr< const std::vector< RangeType > > values_;
  std::string name_;
}; // class Checkerboard


template< class EntityImp, class DomainFieldImp, size_t domainDim,
          class RangeEntityImp, class RangeDomainFieldImp, size_t rangeDomainDim,
          class RangeRangeFieldImp, size_t rangeRangeDim, size_t rangeRangeDimCols = 1 >
class ExpressionCheckerboard
  : public GlobalFunctionValuedFunctionInterface< EntityImp, DomainFieldImp, domainDim,
                                                  RangeEntityImp, RangeDomainFieldImp, rangeDomainDim,
                                                  RangeRangeFieldImp, rangeRangeDim, rangeRangeDimCols >
{
  typedef GlobalFunctionValuedFunctionInterface< EntityImp, DomainFieldImp, domainDim,
                                                 RangeEntityImp, RangeDomainFieldImp, rangeDomainDim,
                                                 RangeRangeFieldImp, rangeRangeDim, rangeRangeDimCols >   BaseType;
  typedef ExpressionCheckerboard< EntityImp, DomainFieldImp, domainDim,
                                  RangeEntityImp, RangeDomainFieldImp, rangeDomainDim,
                                  RangeRangeFieldImp, rangeRangeDim, rangeRangeDimCols >                  ThisType;
public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeDomainType;
  using typename BaseType::RangeRangeType;
  using typename BaseType::RangeJacobianRangeType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;
  using typename BaseType::EntityType;
  using typename BaseType::RangeEntityType;
  using typename BaseType::LocalfunctionType;

  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;

  using typename BaseType::RangeDomainFieldType;
  using typename BaseType::RangeRangeFieldType;
  using BaseType::rangeDimDomain;
  using BaseType::rangeDimRange;
  using BaseType::rangeDimRangeCols;

  static const bool available = true;

  typedef Expression< RangeEntityType,
                      RangeDomainFieldType, rangeDimDomain,
                      RangeRangeFieldType, rangeDimRange, rangeDimRangeCols >         ExpressionFunctionType;
  typedef typename ExpressionFunctionType::LocalfunctionType                          ExpressionLocalfunctionType;

private:
  class Localfunction
    : public LocalfunctionType
  {
  public:
    Localfunction(const EntityType& ent, const ExpressionFunctionType& value)
      : LocalfunctionType(ent)
      , value_(std::make_shared< ExpressionFunctionType >(value))
    {}

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const override
    {
      return 10;
    }

    virtual void evaluate(const DomainType& /*xx*/, std::shared_ptr< const RangeType >& ret) const override
    {
      ret = std::shared_ptr< const RangeType >(value_);
    }

    virtual void jacobian(const DomainType& /*xx*/, std::shared_ptr< const JacobianRangeType >& /*ret*/) const override
    {
      DUNE_THROW(Dune::NotImplemented, "Not implemented");
    }

  private:
    const std::shared_ptr< const ExpressionFunctionType > value_;
  }; // class Localfunction

public:

  static std::string static_id()
  {
    return BaseType::static_id() + ".expressioncheckerboard";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["num_elements"] = "[2 2 2]";
    config["variable"] = "u";
    config["values"] = "[1.0*u[0] 2.0*u[0] 3.0*sin(u[0]) 4.0 5.0 6.0*cos(u[0]) 7.0 8.0]";
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
    // calculate number of values and get values
    auto num_elements = cfg.get("num_elements",
                                default_cfg.get< Common::FieldVector< size_t, dimDomain > >("num_elements"), dimDomain);
    size_t num_values = 1;
    for (size_t ii = 0; ii < num_elements.size(); ++ii)
      num_values *= num_elements[ii];
    const std::string variable = cfg.get("variable",   default_cfg.get< std::string >("variable"));
    std::vector< ExpressionFunctionType > values;
    if (config.has_key("values.0")) { // get every value from its own config entry
      try { // get value as matrix
        std::vector< std::string > row_as_std_vector;
        for (size_t ii = 0; ii < num_values; ++ii) {
          const auto values_matrix = cfg.get< Dune::DynamicMatrix< std::string > >("values." + DSC::toString(ii),
                                                                                   rangeDimRange,
                                                                                   rangeDimRangeCols);
          typename ExpressionFunctionType::ExpressionStringVectorType expression_vector;
          for (size_t rr = 0; rr < rangeDimRange; ++rr) {
            const auto row = values_matrix[rr];
            for (size_t cc = 0; cc < rangeDimRangeCols; ++cc)
              row_as_std_vector[cc] = row[cc];
            expression_vector.emplace_back(row_as_std_vector);
          }
         values.emplace_back(ExpressionFunctionType(variable, expression_vector));
        }
      } catch (const Exceptions::conversion_error& e) {
        if (rangeDimRangeCols == 1) { // get value as vector
          for (size_t ii = 0; ii < num_values; ++ii) {
            const auto values_vector = cfg.get< std::vector< std::string > >("values." + DSC::toString(ii),
                                                                             rangeDimRange);
            typename ExpressionFunctionType::ExpressionStringVectorType expression_vector(1, values_vector);
            values.emplace_back(ExpressionFunctionType(variable, expression_vector));
          }
        } else {
          std::cout << e.what() << std::endl;
        }
      }
    } else {
      // get values as a vector of scalars
      auto values_rf = cfg.get("values", default_cfg.get< std::vector< std::string > >("values"), num_values);
      for (size_t ii = 0; ii < values_rf.size(); ++ii)
        values.emplace_back(ExpressionFunctionType(variable, values_rf[ii]));
    }
    // create
    return Common::make_unique< ThisType >(
            cfg.get("lower_left",
                    default_cfg.get< Common::FieldVector< DomainFieldType, dimDomain > >("lower_left"), dimDomain),
            cfg.get("upper_right",
                    default_cfg.get< Common::FieldVector< DomainFieldType, dimDomain > >("upper_right"), dimDomain),
            std::move(num_elements),
            std::move(values),
            cfg.get("name", default_cfg.get< std::string > ("name")));
  } // ... create(...)

  ExpressionCheckerboard(const Common::FieldVector< DomainFieldType, dimDomain >& lowerLeft,
               const Common::FieldVector< DomainFieldType, dimDomain >& upperRight,
               const Common::FieldVector< size_t, dimDomain >& numElements,
               const std::vector< ExpressionFunctionType >& values,
               const std::string nm = static_id())
    : lowerLeft_(new Common::FieldVector< DomainFieldType, dimDomain >(lowerLeft))
    , upperRight_(new Common::FieldVector< DomainFieldType, dimDomain >(upperRight))
    , numElements_(new Common::FieldVector< size_t, dimDomain >(numElements))
    , values_(new std::vector< ExpressionFunctionType >(values))
    , name_(nm)
  {
    // checks
    size_t totalSubdomains = 1;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      const auto& ll = (*lowerLeft_)[dd];
      const auto& ur = (*upperRight_)[dd];
      const auto& ne = (*numElements_)[dd];
      if (!(ll < ur))
        DUNE_THROW(Dune::RangeError, "lowerLeft has to be elementwise smaller than upperRight!");
      totalSubdomains *= ne;
    }
    if (values_->size() < totalSubdomains)
      DUNE_THROW(Dune::RangeError,
                 "values too small (is " << values_->size() << ", should be " << totalSubdomains << ")");
  } // Checkerboard(...)

  ExpressionCheckerboard(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".expressioncheckerboard";
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::unique_ptr< LocalfunctionType > local_global_function(const EntityType& entity) const override
  {
    // decide on the subdomain the center of the entity belongs to
    const auto center = entity.geometry().center();
    std::vector< size_t > whichPartition(dimDomain, 0);
    const auto& ll = *lowerLeft_;
    const auto& ur = *upperRight_;
    const auto& ne = *numElements_;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // for points that are on upperRight_[d], this selects one partition too much
      // so we need to cap this
      whichPartition[dd] = std::min(size_t(std::floor(ne[dd]*((center[dd] - ll[dd])/(ur[dd] - ll[dd])))),
                                    ne[dd] - 1);
    }
    size_t subdomain = 0;
    if (dimDomain == 1)
      subdomain = whichPartition[0];
    else if (dimDomain == 2)
      subdomain = whichPartition[0] + whichPartition[1]*ne[0];
    else
      subdomain = whichPartition[0] + whichPartition[1]*ne[0] + whichPartition[2]*ne[1]*ne[0];
    // return the component that belongs to the subdomain
    return std::unique_ptr< Localfunction >(new Localfunction(entity, (*values_)[subdomain]));
  } // ... local_function(...)

private:
  std::shared_ptr< const Common::FieldVector< DomainFieldType, dimDomain > > lowerLeft_;
  std::shared_ptr< const Common::FieldVector< DomainFieldType, dimDomain > > upperRight_;
  std::shared_ptr< const Common::FieldVector< size_t, dimDomain > > numElements_;
  std::shared_ptr< const std::vector< ExpressionFunctionType > > values_;
  std::string name_;
}; // class ExpressionCheckerboard


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
