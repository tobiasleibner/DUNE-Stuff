﻿// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTION_INTERFACE_HH
#define DUNE_STUFF_FUNCTION_INTERFACE_HH

#include <vector>
#include <memory>
#include <string>
#include <ostream>
#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#if HAVE_DUNE_GRID
# include <dune/grid/io/file/vtk.hh>
# include <dune/stuff/common/filesystem.hh>
#endif

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

#if HAVE_DUNE_PDELAB
#include <dune/typetree/nodetags.hh>
#include <dune/pdelab/common/function.hh>
#endif

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/type_utils.hh>

namespace Dune {
namespace Stuff {
namespace internal {


template< class F >
struct is_localizable_function_helper
{
  DSC_has_typedef_initialize_once(EntityType)
  DSC_has_typedef_initialize_once(DomainFieldType)
  DSC_has_typedef_initialize_once(RangeFieldType)
  DSC_has_static_member_initialize_once(dimDomain)
  DSC_has_static_member_initialize_once(dimRange)
  DSC_has_static_member_initialize_once(dimRangeCols)

  static const bool is_candidate = DSC_has_typedef(EntityType)< F >::value
                                   && DSC_has_typedef(DomainFieldType)< F >::value
                                   && DSC_has_typedef(RangeFieldType)< F >::value
                                   && DSC_has_static_member(dimDomain)< F >::value
                                   && DSC_has_static_member(dimRange)< F >::value
                                   && DSC_has_static_member(dimRangeCols)< F >::value;
}; // class is_localizable_function_helper


} // namespace internal


// forwards, includes are below
template< class F, bool candidate = internal::is_localizable_function_helper< F >::is_candidate >
struct is_localizable_function;


namespace Functions {
namespace internal {


// additional argument for member functions to differentiate between dimRangeCols = 1 and dimRangeCols > 1 by overloading
template< size_t rangeDimCols >
struct ChooseVariant {};


} // namespace internal


#if HAVE_DUNE_GRID


template< class GridViewType, size_t dimRange, size_t dimRangeCols = 1 >
class VisualizationAdapter;


#endif // HAVE_DUNE_GRID


template< class MinuendType, class SubtrahendType >
class Difference;

template< class LeftSummandType, class RightSummandType >
class Sum;

template< class LeftSummandType, class RightSummandType >
class Product;

template< class FunctionImp >
class Divergence;


} // namespace Functions
namespace Tags {


class LocalizableFunction {};


} // namespace Tags


/**
 *  \brief Interface for a set of globalvalued functions, which can be evaluated locally on one Entity.
 */
template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class LocalfunctionSetInterface
{
  static_assert(EntityImp::dimension == domainDim, "Dimensions do not match!");

  template< class RangeFieldType, size_t dimRange, size_t dimRangeCols>
  struct RangeTypeSelector
  {
    typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimRangeCols > type;
  };

  template< class RangeFieldType, size_t dimRange >
  struct RangeTypeSelector< RangeFieldType, dimRange, 1 >
  {
    typedef Dune::FieldVector< RangeFieldType, dimRange > type;
  };

  template< size_t dimDomain, class RangeFieldType, size_t dimRange, size_t dimRangeCols >
  struct JacobianRangeTypeSelector
  {
    typedef Dune::FieldVector< Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain >, dimRangeCols > type;
  };

  template< size_t dimDomain, class RangeFieldType, size_t dimRange >
  struct JacobianRangeTypeSelector< dimDomain, RangeFieldType, dimRange, 1 >
  {
    typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > type;
  };

public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp                                  DomainFieldType;
  static const size_t                                     dimDomain = domainDim;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;

  typedef RangeFieldImp                                                              RangeFieldType;
  static const size_t                                                                dimRange = rangeDim;
  static const size_t                                                                dimRangeCols = rangeDimCols;
  typedef typename RangeTypeSelector< RangeFieldType, dimRange, dimRangeCols >::type RangeType;
  typedef typename JacobianRangeTypeSelector
      < dimDomain, RangeFieldType, dimRange, dimRangeCols >::type                    JacobianRangeType;

  LocalfunctionSetInterface(const EntityType& ent)
    : entity_(ent)
  {}

  virtual ~LocalfunctionSetInterface() {}

  virtual const EntityType& entity() const
  {
    return entity_;
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual size_t size() const = 0;

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& /*xx*/, std::vector< RangeType >& /*ret*/) const = 0;

  virtual void jacobian(const DomainType& /*xx*/, std::vector< JacobianRangeType >& /*ret*/) const = 0;
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  std::vector< RangeType > evaluate(const DomainType& xx) const
  {
    std::vector< RangeType > ret(size(), RangeType(0));
    evaluate(xx, ret);
    return ret;
  }

  std::vector< JacobianRangeType > jacobian(const DomainType& xx) const
  {
    std::vector< JacobianRangeType > ret(size(), JacobianRangeType(0));
    jacobian(xx, ret);
    return ret;
  }
  /* @} */

protected:
  bool is_a_valid_point(const DomainType&
#ifndef DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS
                                          xx
#else
                                          /*xx*/
#endif
                                                ) const
  {
#ifndef DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS
    const auto& reference_element = ReferenceElements< DomainFieldType, dimDomain >::general(entity().type());
    return reference_element.checkInside(xx);
#else // DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS
    return true;
#endif
  }

  const EntityType& entity_;
}; // class LocalfunctionSetInterface


/**
 *  \brief  Interface for functions, which can be evaluated locally on one Entity.
 */
template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class LocalfunctionInterface
  : public LocalfunctionSetInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
{
  typedef LocalfunctionSetInterface
      < EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols > BaseType;
public:
  typedef EntityImp EntityType;

  typedef typename BaseType::DomainFieldType   DomainFieldType;
  static const size_t                          dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType        DomainType;
  typedef typename BaseType::RangeType         RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  LocalfunctionInterface(const EntityType& ent)
    : BaseType(ent)
  {}

  virtual ~LocalfunctionInterface() {}

  /**
   * \defgroup haveto ´´These methods have to be implemented in addition to the ones required from the BaseType.''
   * @{
   **/
  virtual void evaluate(const DomainType& /*xx*/, RangeType& /*ret*/) const = 0;

  virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/) const = 0;
  /* @} */

  /**
   * \defgroup providedbase ´´These methods are provided by the interface to please LocalfunctionSetInterface.''
   * @{
   **/
  virtual size_t size() const final
  {
    return 1;
  }

  virtual void evaluate(const DomainType& xx, std::vector< RangeType >& ret) const final
  {
    assert(ret.size() >= 1);
    evaluate(xx, ret[0]);
  }

  virtual void jacobian(const DomainType& xx, std::vector< JacobianRangeType >& ret) const final
  {
    assert(ret.size() >= 1);
    jacobian(xx, ret[0]);
  }
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  RangeType evaluate(const DomainType& xx) const
  {
    RangeType ret(0);
    evaluate(xx, ret);
    return ret;
  }

  JacobianRangeType jacobian(const DomainType& xx) const
  {
    JacobianRangeType ret(0);
    jacobian(xx, ret);
    return ret;
  }

  //! evaluate at N quadrature points into vector of size >= N
  void evaluate(const Dune::QuadratureRule< DomainFieldType, dimDomain >& quadrature, std::vector< RangeType >& ret)
  {
    assert(ret.size() >= quadrature.size());
    std::size_t i = 0;
    for (const auto& point : quadrature)
      evaluate(point.position(), ret[i++]);
  }

  //! jacobian at N quadrature points into vector of size >= N
  void jacobian(const Dune::QuadratureRule< DomainFieldType, dimDomain >& quadrature,
                std::vector< JacobianRangeType >& ret)
  {
    assert(ret.size() >= quadrature.size());
    std::size_t i = 0;
    for (const auto& point : quadrature)
      jacobian(point.position(), ret[i++]);
  }
  /* @} */
}; // class LocalfunctionInterface


class IsLocalizableFunction {};


/**
 * \brief Interface for functions which provide a LocalfunctionInterface for an entity.
 */
template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class LocalizableFunctionInterface
  : public IsLocalizableFunction
  , public Tags::LocalizableFunction
{
  typedef LocalizableFunctionInterface
      < EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols > ThisType;
public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const size_t    dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;
  static const size_t   dimRange = rangeDim;
  static const size_t   dimRangeCols = rangeDimCols;

  typedef LocalfunctionInterface
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols > LocalfunctionType;

  typedef typename LocalfunctionType::DomainType        DomainType;
  typedef typename LocalfunctionType::RangeType         RangeType;
  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;

  static const bool available = false;

  typedef Functions::Difference< ThisType, ThisType > DifferenceType;
  typedef Functions::Sum< ThisType, ThisType >        SumType;
  typedef Functions::Divergence< ThisType >           DivergenceType;

  virtual ~LocalizableFunctionInterface() {}

  static std::string static_id()
  {
    return "stuff.function";
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityType& /*entity*/) const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string type() const
  {
    return "stuff.function";
  }

  virtual std::string name() const
  {
    return "stuff.function";
  }
  /* @} */

  DifferenceType operator-(const ThisType& other) const
  {
    return DifferenceType(*this, other);
  }

  SumType operator+(const ThisType& other) const
  {
    return SumType(*this, other);
  }

  template< class OtherType >
  typename std::enable_if< is_localizable_function< OtherType >::value,
                           Functions::Product< ThisType, OtherType > >::type
  operator*(const OtherType& other) const
  {
    return Functions::Product< ThisType, OtherType >(*this, other);
  }

  DivergenceType divergence() const
  {
    return DivergenceType(*this);
  }

#if HAVE_DUNE_GRID
  /**
   * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default. This means that the grid you
   *        see in the visualization is a refinement of the actual grid!
   */
  template< class GridViewType >
  void visualize(const GridViewType& grid_view,
                 const std::string path,
                 const bool subsampling = true,
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    if (path.empty()) DUNE_THROW(RangeError, "Empty path given!");
    const auto directory = DSC::directoryOnly(path);
    const auto filename = DSC::filenameOnly(path);
    auto adapter = std::make_shared< Functions::VisualizationAdapter< GridViewType, dimRange, dimRangeCols > >(*this);
    std::unique_ptr< VTKWriter< GridViewType > > vtk_writer =
        subsampling ? DSC::make_unique< SubsamplingVTKWriter< GridViewType > >(grid_view, VTK::nonconforming)
                    : DSC::make_unique< VTKWriter< GridViewType > >(grid_view, VTK::nonconforming);
    vtk_writer->addVertexData(adapter);
    DSC::testCreateDirectory(directory);
    if (MPIHelper::getCollectiveCommunication().size() == 1)
      vtk_writer->write(path, vtk_output_type);
    else
      vtk_writer->pwrite(filename, directory, "", vtk_output_type);
  } // ... visualize(...)
#endif // HAVE_DUNE_GRID

  virtual void report(std::ostream& out, const std::string prefix = "") const
  {
    out << prefix << "function '" << name() << "' (of type " << type() << ")";
  }

private:
  template< class T >
  friend std::ostream& operator<<(std::ostream& /*out*/, const ThisType& /*function*/);
}; // class LocalizableFunctionInterface


template< class E, class D, size_t d, class R, size_t r, size_t rC >
std::ostream& operator<<(std::ostream& out, const LocalizableFunctionInterface< E, D, d, R, r, rC >& function)
{
  function.report(out);
  return out;
} // ... operator<<(...)


template < class OtherEntityImp, class GlobalFunctionImp >
class TransferredGlobalFunction;


/**
 * base class for global matrix-valued valued functions that provides automatic local functions via LocalizableFunctionInterface
 */
template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class GlobalFunctionInterface
  : public LocalizableFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
{
  typedef LocalizableFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
      BaseType;
  typedef GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols >
      ThisType;
public:
  typedef typename BaseType::LocalfunctionType LocalfunctionType;
  typedef typename BaseType::DomainType        DomainType;
  typedef typename BaseType::RangeType         RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  virtual ~GlobalFunctionInterface() {}

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& xx, RangeType& ret) const = 0;

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "This does not make sense yet for matrix-valued functions!");
  }

  virtual RangeType evaluate(const DomainType& xx) const
  {
    RangeType ret;
    evaluate(xx, ret);
    return ret;
  }

  virtual JacobianRangeType jacobian(const DomainType& xx) const
  {
    JacobianRangeType ret;
    jacobian(xx, ret);
    return ret;
  }

  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityImp& entity) const override final
  {
    return Common::make_unique< Localfunction >(entity, *this);
  }

  virtual std::string type() const override
  {
    return "stuff.globalfunction";
  }

  virtual std::string name() const override
  {
    return "stuff.globalfunction";
  }

private:
  class Localfunction
    : public LocalfunctionType
  {
  public:
    Localfunction(const EntityImp& entity_in, const ThisType& global_function)
      : LocalfunctionType(entity_in)
      , geometry_(entity_in.geometry())
      , global_function_(global_function)
    {}

    virtual ~Localfunction() {}

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian(xx_global, ret);
    }

    virtual size_t order() const override final
    {
      return global_function_.order();
    }

  private:
      const typename EntityImp::Geometry geometry_;
      const ThisType& global_function_;
  }; //class Localfunction

public:
  template < class OtherEntityImp >
  struct Transfer {
    typedef TransferredGlobalFunction< OtherEntityImp, ThisType> Type;
  };

  template < class OtherEntityImp >
  typename Transfer<OtherEntityImp>::Type transfer() const {
    return typename Transfer<OtherEntityImp>::Type(*this);
  }
}; // class GlobalFunctionInterface


/**
 * base class for global valued functions that provides automatic local functions via LocalizableFunctionInterface
 */
template< class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public LocalizableFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
#if HAVE_DUNE_FEM
  , public Dune::Fem::Function< Dune::Fem::FunctionSpace< DomainFieldImp, RangeFieldImp, domainDim, rangeDim >,
                                GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > >
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
  , public TypeTree::LeafNode
  , public PDELab::FunctionInterface< PDELab::FunctionTraits< DomainFieldImp, domainDim, FieldVector<DomainFieldImp, domainDim >
                                                            , RangeFieldImp, rangeDim, FieldVector<RangeFieldImp, rangeDim > >
                                    , GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > >
#endif
{
  typedef LocalizableFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > BaseType;
  typedef GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >      ThisType;
public:
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::RangeType;
  using typename BaseType::LocalfunctionType;
  using typename BaseType::EntityType;
#if HAVE_DUNE_FEM
  typedef typename Dune::Fem::Function
      < Dune::Fem::FunctionSpace< DomainFieldImp, RangeFieldImp, domainDim, rangeDim >
      , GlobalFunctionInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > >::JacobianRangeType
                                               JacobianRangeType;
#else
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
#endif

  virtual ~GlobalFunctionInterface() {}

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& x, RangeType& ret) const = 0;

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "You have to imlement it if you intend to use it!");
  }

  virtual RangeType evaluate(const DomainType& xx) const
  {
    RangeType ret;
    evaluate(xx, ret);
    return ret;
  }

  virtual JacobianRangeType jacobian(const DomainType& xx) const
  {
    JacobianRangeType ret;
    jacobian(xx, ret);
    return ret;
  }

  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityImp& entity) const override final
  {
    return Common::make_unique< Localfunction >(entity, *this);
  }

  virtual std::string type() const override
  {
    return "stuff.globalfunction";
  }

  virtual std::string name() const override
  {
    return "stuff.globalfunction";
  }

private:
  class Localfunction
    : public LocalfunctionType
  {
  public:
    Localfunction(const EntityImp& entity_in, const ThisType& global_function)
      : LocalfunctionType(entity_in)
      , geometry_(entity_in.geometry())
      , global_function_(global_function)
    {}

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian(xx_global, ret);
    }

    virtual size_t order() const override final
    {
      return global_function_.order();
    }

  private:
      const typename EntityImp::Geometry geometry_;
      const ThisType& global_function_;
  }; //class Localfunction

public:
  template < class OtherEntityImp >
  struct Transfer {
    typedef TransferredGlobalFunction< OtherEntityImp, ThisType> Type;
  };

  template < class OtherEntityImp >
  typename Transfer<OtherEntityImp>::Type transfer() const {
    return typename Transfer<OtherEntityImp>::Type(*this);
  }
}; // class GlobalFunctionInterface< ..., 1 >

template < class OtherEntityImp, class GlobalFunctionImp >
class TransferredGlobalFunction
  : public GlobalFunctionInterface< OtherEntityImp, typename GlobalFunctionImp::DomainFieldType,
                                    GlobalFunctionImp::dimDomain, typename GlobalFunctionImp::RangeFieldType,
                                    GlobalFunctionImp::dimRange, GlobalFunctionImp::dimRangeCols >
{
  typedef GlobalFunctionInterface< OtherEntityImp, typename GlobalFunctionImp::DomainFieldType,
                                   GlobalFunctionImp::dimDomain, typename GlobalFunctionImp::RangeFieldType,
                                   GlobalFunctionImp::dimRange, GlobalFunctionImp::dimRangeCols > BaseType;

public:
  TransferredGlobalFunction(const GlobalFunctionImp& function)
    : function_(function)
  {}

  virtual size_t order() const
  {
    return function_.order();
  }

  virtual void evaluate(const typename BaseType::DomainType& x, typename BaseType::RangeType& ret) const
  {
    function_.evaluate(x, ret);
  }

private:
  const GlobalFunctionImp& function_;
}; // class TransferredGlobalFunction


/**
 *  \brief  Interface for functions, which can be partially evaluated locally on one Entity
 */
template< class EntityImp, class DomainFieldImp, size_t domainDim,
          class RangeEntityImp, class RangeDomainFieldImp, size_t rangeDomainDim,
          class RangeRangeFieldImp, size_t rangeRangeDim, size_t rangeRangeDimCols >
class FunctionValuedLocalfunctionInterface
{
  typedef FunctionValuedLocalfunctionInterface< EntityImp, DomainFieldImp, domainDim,
                                                RangeEntityImp, RangeDomainFieldImp, rangeDomainDim,
                                                RangeRangeFieldImp, rangeRangeDim, rangeRangeDimCols >         ThisType;

  template< class R, size_t r, size_t rC>
  struct RangeRangeTypeSelector
  {
    typedef Dune::FieldMatrix< R, r, rC > type;
  };

  template< class R, size_t r >
  struct RangeRangeTypeSelector< R, r, 1 >
  {
    typedef Dune::FieldVector< R, r > type;
  };

  template< size_t d, class R, size_t r, size_t rC >
  struct RangeJacobianRangeTypeSelector
  {
    typedef Dune::FieldVector< Dune::FieldMatrix< R, r, d >, rC > type;
  };

  template< size_t d, class R, size_t r >
  struct RangeJacobianRangeTypeSelector< d, R, r, 1 >
  {
    typedef Dune::FieldMatrix< R, r, d > type;
  };

public:
  typedef EntityImp EntityType;
  typedef RangeEntityImp RangeEntityType;

  typedef DomainFieldImp DomainFieldType;
  static const size_t dimDomain = domainDim;
  typedef RangeDomainFieldImp RangeDomainFieldType;
  static const size_t rangeDimDomain = rangeDomainDim;
  static const size_t entireDimDomain = dimDomain + rangeDimDomain;

  typedef RangeRangeFieldImp RangeRangeFieldType;
  static const size_t   rangeDimRange = rangeRangeDim;
  static const size_t   rangeDimRangeCols = rangeRangeDimCols;
  typedef Dune::FieldVector< DomainFieldType, dimDomain >            DomainType;
  typedef Dune::FieldVector< RangeDomainFieldType, rangeDimDomain >  RangeDomainType;
  typedef RangeRangeTypeSelector< RangeRangeFieldType, rangeDimRange, rangeDimRangeCols > RangeRangeType;
  typedef RangeJacobianRangeTypeSelector< rangeDimDomain, RangeRangeFieldType, rangeDimRange, rangeDimRangeCols > RangeJacobianRangeType;
  typedef LocalizableFunctionInterface< RangeEntityType,
                                        RangeDomainFieldType, rangeDimDomain,
                                        RangeRangeFieldType, rangeDimRange, rangeDimRangeCols >   RangeType;
  typedef LocalizableFunctionInterface< RangeEntityType,
                                        RangeDomainFieldType, rangeDimDomain,
                                        RangeRangeFieldType, rangeDimRangeCols, dimDomain >       JacobianRangeType;

  FunctionValuedLocalfunctionInterface(const EntityType& ent)
    : entity_(ent)
  {}

  virtual ~FunctionValuedLocalfunctionInterface() {}

  virtual const EntityType& entity() const
  {
    return entity_;
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& /*xx*/, std::shared_ptr< const RangeType >& ret) const = 0;

  virtual void jacobian(const DomainType& /*xx*/, std::shared_ptr< const JacobianRangeType >& ret) const = 0;

  /** @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/

  std::shared_ptr< const RangeType > evaluate(const DomainType& xx) const
  {
    std::shared_ptr< const RangeType > ret;
    evaluate(xx, ret);
    return ret;
  }

  std::shared_ptr< const JacobianRangeType > jacobian(const DomainType& xx) const
  {
    std::shared_ptr< const JacobianRangeType > ret;
    jacobian(xx, ret);
    return ret;
  }

  /** @} */

protected:
  bool is_a_valid_point(const DomainType&
#ifndef DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS
                                          xx
#else
                                          /*xx*/
#endif
                                                ) const
  {
#ifndef DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS
    const auto& reference_element = ReferenceElements< DomainFieldType, dimDomain >::general(entity().type());
    return reference_element.checkInside(xx);
#else // DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS
    return true;
#endif
  }

private:
    const EntityType& entity_;
}; // class FunctionValuedLocalfunctionInterface

/**
 * \brief Interface for functions f(x,u) that are localizable in x but not in u
 */
template< class EntityImp, class DomainFieldImp, size_t domainDim,
          class RangeEntityImp, class RangeDomainFieldImp, size_t rangeDomainDim,
          class RangeRangeFieldImp, size_t rangeRangeDim, size_t rangeRangeDimCols >
class FunctionValuedFunctionInterface
{
  typedef FunctionValuedFunctionInterface< EntityImp, DomainFieldImp, domainDim,
                                                RangeEntityImp, RangeDomainFieldImp, rangeDomainDim,
                                                RangeRangeFieldImp, rangeRangeDim, rangeRangeDimCols >         ThisType;
public:
  typedef FunctionValuedLocalfunctionInterface< EntityImp,
                                          DomainFieldImp, domainDim,
                                          RangeEntityImp, RangeDomainFieldImp, rangeDomainDim,
                                          RangeRangeFieldImp, rangeRangeDim, rangeRangeDimCols >      LocalfunctionType;

  typedef typename LocalfunctionType::DomainType                   DomainType;
  typedef typename LocalfunctionType::RangeDomainType              RangeDomainType;
  typedef typename LocalfunctionType::RangeRangeType               RangeRangeType;
  typedef typename LocalfunctionType::RangeJacobianRangeType       RangeJacobianRangeType;
  typedef typename LocalfunctionType::RangeType                    RangeType;
  typedef typename LocalfunctionType::JacobianRangeType            JacobianRangeType;
  typedef typename LocalfunctionType::EntityType                   EntityType;
  typedef typename LocalfunctionType::RangeEntityType              RangeEntityType;

  typedef typename LocalfunctionType::DomainFieldType              DomainFieldType;
  static const size_t dimDomain = LocalfunctionType::dimDomain;

  typedef typename LocalfunctionType::RangeDomainFieldType         RangeDomainFieldType;
  typedef typename LocalfunctionType::RangeRangeFieldType          RangeRangeFieldType;
  static const size_t rangeDimDomain = LocalfunctionType::rangeDimDomain;
  static const size_t rangeDimRange = LocalfunctionType::rangeDimRange;
  static const size_t rangeDimRangeCols = LocalfunctionType::rangeDimRangeCols;

  static const bool available = false;

  typedef Functions::Difference< ThisType, ThisType > DifferenceType;
  typedef Functions::Sum< ThisType, ThisType >        SumType;
  typedef Functions::Divergence< ThisType >           DivergenceType;

  virtual ~FunctionValuedFunctionInterface() {}

  static std::string static_id()
  {
    return "stuff.function";
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityType& /*entity*/) const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string type() const
  {
    return "stuff.function";
  }

  virtual std::string name() const
  {
    return "stuff.function";
  }
  /* @} */

  DifferenceType operator-(const ThisType& other) const
  {
    return DifferenceType(*this, other);
  }

  SumType operator+(const ThisType& other) const
  {
    return SumType(*this, other);
  }

  template< class OtherType >
  typename std::enable_if< is_localizable_function< OtherType >::value,
                           Functions::Product< ThisType, OtherType > >::type
  operator*(const OtherType& other) const
  {
    return Functions::Product< ThisType, OtherType >(*this, other);
  }

  DivergenceType divergence() const
  {
    return DivergenceType(*this);
  }

  virtual void report(std::ostream& out, const std::string prefix = "") const
  {
    out << prefix << "function '" << name() << "' (of type " << type() << ")";
  }

private:
  template< class T >
  friend std::ostream& operator<<(std::ostream& /*out*/, const ThisType& /*function*/);
}; // class FunctionValuedFunctionInterface


/**
 * \brief Interface for functions f(x,u) that are localizable in x and global evaluable in u
 */
template< class EntityImp, class DomainFieldImp, size_t domainDim,
          class RangeEntityImp, class RangeDomainFieldImp, size_t rangeDomainDim,
          class RangeRangeFieldImp, size_t rangeRangeDim, size_t rangeRangeDimCols >
class GlobalFunctionValuedFunctionInterface
    : public FunctionValuedFunctionInterface < EntityImp, DomainFieldImp, domainDim,
                                                    RangeEntityImp, RangeDomainFieldImp, rangeDomainDim,
                                                    RangeRangeFieldImp, rangeRangeDim, rangeRangeDimCols >
{
  typedef GlobalFunctionValuedFunctionInterface< EntityImp, DomainFieldImp, domainDim,
                                                 RangeEntityImp, RangeDomainFieldImp, rangeDomainDim,
                                                 RangeRangeFieldImp, rangeRangeDim, rangeRangeDimCols >   ThisType;
  typedef FunctionValuedFunctionInterface< EntityImp, DomainFieldImp, domainDim,
                                           RangeEntityImp, RangeDomainFieldImp, rangeDomainDim,
                                           RangeRangeFieldImp, rangeRangeDim, rangeRangeDimCols >         BaseType;
public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeDomainType;
  using typename BaseType::RangeRangeType;
  using typename BaseType::RangeJacobianRangeType;
  using typename BaseType::EntityType;
  using typename BaseType::RangeEntityType;
  typedef typename BaseType::LocalfunctionType BaseLocalFunctionType;

  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;

  using typename BaseType::RangeDomainFieldType;
  using typename BaseType::RangeRangeFieldType;
  using BaseType::rangeDimDomain;
  using BaseType::rangeDimRange;
  using BaseType::rangeDimRangeCols;

  typedef GlobalFunctionInterface< RangeEntityType,
                                   RangeDomainFieldType, rangeDimDomain,
                                   RangeRangeFieldType, rangeDimRange, rangeDimRangeCols >   RangeType;
  typedef GlobalFunctionInterface< RangeEntityType,
                                   RangeDomainFieldType, rangeDimDomain,
                                   RangeRangeFieldType, rangeDimRangeCols, dimDomain >       JacobianRangeType;

  static const bool available = false;

  virtual ~GlobalFunctionValuedFunctionInterface() {}

private:
  class Localfunction
      : public BaseLocalFunctionType
  {
    typedef typename BaseLocalFunctionType::RangeType         BaseRangeType;
    typedef typename BaseLocalFunctionType::JacobianRangeType BaseJacobianRangeType;

    typedef GlobalFunctionInterface< RangeEntityType,
                                     RangeDomainFieldType, rangeDimDomain,
                                     RangeRangeFieldType,  rangeDimRange, rangeDimRangeCols >   RangeType;
    typedef GlobalFunctionInterface< RangeEntityType,
                                     RangeDomainFieldType, rangeDimDomain,
                                     RangeRangeFieldType,  rangeDimRange, dimDomain >       JacobianRangeType;
  public:
    Localfunction(const EntityImp& entity_in)
      : BaseLocalFunctionType(entity_in)
    {}

    virtual void evaluate(const DomainType& /*xx*/, std::shared_ptr< const RangeType >& ret) const = 0;

    virtual void jacobian(const DomainType& /*xx*/, std::shared_ptr< const JacobianRangeType >& ret) const = 0;

    virtual void evaluate(const DomainType& xx, std::shared_ptr< const BaseRangeType >& ret) const override final
    {
      std::shared_ptr< const RangeType > derived;
      evaluate(xx, derived);
      ret = std::shared_ptr< const BaseRangeType >(derived);
    }

    virtual void jacobian(const DomainType& xx, std::shared_ptr< const BaseJacobianRangeType >& ret) const override final
    {
      std::shared_ptr< const JacobianRangeType > derived;
      jacobian(xx, derived);
      ret = std::shared_ptr< const BaseJacobianRangeType >(derived);
    }

    std::shared_ptr< const RangeType > evaluate(const DomainType& xx) const
    {
      std::shared_ptr< const RangeType > ret;
      evaluate(xx, ret);
      return ret;
    }

    std::shared_ptr< const JacobianRangeType > jacobian(const DomainType& xx) const
    {
      std::shared_ptr< const JacobianRangeType > ret;
      jacobian(xx, ret);
      return ret;
    }

    void evaluate(const DomainType& xx, const RangeDomainType& uu, RangeRangeType& ret) const
    {
      std::shared_ptr< const RangeType > range_func = evaluate(xx);
      ret = range_func->evaluate(uu);
    }

    RangeRangeType evaluate(const DomainType& xx, const RangeDomainType& uu) const
    {
      RangeRangeType ret(RangeRangeFieldType(0));
      evaluate(xx, uu, ret);
      return ret;
    }

    using BaseLocalFunctionType::jacobian;

    void jacobian(const DomainType& xx, const RangeDomainType& uu, RangeJacobianRangeType& ret) const
    {
      std::shared_ptr< const JacobianRangeType > range_func = jacobian(xx);
      ret = range_func->jacobian(uu);
    }

    RangeJacobianRangeType jacobian(const DomainType& xx, const RangeDomainType& uu) const
    {
      RangeJacobianRangeType ret(RangeRangeFieldType(0));
      jacobian(xx, uu, ret);
      return ret;
    }
  }; //class Localfunction

public:
  typedef Localfunction LocalfunctionType;

  virtual std::unique_ptr< LocalfunctionType > local_global_function(const EntityType& /*entity*/) const = 0;

  virtual std::unique_ptr< BaseLocalFunctionType > local_function(const EntityType& entity) const override final
  {
    return std::unique_ptr< BaseLocalFunctionType >(std::move(local_global_function(entity)));
  }


  virtual std::string type() const override
  {
    return "stuff.globalfunctionvaluedfunction";
  }

  virtual std::string name() const override
  {
    return "stuff.globalfunctionvaluedfunction";
  }

}; // class GlobalFunctionValuedFunctionInterface



//! Utility to generate a complete Function Type from an existing one and a template
template <class FunctionImp, template <class, class, size_t, class, size_t, size_t> class OutTemplate>
struct FunctionTypeGenerator {
  typedef OutTemplate<typename FunctionImp::EntityType, typename FunctionImp::DomainFieldType,
  FunctionImp::dimDomain, typename FunctionImp::RangeFieldType,
  FunctionImp::dimRange, FunctionImp::dimRangeCols> type;
};


template< class F >
struct is_localizable_function< F, true >
  : public std::is_base_of< LocalizableFunctionInterface< typename F::EntityType,
                                                          typename F::DomainFieldType, F::dimDomain,
                                                          typename F::RangeFieldType,  F::dimRange, F::dimRangeCols >
                          , F >
{};


template< class F >
struct is_localizable_function< F, false >
  : public std::false_type
{};


} // namespace Stuff
} // namespace Dune

#include "default.hh"
#include "combined.hh"
#include "derived.hh"

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
