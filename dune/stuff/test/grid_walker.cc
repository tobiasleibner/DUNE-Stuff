// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#if HAVE_DUNE_GRID

# include <dune/stuff/grid/walker.hh>
# include <dune/stuff/grid/provider/cube.hh>
# include <dune/stuff/common/parallel/partitioner.hh>
# include <dune/stuff/common/logstreams.hh>

# if DUNE_VERSION_NEWER(DUNE_COMMON,3,9) && HAVE_TBB // EXADUNE
#   include <dune/grid/utility/partitioning/seedlist.hh>
# endif

using namespace Dune::Stuff;
using namespace Dune::Stuff::Common;
using namespace Dune::Stuff::Grid;
using namespace std;

typedef testing::Types< Int<1>, Int<2>, Int<3> > GridDims;

template < class T >
struct GridWalkerTest : public ::testing::Test
{
  static const int griddim  = T::value;
  static const unsigned int level = 4;
  typedef Dune::YaspGrid<griddim> GridType;
  typedef typename GridType::LeafGridView GridViewType;
  typedef typename GridType::template Codim<0>::Entity EntityType;
  const DSG::Providers::Cube<GridType> grid_prv;
  GridWalkerTest()
    :grid_prv(0.f,1.f,level)
  {}

  void check() {
    const auto gv = grid_prv.grid().leafGridView();
    Walker<GridViewType> walker(gv);
    const auto correct_size = gv.size(0);
    atomic<size_t> count(0);
    auto counter = [&](const EntityType&){count++;};
    auto test1 = [&]{ walker.add(counter); walker.walk(false); };
    auto test2 = [&]{ walker.add(counter); walker.walk(true); };
    list<function<void()>> tests({ test1, test2 });
# if DUNE_VERSION_NEWER(DUNE_COMMON,3,9) // EXADUNE
    auto test3 = [&]{
                    IndexSetPartitioner<GridViewType> partitioner(gv.grid().leafIndexSet());
                    Dune::SeedListPartitioning<GridType, 0> partitioning(gv, partitioner);
                    walker.add(counter);
                    walker.walk(partitioning);
                  };
    tests.push_back(test3);
# endif // DUNE_VERSION_NEWER(DUNE_COMMON,3,9) // EXADUNE

    for (const auto& test : tests) {
      count = 0;
      test();
      EXPECT_EQ(count, correct_size);
    }
  }
};

TYPED_TEST_CASE(GridWalkerTest, GridDims);
TYPED_TEST(GridWalkerTest, Misc) {
  this->check();
}

# endif // HAVE_DUNE_GRID


