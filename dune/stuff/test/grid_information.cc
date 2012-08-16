#include "test_common.hh"

#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/math.hh>
#include <dune/stuff/common/logstreams.hh>

using namespace Dune::Stuff;
using namespace Dune::Stuff::Grid::Information;
using namespace std;

template < int i >
struct Int {
  static const int value = i;
};

typedef testing::Types< Int<1>, Int<2>, Int<3>, Int<4>, Int<5>, Int<6>, Int<7>, Int<8>, Int<9>, Int<10> > GridDims;

template < class T >
struct GridInfoTest : public ::testing::Test {
  static const int griddim  = T::value;
  static const int level = 1;
  typedef Dune::YaspGrid<griddim> GridType;
  Dune::shared_ptr<GridType> gridPtr;
  GridInfoTest()
    :gridPtr(Grid::Provider::UnitCube<GridType>(level).gridPtr())
  {}

  void check() {
    const Dimensions<GridType> dim(*gridPtr);
    const auto gv = gridPtr->leafView();
    const int entities = gv.size(0);
    EXPECT_DOUBLE_EQ(1.0/double(entities),dim.entity_volume.min());
    EXPECT_DOUBLE_EQ(dim.entity_volume.min(),dim.entity_volume.max());
    EXPECT_DOUBLE_EQ(dim.entity_volume.min(),dim.entity_volume.average());
    EXPECT_DOUBLE_EQ(1.0,dim.volumeRelation());
    const auto& dl = dim.coord_limits;
    for( int i : Common::Math::range(griddim) )
    {
      EXPECT_DOUBLE_EQ(dl[i].max(),1.0);
      EXPECT_DOUBLE_EQ(dl[i].min(),0.0);
      EXPECT_DOUBLE_EQ(dl[i].average(),0.5);
    }
    const Statistics st(gv);
    const int line = std::pow(2,level);
    EXPECT_EQ(std::pow(line,griddim-1)*2*griddim, st.numberOfBoundaryIntersections);
    EXPECT_EQ(entities*(2*griddim), st.numberOfIntersections);
    EXPECT_EQ(st.numberOfIntersections - st.numberOfBoundaryIntersections, st.numberOfInnerIntersections );
    EXPECT_EQ(griddim*2, maxNumberOfNeighbors(gv));
  }

  void print(std::ostream& out) {
    const auto& gv = gridPtr->leafView();
    ::print(gv, out);
  }
};

TYPED_TEST_CASE(GridInfoTest, GridDims);
TYPED_TEST(GridInfoTest, Misc) {
  this->check();
  this->print(Common::dev_null);

}

//TEST_F(GridInfoTest, Dimension) {
//  Grid::Provider::UnitCube cube;
//}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
