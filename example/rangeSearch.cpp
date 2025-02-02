#include <assert.h>
#include "dataset/uniform.h"
#include <iostream>
#include "kdTree/kdTree.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include <random>
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"

int main(int argc, char* argv[]) {

  /* Parameters */

  static const int dim = 2; // Data set dimensionality

  size_t n = 1000; // Number of data points

  std::random_device rd;
  std::mt19937 gen(rd());

  auto P = pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, 0, 1.0);
  pargeo::pointIO::writePointsToFile(P, "output.txt");

  /* Build a tree */

  pargeo::kdTree::node<dim, pargeo::point<dim>>* tree =
    pargeo::kdTree::build<dim, pargeo::point<dim>>(P, true);

  /* Spherical range query example
     surrounding P[0] with radius 0.1 */

   parlay::sequence<pargeo::point<dim>*> elems1 =
     pargeo::kdTree::rangeSearch(tree, P[0], 0.1);

  /* Rectangular range query example
     surrounding P[0] with half-length 0.1 */

  /* Orthogonal Range Search */

  std::cout << "\nOrthogonal Range Search wiht 0.1\n";

  parlay::sequence<pargeo::point<dim>*> elems2 =
    pargeo::kdTree::orthogonalRangeSearch(tree, P[0], 0.1);

  for (pargeo::point<dim>* ptr : elems2) {
        for (int i = 0; i < dim; ++i) {
            std::cout << (*ptr)[i];
            if (i < dim - 1) {
                std::cout << ", ";
            }
        }
        std::cout << std::endl;
    }

  std::cout << "\nOrthogonal Range Count";

  size_t count1 =
    pargeo::kdTree::orthogonalRangeSearchCount(tree, P[0], 0.1);

  // size_t count2 =
  //   pargeo::kdTree::orthogonalRangeSearchCount(tree, P[0], 0.5);

  // size_t count3 =
  //   pargeo::kdTree::orthogonalRangeSearchCount(tree, P[0], 0.8);
  std::cout << "\nNo of points with 0.1: " << count1;
  //  << "\nNo of points with 0.5: " << count2 << "\nNo of points with 0.8: " << count3 << std::endl;

  /* Orthogonal range sample*/
  pargeo::kdTree::node<dim, pargeo::point<dim>>* sampleNode =  pargeo::kdTree::orthogonalRangeSample(tree, P[0], 0.1);

  if (sampleNode->size() > 1) {
    std::uniform_int_distribution<int> distribution(1, sampleNode->size());
    int sampleIndex = distribution(gen);
    pargeo::point<dim>* samplePoint = sampleNode->getItem(sampleIndex);
    std::cout << "\nsample point: " << *samplePoint;
  } else {
    std::cout << "\nsample point with leafsize 1: " << *sampleNode->getItem(1);
  }

  pargeo::kdTree::del(tree);
}
