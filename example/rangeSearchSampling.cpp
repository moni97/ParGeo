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

  static const int dim = 2;
  size_t n = 10;

  /* Random number initializer */

  std::random_device rd;
  std::mt19937 gen(rd());

  /* Data generator */
  
  auto P = pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, 0, 1.0);
  pargeo::pointIO::writePointsToFile(P, "output.txt");

  /* Build a tree */

  pargeo::kdTree::node<dim, pargeo::point<dim>>* tree =
    pargeo::kdTree::build<dim, pargeo::point<dim>>(P, true);

  /* Rectangular range query example
     surrounding P[0] with half-length 0.1 */

  /* Orthogonal Range Search */
  std::cout << "\nOrthogonal Range Search\n";

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

  /* Orthogonal Range Count */
  std::cout << "\nOrthogonal Range Count\n";

  size_t count1 =
    pargeo::kdTree::orthogonalRangeSearchCount(tree, P[0], 0.1);

  std::cout << count1 << std::endl;

  /* Orthogonal range sample*/
  std::cout << "\nOrthogonal Range Sample\n";

  pargeo::kdTree::node<dim, pargeo::point<dim>>* sampleNode =  pargeo::kdTree::orthogonalRangeSample(tree, P[0], 0.1);

  if (sampleNode->size() > 1) {
    std::uniform_int_distribution<int> distribution(1, sampleNode->size());
    int sampleIndex = distribution(gen);
    pargeo::point<dim>* samplePoint = sampleNode->getItem(sampleIndex);
    std::cout << "\nsample point: " << *samplePoint;
  } else {
    int i;
    std::cout << sampleNode->size() << std::endl;
    for(i = 0; i < sampleNode->size(); ++i) {
        std::cout << *sampleNode->getItem(i) << std::endl;
    }
    std::cout << "\nsample point with leafsize 1: " << *sampleNode->getItem(0);
  }

  pargeo::kdTree::del(tree);
}
