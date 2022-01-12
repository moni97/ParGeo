#include "dataset/uniform.h"
#include <iostream>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "kdTree/kdTree.h"

int main(int argc, char* argv[]) {

  /* Parameters */

  static const int dim = 4; // Data set dimensionality

  size_t n = 10000; // Number of data points

  /* Create a data set */

  auto P = pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, 0, 1.0);

  /* Build a kd-tree with leaf size 1 for the data set */

  pargeo::kdNode<dim, pargeo::point<dim>>* T =
    pargeo::buildKdTree<dim, pargeo::point<dim>>(P, true, 1);

  /* Compute WSPD with s = 2 */

  parlay::sequence<pargeo::wsp<pargeo::kdNode<dim, pargeo::point<dim>>>> pairs =
    pargeo::wspdParallel(T, 2);

  free(T);
}
