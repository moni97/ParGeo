#include <assert.h>
#include "dataset/uniform.h"
#include <iostream>
#include "kdTree/kdTree.h"
#include "pargeo/getTime.h"
#include "pargeo/point.h"
#include <random>
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"
#include <cstdlib>
#include <unordered_map>

int main(int argc, char* argv[]) {

  /* Parameters */

  static const int dim = 2;
  size_t n = 100;
  pargeo::point<dim> queryPoint;
  size_t queryRadius;
  size_t no_of_colors = 3;
  int colors[no_of_colors], i, j;
  pargeo::kdTree::node<dim, pargeo::point<dim>>* treeMap[no_of_colors+1];
  
  /* Random number initializer */

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> distribution(0, no_of_colors-1);

  /* Data generator */
  /* uniformInPolyPointsWithAttr => generates points with additional attribute */
  
  auto P = pargeo::uniformInPolyPointsWithAttr<dim, pargeo::point<dim>>(n, 0, no_of_colors, 1.0);

  /* Group the points based on the color */
  
  std::unordered_map<int, parlay::sequence<pargeo::point<dim>>> colorMap = pargeo::groupPoints<dim, pargeo::point<dim>>(P, no_of_colors);
  
  /* Create separate tree for each color */

  for(i = 0; i < no_of_colors; ++i) {
    auto p = colorMap[i];
    treeMap[i] = pargeo::kdTree::build<dim, pargeo::point<dim>>(p, true);
   }

  /* Query Parameters */ 
  
  queryPoint = P[8];
  queryRadius = 1.0;
  
  /* Build a tree */

  pargeo::kdTree::node<dim, pargeo::point<dim>>* tree =
    pargeo::kdTree::build<dim, pargeo::point<dim>>(P, true);

  /* Rectangular range query example
     surrounding P[0] with half-length 1.0 */

  /* Orthogonal Range Search */
  std::cout << "\nOrthogonal Range Search\n";

  parlay::sequence<pargeo::point<dim>*> elems2 =
    pargeo::kdTree::orthogonalRangeSearch(tree, queryPoint, queryRadius);
  
  /* Prints all the sample points */
  // j = 0;
  // for (pargeo::point<dim>* ptr : elems2) {
  //       std::cout << "Sample point " << j << ": ";
  //       for (int i = 0; i < dim; ++i) {
  //           std::cout << (*ptr)[i];
  //           if (i < dim - 1) {
  //               std::cout << ", ";
  //           }
  //       }
  //       ++j;
  //       std::cout << std::endl;
  //   }

  /* Orthogonal Range Count */
  std::cout << "\nOrthogonal Range Count\n";

  size_t count =
    pargeo::kdTree::orthogonalRangeSearchCount(tree, queryPoint, queryRadius);

  std::cout << count << std::endl;

  /* Orthogonal range sample*/

  std::cout << "\nOrthogonal Range Sample\n";
  pargeo::point<dim> samplPoint = pargeo::kdTree::orthogonalRangeSample(tree, queryPoint, queryRadius);
  std::cout << "sample point : "<< samplPoint << std::endl; 

  /* Orthogonal range entropy */

  double e = pargeo::kdTree::orthogonalRangeEntropy(tree, queryPoint, queryRadius, treeMap);
  std::cout << "Entropy e : "<< e << std::endl; 
  
  /* Delete all the trees */
  
  pargeo::kdTree::del(tree);
  for(i=0; i< no_of_colors; ++i) {
    pargeo::kdTree::del(treeMap[i]);
  }
}
