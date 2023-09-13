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
#include <fstream>
#include <sstream>
#include <chrono>

int main(int argc, char* argv[]) {

  /* Parameters */

  static const int dim = 2;
  size_t n = 100000;
  pargeo::point<dim> queryPoint;
  size_t no_of_colors = 10;
  int colors[no_of_colors], i, j;
  pargeo::kdTree::node<dim, pargeo::point<dim>>* treeMap[no_of_colors+1];
  std::ofstream outputFile ("../test/datasets/output.txt");

  /* Random number initializer */
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> distribution(0, no_of_colors-1);

  /* Data generator */
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
  
  /* Build a tree */
  pargeo::kdTree::node<dim, pargeo::point<dim>>* tree =
    pargeo::kdTree::build<dim, pargeo::point<dim>>(P, true);

  /* Orthogonal range entropy */
  auto points = P;
  std::cout << tree->size() << std::endl;
  // std::cout << "No. of points: " <<  << " No. of colors: " << no_of_groups << std::endl; 
  int currRadius;
  double diff;
  // std::cout << "Rad,Count,Canonical size,Calculated Ent,Actual Ent,Difference,Time for cal ent,Time for act ent" << std::endl;
  if (outputFile.is_open()) {
  for(j = 0; j < 5; ++j) {
      int queryIdx = distribution(gen);
      queryPoint = P[queryIdx];
      std::cout<<j<<"\n";
      for(i=0; i < 1; ++i) {
          currRadius = 7;
          outputFile << currRadius << ", ";
          auto start = std::chrono::high_resolution_clock::now();
          double calEntropy = pargeo::kdTree::orthogonalRangeEntropy(tree, queryPoint, currRadius, colorMap, no_of_colors);
          auto end = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double> duration = end - start;

          auto start1 = std::chrono::high_resolution_clock::now();
          double actualEntropy = pargeo::kdTree::rangeEntropyBruteForce(tree, queryPoint, currRadius, n);
          auto end1 = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double> duration1 = end1 - start1;
          outputFile<< calEntropy << ", "<< actualEntropy << ", " << calEntropy - actualEntropy << ", " << duration.count() << ", " << duration1.count() << std::endl;
      }
    }
  }
  
  /* Delete all the trees */
  colorMap.clear();
  pargeo::kdTree::del(tree);
  for(i=0; i< no_of_colors; ++i) {
    pargeo::kdTree::del(treeMap[i]);
  }
}
