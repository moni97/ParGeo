#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "kdTree/kdTree.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"
#include "dataset/uniform.h"


typedef pargeo::kdTree::node<2, pargeo::point<2>> nodeT;

void printPargeoNodeDetails(nodeT* root) {
  std::cout << "\nNode details" << "\nNode size: " << root->size() << "\nNode min: " << root->getMin() << "\nNode max: " << root->getMax() << "\nisLeaf: " << root->isLeaf();
  int i;
  std::cout << "\nItems:";
  for (i = 0; i < root->size(); ++i) {
    std::cout << "\n" << *root->getItem(i);
  }
}
// In-order traversal: Left -> Root -> Right
void inOrderTraversal(nodeT* root) {
    if (root == nullptr) return;
    inOrderTraversal(root->L());
    printPargeoNodeDetails(root);
    inOrderTraversal(root->R());
}

int main(int argc, char* argv[]) {
  int dim = 2;
  double coordinates[2] = {2.0, 3.0}; // Coordinates: x=2.0, y=3.0
  int attributeValue = 42;
  pargeo::apoint<2> p(coordinates, attributeValue);
  std::cout << p[0] << p[1] << p.attribute;

  // if (dim == 2) {
    
  //   // parlay::sequence<pargeo::point<2>> P =
  //   //   pargeo::pointIO::readPointsFromFile<pargeo::point<2>>(iFile);

  //   int n = 30;
  //   auto P = pargeo::uniformInPolyPoints<2, pargeo::point<2>>(n, 0, 1.0);
    
  //   nodeT *tree = pargeo::kdTree::build<2, pargeo::point<2>>(P, true);
  //   inOrderTraversal(tree);
    
  // }
}

