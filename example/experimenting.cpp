#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "pargeo/point.h"
#include "kdTree/kdTree.h"
#include "pargeo/getTime.h"
#include "pargeo/pointIO.h"
#include "pargeo/parseCommandLine.h"


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
  pargeo::commandLine P(argc,argv,"[-k <param>] [-o <outFile>] <inFile>\n k = 1 will just return self.");
  char* iFile = P.getArgument(0);
  size_t k = P.getOptionIntValue("-k",1);
  char* oFile = P.getOptionValue("-o");

  int dim = pargeo::pointIO::readHeader(iFile);

  if (dim == 2) {
    
    parlay::sequence<pargeo::point<2>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<2>>(iFile);
    
    nodeT *tree = pargeo::kdTree::build<2, pargeo::point<2>>(P, true);
    inOrderTraversal(tree);
    // std::cout << "root node size: " << tree->size() << "\n";
    // if (tree->L() != 0) {
    //     nodeT *lefttNode = tree->L();
    //     std::cout << "Left node size: " << lefttNode->size() << "\n";
    //     if (lefttNode->L() != 0) {
    //       std::cout << "Left node size: " << lefttNode->L()->size() << "\n";
    //     }   
    // }
    // if (tree->R() != 0) {
    //     nodeT *rightNode = tree->L();
    //     std::cout << "Right node size: " << rightNode->size() << "\n";
    //     if (rightNode->R() != 0) {
    //       std::cout << "Right node size: " << rightNode->R()->size() << "\n";
    //     }
    // }
  }
}

