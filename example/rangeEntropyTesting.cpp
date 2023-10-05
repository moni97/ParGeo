#include <assert.h>
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
#include <unordered_set>

int getDimension(std::string filename) {
    std::ifstream inputFile(filename);
    std::string line;
    bool isHeader = true;
    int numParameters = 0, dim;
    std::vector<std::string> headers;
    if (!inputFile.is_open())
      throw std::runtime_error("Unable to open file");
    else {
        while ( getline (inputFile,line) )
        {
            std::string tmpString;
            if(isHeader) {
                tmpString = "";
                std::stringstream inputString(line);
                while(getline(inputString, tmpString, ',')) {
                    ++numParameters;
                }
                isHeader = false;
                dim = numParameters - 1;
            } else {
                inputFile.close();
            }
        }
    }
    return dim;
}

template<int dim, class pointT = pargeo::point<dim>>
void rangeEntropy(std::string filename, bool printValues) {
    std::ifstream inputFile(filename);
    std::ofstream outputFile ("../test/datasets/popsimOutput1.txt");
    std::string line;
    std::vector<pointT> points;
    std::unordered_set<int> groupValues;
    int i, no_of_points=0;
    double coordinates[] = {0.1, 0.2, 0.3};

    if (!inputFile.is_open())
      throw std::runtime_error("Unable to open file");
    else {
        while ( getline (inputFile,line) )
        {
            std::string tmpString;
            std::stringstream inputString(line);
            double parameters[dim];

            tmpString = "";
            getline(inputString, tmpString, ',');
            int group = atoi(tmpString.c_str());


            for(i=0; i < dim; ++i) {
                tmpString = "";
                getline(inputString, tmpString, ',');
                parameters[i] = atof(tmpString.c_str());
            }

            pargeo::point<dim> p(parameters, group);
            points.push_back(p);
            ++no_of_points;
            groupValues.insert(group);
        }
    }

    auto P = parlay::sequence<pointT>(no_of_points);
    parlay::parallel_for (0, no_of_points, [&](size_t i) {
        P[i] = points[i];
    });

    if(printValues) std::cout << "Points read from file: " << P.size() << std::endl;
    //  std::cout << "sample point: " << P[0] << std::endl;

    pargeo::kdTree::node<dim, pointT>* tree =
        pargeo::kdTree::build<dim, pointT>(P, true);
    
    if(printValues) std::cout << "Tree created\n";
    int no_of_groups = groupValues.size();

    std::unordered_map<int, parlay::sequence<pointT>> groupsToPointsMap(no_of_groups);
    pargeo::kdTree::node<dim, pargeo::point<dim>>* treeForGroupMap[no_of_groups];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(0, no_of_points + 1);

    parlay::parallel_for (0, no_of_points, [&](size_t i) {
        groupsToPointsMap[points[i].attribute].push_back(points[i]);
    });

    // for (i = 0; i < no_of_points; ++i) {
    //     groupsToPointsMap[points[i].attribute].push_back(points[i]);
    // }

    // Create trees for colors
    for(i = 0; i < no_of_groups; ++i) {
        auto p = groupsToPointsMap[i];
        treeForGroupMap[i] = pargeo::kdTree::build<dim, pargeo::point<dim>>(p, true);
    }

    if(printValues) std::cout << "Trees for colors created\n";

    /* Running experiments */
    // double queryRadius[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25};
    double delta[] = {0.01, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5};
    // int arrSize = sizeof(queryRadius) / sizeof(queryRadius[0]);
    int deltaArraySize = sizeof(delta) / sizeof(delta[0]);
    
    /* Orthogonal range entropy */
    std::cout << "No. of points: " << no_of_points << " No. of colors: " << no_of_groups << std::endl; 
    int j=0, k = 0, count, numCanNodes = 0;
    double diff, currRadius = 0.01, currDelta, m;
    pointT queryPoint;
    if (outputFile.is_open()) {
        for(k = 0; k < 1; ++k) {
            // for(j = 1; j < 2; ++j) {
                queryPoint = P[k];
                for(i=0; i < 1; ++i) {
                    currDelta = delta[i];
                    std::cout << currDelta << ", ";
                    auto start = std::chrono::high_resolution_clock::now();
                    double calEntropy = pargeo::kdTree::orthogonalRangeEntropyAdditive(tree, queryPoint, currRadius, treeForGroupMap, no_of_groups, count, numCanNodes, currDelta, m);
                    auto end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> duration = end - start;

                    // double calMulEntropy = pargeo::kdTree::orthogonalRangeEntropyMultiplicative(tree, queryPoint, currRadius, treeForGroupMap, no_of_groups, count, numCanNodes, currDelta);
                    
                    auto start1 = std::chrono::high_resolution_clock::now();
                    double actualEntropy = pargeo::kdTree::rangeEntropyBruteForce(tree, queryPoint, currRadius, no_of_points);
                    auto end1 = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> duration1 = end1 - start1;
                    std::cout << count << ", " << m << ", " << calEntropy <<  ", " << actualEntropy << ", " << duration.count() << ", " << duration1.count() <<", " << duration.count() - duration1.count()<< std::endl;
                }
                // std::cout << "avg entropy: " << sumEnt / 8 << std::endl;
            // }
        }
    }

    /* Delete all the trees */
    groupsToPointsMap.clear();
    pargeo::kdTree::del(tree);
    for(i=0; i< no_of_groups; ++i) {
        pargeo::kdTree::del(treeForGroupMap[i]);
    }

    /* Delete all the trees */
    groupsToPointsMap.clear();
    pargeo::kdTree::del(tree);
    for(i=0; i< no_of_groups; ++i) {
        pargeo::kdTree::del(treeForGroupMap[i]);
    }

    inputFile.close();
}


int main(int argc, char* argv[]) {
    std::string filename = "../test/datasets/popsim.csv";
    int dimension = getDimension(filename);
    if(dimension == 2) rangeEntropy<2>(filename, true);
    else if(dimension == 3) rangeEntropy<3>(filename, true);
    else if(dimension == 4) rangeEntropy<4>(filename, true);
    else if(dimension == 5) rangeEntropy<5>(filename, true);
    else if(dimension == 6) rangeEntropy<6>(filename, true);
    else if(dimension == 7) rangeEntropy<7>(filename, true);
    else if(dimension == 8) rangeEntropy<8>(filename, true);
    else if(dimension == 9) rangeEntropy<9>(filename, true);
    return 0;
}
