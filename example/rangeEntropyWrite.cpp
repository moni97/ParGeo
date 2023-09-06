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

    // Variables
    std::string line;
    std::ifstream file ("../test/datasets/adult.csv");
    std::ofstream outputFile ("../test/datasets/output.txt");
    const int dim = 6;
    int no_of_points = 0,
        no_of_groups = 10,
        i, j;
    using pointT = pargeo::point<6>;
    pointT queryPoint;
    std::vector<pointT> points;
    std::unordered_map<int, parlay::sequence<pointT>> groupsToPointsMap(no_of_groups);
    pargeo::kdTree::node<dim, pargeo::point<dim>>* treeForGroupMap[no_of_groups];
    
    // File reading
    if (!file.is_open())
      throw std::runtime_error("Unable to open file");
    else {
        while ( getline (file,line) )
        {
            std::string tmpString;
            double param1, param2, param3, param4, param5, param6;
            std::stringstream inputString(line);

            getline(inputString, tmpString, ',');
            getline(inputString, tmpString, ',');
            getline(inputString, tmpString, ',');

            tmpString = "";
            getline(inputString, tmpString, ',');
            int group = atoi(tmpString.c_str());

            tmpString = "";
            getline(inputString, tmpString, ',');
            param1 = atof(tmpString.c_str());

            tmpString = "";
            getline(inputString, tmpString, ',');
            param2 = atof(tmpString.c_str());

            tmpString = "";
            getline(inputString, tmpString, ',');
            param3 = atof(tmpString.c_str());

            tmpString = "";
            getline(inputString, tmpString, ',');
            param4 = atof(tmpString.c_str());

            tmpString = "";
            getline(inputString, tmpString, ',');
            param5 = atof(tmpString.c_str());

            tmpString = "";
            getline(inputString, tmpString, ',');
            param6 = atof(tmpString.c_str());

            double coordinates[6] = {param1, param2, param3, param4, param5, param6};
            pargeo::point<6> p(coordinates, group);
            points.push_back(p);
            groupsToPointsMap[p.attribute].push_back(p);
            ++no_of_points;
        }

        auto P = parlay::sequence<pointT>(no_of_points);
        parlay::parallel_for (0, no_of_points, [&](size_t i) {
            P[i] = points[i];
        });
        
        // Build tree
        pargeo::kdTree::node<dim, pointT>* tree =
            pargeo::kdTree::build<dim, pointT>(P, true);
        
        std::cout << "Building tree\n"; 

        // Create trees for colors
        for(i = 0; i < no_of_groups; ++i) {
            auto p = groupsToPointsMap[i];
            treeForGroupMap[i] = pargeo::kdTree::build<dim, pargeo::point<dim>>(p, true);
        }

        std::cout << "Tree for each color created\n";

        // Random number generator
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> distribution(0, no_of_points);

        /* Running experiments */
        int queryRadius[] = {1,5,7,9,11,13,15,17,19,21};
        int arrSize = sizeof(queryRadius) / sizeof(queryRadius[0]);

        std::cout << "No. of points: " << no_of_points << " No. of colors: " << no_of_groups << std::endl; 
        
        /* Orthogonal range entropy */
        int i=0, j=0, rad;
        double diff;
        pointT queryPoint;

        for(j = 0; j < 51; ++j) {
            queryPoint = P[j];
            // rad = 3;

            for(i=0; i < arrSize; ++i) {
                rad = queryRadius[i];
                std::cout << rad << ", ";
                auto start = std::chrono::high_resolution_clock::now();
                double calEntropy = pargeo::kdTree::orthogonalRangeEntropy(tree, queryPoint, rad, treeForGroupMap);
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;

                auto start1 = std::chrono::high_resolution_clock::now();
                double actualEntropy = pargeo::kdTree::rangeEntropyBruteForce(tree, queryPoint, rad, no_of_points);
                auto end1 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration1 = end1 - start1;
                diff += calEntropy - actualEntropy;
                std::cout << calEntropy << ", "<< actualEntropy << ", " << calEntropy - actualEntropy << ", " << duration.count() << ", " << duration1.count() << std::endl;
            }
        // outputFile << "difference: " << diff/50 << std::endl;
        // outputFile.close();
        }

        /* Delete all the trees */
        groupsToPointsMap.clear();
        pargeo::kdTree::del(tree);
        for(i=0; i< no_of_groups; ++i) {
            pargeo::kdTree::del(treeForGroupMap[i]);
        }
        file.close();
    }
}
