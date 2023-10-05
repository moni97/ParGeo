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
#include <unordered_set>

int main(int argc, char* argv[]) {

    // Variables
    static const int dim = 3;
    using pointT = pargeo::point<dim>;
    std::string line;
    std::ifstream file ("../test/datasets/diabetesData.csv");
    std::ofstream outputFile ("../test/datasets/mAanlysisDiabetes8.txt");
    std::ofstream pointsFile ("../test/datasets/pointsFile.txt");
    int no_of_points = 0,
        i, j;
    pointT queryPoint;
    std::vector<pointT> points;
    std::unordered_set<int> groupValues;

    // File reading
    if (!file.is_open())
      throw std::runtime_error("Unable to open file");
    else {
        while ( getline (file,line) )
        {
            std::string tmpString;
            double param1, param2, param3, param4, param5, param6;
            std::stringstream inputString(line);
            
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
            int group = atoi(tmpString.c_str());

            // tmpString = "";
            // getline(inputString, tmpString, ',');
            // param4 = atof(tmpString.c_str());

            // tmpString = "";
            // getline(inputString, tmpString, ',');
            // param5 = atof(tmpString.c_str());

            // tmpString = "";
            // getline(inputString, tmpString, ',');
            // param6 = atof(tmpString.c_str());

            double coordinates[dim] = {param1, param2, param3};
            pargeo::point<dim> p(coordinates, group);
            points.push_back(p);
            if(pointsFile.is_open()) {
                pointsFile<<p<<std::endl;
            }
            ++no_of_points;
            groupValues.insert(group);
        }
        
        auto P = parlay::sequence<pointT>(no_of_points);
        parlay::parallel_for (0, no_of_points, [&](size_t i) {
            P[i] = points[i];
        });
        
        pargeo::kdTree::node<dim, pointT>* tree =
            pargeo::kdTree::build<dim, pointT>(P, true);
        
        int no_of_groups = groupValues.size();

        std::unordered_map<int, parlay::sequence<pointT>> groupsToPointsMap(no_of_groups);
        pargeo::kdTree::node<dim, pargeo::point<dim>>* treeForGroupMap[no_of_groups];

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> distribution(0, no_of_points + 1);

        for (i = 0; i < no_of_points; ++i) {
            groupsToPointsMap[points[i].attribute].push_back(points[i]);
        }

        // Create trees for colors
        for(i = 0; i < no_of_groups; ++i) {
            auto p = groupsToPointsMap[i];
            treeForGroupMap[i] = pargeo::kdTree::build<dim, pargeo::point<dim>>(p, true);
        }
       
        /* Running experiments */
        double delta[] = {0.01, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5};
        int m[] = {10, 100, 500, 1000, 5000, 10000};
        int mArraySize = sizeof(m) / sizeof(m[0]);
        
        /* Orthogonal range entropy */
        auto points = P;
        std::cout << "No. of points: " << no_of_points << " No. of colors: " << no_of_groups << std::endl; 
        int i=0, j=0, count, numCanNodes = 0;
        double diff, currRadius = 3, currDelta;
        pointT queryPoint;
        int k=0;
        if (outputFile.is_open()) {
            for(k = 0; k < 5; ++k) {
                std::cout << "In interation: " << k << std::endl;
                queryPoint = P[9];
                for(i=0; i < mArraySize; ++i) {
                    currDelta = 0.1;
                    outputFile << currDelta << ", ";
                    // auto start = std::chrono::high_resolution_clock::now();
                    double calEntropy = pargeo::kdTree::orthogonalRangeEntropyAdditive(tree, queryPoint, currRadius, treeForGroupMap, no_of_groups, count, numCanNodes, currDelta, m[i]);
                    // auto end = std::chrono::high_resolution_clock::now();
                    // std::chrono::duration<double> duration = end - start;

                    // double calMulEntropy = pargeo::kdTree::orthogonalRangeEntropyMultiplicative(tree, queryPoint, currRadius, treeForGroupMap, no_of_groups, count, numCanNodes, currDelta);
                    
                    // auto start1 = std::chrono::high_resolution_clock::now();
                    double actualEntropy = pargeo::kdTree::rangeEntropyBruteForce(tree, queryPoint, currRadius, no_of_points);
                    // auto end1 = std::chrono::high_resolution_clock::now();
                    // std::chrono::duration<double> duration1 = end1 - start1;
                    outputFile << count << ", " << m[i]  << ", " << numCanNodes << ", " << calEntropy << ", "<< actualEntropy << ", " << calEntropy - actualEntropy << std::endl;
                    // ", " << duration.count() << ", " << duration1.count() <<", " << duration.count() - duration1.count()<< std::endl;
                }
            }
        }

        /* Delete all the trees */
        groupsToPointsMap.clear();
        pargeo::kdTree::del(tree);
        for(i=0; i< no_of_groups; ++i) {
            pargeo::kdTree::del(treeForGroupMap[i]);
        }

        file.close();
        outputFile.close();
    }
}
