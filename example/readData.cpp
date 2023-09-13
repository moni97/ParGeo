#include <unordered_map>
#include <fstream>
#include <sstream>
#include "pargeo/point.h"
#include <iostream>
#include "kdTree/kdTree.h"
#include "pargeo/pointIO.h"
#include <string>

int main(int argc, char* argv[]) {
    std::string line, word;
    std::ifstream file ("../test/datasets/diabetesData.csv");
    bool notInitialized = true;

    int no_of_points, dim = 0;
    if (!file.is_open())
      throw std::runtime_error("Unable to open file");
    else {
        static const int numParameters = 3;
        using pointT = pargeo::point<numParameters>;
        double params[numParameters];
        std::string tmpString = "";
        std::vector<pointT> points;
        while ( getline(file,line) )
        {
            int i;
            std::stringstream inputString(line);
            for(i = 0;i < numParameters; ++i) {
                getline(inputString, tmpString, ',');
                params[i] = atof(tmpString.c_str());
                std::cout << "point at " << i << " " << params[i] << std::endl;
            }
            tmpString = "";
            getline(inputString, tmpString, ',');
            int group = atoi(tmpString.c_str());
            std::cout << "group: " <<group << std::endl;
            pargeo::point<numParameters> p(params, group);
            points.push_back(p);
            ++no_of_points;
        }
        std::cout << "num point: " << no_of_points << " point: " << points[2]; 
    }
}