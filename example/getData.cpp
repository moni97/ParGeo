#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include "pargeo/pointIO.h"
#include "pargeo/point.h"
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
pointT processData(std::string filename) {
    std::ifstream inputFile(filename);
    std::string line;
    std::vector<pointT> points;
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
            ++no_of_points;
            groupValues.insert(group);
        }
    }

    auto P = parlay::sequence<pointT>(no_of_points);
    parlay::parallel_for (0, no_of_points, [&](size_t i) {
        P[i] = points[i];
    });

    inputFile.close();
    return P;
}

auto main(int argc, char* argv[]) {
    static const int dimension = getDimension("../test/datasets/diabetesData2.csv");
    if (dimension == 2) {
        auto P = processData<2>("../test/datasets/diabetesData2.csv");
        return P;
    }
    else if (dimension == 3) { 
        auto P = processData<3>("../test/datasets/diabetesData2.csv");
        return P;
    }
    // else if (dimension == 4) { auto P = processData<4>("../test/datasets/diabetesData2.csv"); }
}