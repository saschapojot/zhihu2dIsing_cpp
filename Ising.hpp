//
// Created by polya on 3/11/24.
//

#ifndef ZHIHU2DISING_ISING_HPP
#define ZHIHU2DISING_ISING_HPP
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <memory>
#include <array>
#include <random>
#include <chrono>
#include <boost/filesystem.hpp>

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>
#include <cstdlib>
#include <regex>
#include <fstream>

namespace fs = boost::filesystem;
class Ising{

public:
    Ising(double temperature,const int &partNum){
        this->T = temperature;
        this->beta=1/temperature;
        this->part=partNum;

    }

public:
    int N = 20;//   length of one direction
    double J = 1;

    double T=0;
    double beta=0;
    int part=0;

    std::vector<double>sRange{-1,1};
    int sweepNumInOneFlush=3000;// flush the results to python every sweepNumInOneFlush*L iterations
    int flushMaxNum=100;
    int dataNumTotal=15000;
    int lastFileNum=0;
public:
    ///
    /// @param cmd python execution string
    /// @return signal from the python
    static std::string execPython(const char* cmd);

///
    /// @param lag decorrelation length
    /// @param loopEq total loop numbers in reaching equilibrium
    void executionMC(const int& lag,const int & loopEq);// mc simulation without inquiring equilibrium after reaching equilibrium

    ///
    /// @param ferro is ferromagnetic
    /// @param lag decorrelation length
    /// @param loopTotal total mc steps
    void reachEqMC(bool& ferro, int &lag, int&loopTotal);// mc simulation while communicating with python to inquire equilibrium


    template<class T>
    static void printVec(const std::vector<T>& vec){
        for(int i=0;i<vec.size()-1;i++){
            std::cout<<vec[i]<<",";
        }
        std::cout<<vec[vec.size()-1]<<std::endl;
    }

    ///
    /// @param sCurr spin values
    /// @return total energy
    double ETot(const std::vector<std::vector<double>> & sCurr);


    ///
    /// @param sCurr spin values
    /// @return average of absolute values of spin in on configuration sCurr
    double sAvg(const std::vector<std::vector<double>> &sCurr);



    ///
    /// @param filename xml file name of vec
    ///@param vec vector to be saved
    static  void saveVecToXML(const std::string &filename,const std::vector<double> &vec);
};
#endif //ZHIHU2DISING_ISING_HPP
