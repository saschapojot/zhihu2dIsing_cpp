//
// Created by polya on 3/11/24.
//


#include "Ising.hpp"


///
/// @param sCurr spin values
/// @return average of absolute values of spin in on configuration sCurr
double Ising::sAvg(const std::vector<std::vector<double>> &sCurr) {
    double sVal=0;
    for(const auto &vec:sCurr){
        for(const auto&val:vec){
            sVal+=val;
        }

    }

    sVal/=static_cast<double >(N*N);
    return std::abs(sVal);


}


///
/// @param sCurr spin values
/// @return total energy
double Ising::ETot(const std::vector<std::vector<double>> & sCurr){
    double HTmp=0;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){

            int leftInd=((j-1)%N+N)%N;
            int rightInd=((j+1)%N+N)%N;
            int upInd=((i-1)%N+N)%N;
            int downInd=((i+1)%N+N)%N;
//            printVec(std::vector<int>{leftInd,rightInd,upInd,downInd});
            double sij= sCurr[i][j];
            double sLeft=sCurr[i][leftInd];
            double sRight=sCurr[i][rightInd];
            double sUp=sCurr[upInd][j];
            double sDown= sCurr[downInd][j];

            HTmp+=sij*sUp+sij*sDown+sij*sLeft+sij*sRight;


        }
    }

    HTmp*=-0.5*J;

    return HTmp;



}



///
/// @param cmd python execution string
/// @return signal from the python
std::string Ising::execPython(const char* cmd){

    std::array<char, 4096> buffer; // Buffer to store command output
    std::string result; // String to accumulate output

    // Open a pipe to read the output of the executed command
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }

    // Read the output a chunk at a time and append it to the result string
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }

    return result; // Return the accumulated output


}

///
/// @param ferro is ferromagnetic
/// @param lag decorrelation length
/// @param loopTotal total mc steps
void Ising::reachEqMC(bool& ferro, int &lag, int&loopTotal) {// mc simulation while communicating with python to inquire equilibrium

    std::random_device rd;
    std::uniform_int_distribution<int> indsAll(0, 1);
    std::uniform_int_distribution<int> flipInds(0, N - 1);

    std::vector<std::vector<double>> sCurr;//init s
    for(int i=0;i<N;i++){
        std::vector<double> sOneRow;
        for(int j=0;j<N;j++){
            sOneRow.push_back(0);
        }
        sCurr.push_back(sOneRow);
    }


    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            sCurr[i][j]=sRange[indsAll(rd)];
        }
    }
//    std::cout<<"init s="<<sAvg(sCurr)<<std::endl;


    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);

    int flipNum = 0;
    int noFlipNum = 0;

    double ECurr=this->ETot(sCurr);
    double sValCurr=this->sAvg(sCurr);


    //output directory
    std::string outDir = "./part" + std::to_string(this->part)+"data"+"/part"+std::to_string(this->part)+ "/T" + std::to_string(this->T) + "/";
    std::string outEAllSubDir = outDir + "EAll/";

    std::string outSAllSubDir = outDir + "sAll/";
    if (!fs::is_directory(outEAllSubDir) || !fs::exists(outEAllSubDir)) {
        fs::create_directories(outEAllSubDir);
    }
    if (!fs::is_directory(outSAllSubDir) || !fs::exists(outSAllSubDir)) {
        fs::create_directories(outSAllSubDir);
    }


    const auto tMCStart{std::chrono::steady_clock::now()};
    int counter = 0;
    int fls = 0;
    bool active = true;
    std::regex stopRegex("stop");
    std::regex wrongRegex("wrong");
    std::regex ErrRegex("Err");
    std::regex lagRegex("\\d+");
    std::regex ferroRegex("ferro");
    std::regex eqRegex("equilibrium");

    std::smatch matchEAvgStop;
    std::smatch matchEAvgWrong;
    std::smatch matchEAvgErr;
    std::smatch matchEAvgLag;
    std::smatch matchEAvgFerro;
    std::smatch matchEAvgEq;

//    std::smatch matchSAvgStop;
//    std::smatch matchSAvgWrong;
//    std::smatch matchSAvgErr;
//    std::smatch matchSAvgLag;
//    std::smatch matchSAvgFerro;
//    std::smatch matchSAvgEq;


    while (fls < this->flushMaxNum and active == true) {
        std::vector<double>EAll;
        std::vector<double> sAll;
        int loopStart = fls * this->sweepNumInOneFlush * N*N;
        for (int i = 0; i < this->sweepNumInOneFlush * N*N; i++) {
            //perform a flip
            auto sNext=std::vector<std::vector<double>>(sCurr);
            int flip_i=flipInds(rd);
            int flip_j=flipInds(rd);
            sNext[flip_i][flip_j]*=-1;
            double ENext=this->ETot(sNext);
            double DeltaE=(ENext-ECurr);
            //decide if flip is accepted
            if(DeltaE<=0){
                sCurr=std::vector<std::vector<double>>(sNext);
                ECurr=ENext;
                flipNum++;
            }else{
                double r = distUnif01(e2);
                if(r <= std::exp(-this->beta * DeltaE)){
                    sCurr=std::vector<std::vector<double>>(sNext);
                    ECurr=ENext;
                    flipNum++;
                }else{
                    noFlipNum++;
                }

            }

            sAll.push_back(sValCurr);
            EAll.push_back(ECurr);
            counter += 1;
        }//end for loop
        int loopEnd = loopStart + this->sweepNumInOneFlush * N*N - 1;

        std::string filenameMiddle = "loopStart" + std::to_string(loopStart) +
                                     "loopEnd" + std::to_string(loopEnd) + "T"  + "J" + std::to_string(this->J);
        std::string outEFileTmp = outEAllSubDir + filenameMiddle + ".EAll.xml";
        this->saveVecToXML(outEFileTmp,EAll);
        std::string outSFileTmp = outSAllSubDir + filenameMiddle + ".sAll.xml";
        this->saveVecToXML(outSFileTmp,sAll);
        const auto tflushEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_seconds{tflushEnd - tMCStart};
        std::cout << "flush " << fls << std::endl;
        std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;

        //communicate with python to inquire equilibrium

        //inquire equilibrium of EAvg
        std::string commandEAvg = "python3 checkVec.py " + outEAllSubDir;
        //inquire equilibrium of s
//        std::string commandSAvg="python3 check_s_Vec.py "+outSAllSubDir;
        std::string resultEAvg;
//        std::string resultSAvg;
        if (fls % 3 == 2) {
            try {
                resultEAvg = this->execPython(commandEAvg.c_str());
//                resultSAvg=this->execPython(commandSAvg.c_str());

                std::cout << "EAvg message from python: " << resultEAvg << std::endl;
//                std::cout << "sAvg message from python: " << resultSAvg << std::endl;

            } catch (const std::exception &e) {
                std::cerr << "Error: " << e.what() << std::endl;
                std::exit(10);
            }
            catch (...) {
                // Handle any other exceptions
                std::cerr << "Error" << std::endl;
                std::exit(11);
            }

            // parse result
            if (std::regex_search(resultEAvg, matchEAvgErr, ErrRegex)) {
                std::cout << "error encountered" << std::endl;
                std::exit(12);
            }

            if (std::regex_search(resultEAvg, matchEAvgWrong, wrongRegex)) {
                std::exit(13);
            }

//        bool ferro= false;

            if (std::regex_search(resultEAvg, matchEAvgStop, stopRegex) ) {
                if (std::regex_search(resultEAvg, matchEAvgFerro, ferroRegex) ) {
                    active = false;
                    ferro = true;
                }


            }


            if (std::regex_search(resultEAvg, matchEAvgEq, eqRegex)) {
                if (std::regex_search(resultEAvg, matchEAvgLag, lagRegex) ) {
                    std::string lagStrEAvg = matchEAvgLag.str(0);

                    int lagEAvg = std::stoi(lagStrEAvg);

                    lag=lagEAvg;

                    active = false;
                }

            }
        }

        fls++;

    }//end of while loop
    std::ofstream outSummary(outDir + "summary.txt");
    loopTotal = flipNum + noFlipNum;
    const auto tMCEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
    outSummary << "total mc time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    outSummary << "total sweep number: " << static_cast<int>(loopTotal / (N*N)) << std::endl;
    outSummary << "total loop number: " << loopTotal << std::endl;

    outSummary << "flip number: " << flipNum << std::endl;
    outSummary << "no flip number: " << noFlipNum << std::endl;

    outSummary << "equilibrium reached: " << !active << std::endl;
    outSummary<<"ferro: "<<ferro<<std::endl;
    outSummary << "lag=" << lag << std::endl;
    outSummary.close();

}



///
/// @param filename xml file name of vec
///@param vec vector to be saved
  void Ising::saveVecToXML(const std::string &filename,const std::vector<double> &vec){

    std::ofstream ofs(filename);
    boost::archive::xml_oarchive oa(ofs);
    oa & BOOST_SERIALIZATION_NVP(vec);

}


///
/// @param lag decorrelation length
/// @param loopEq total loop numbers in reaching equilibrium
void Ising::executionMC(const int& lag,const int & loopEq) {// mc simulation without inquiring equilibrium after reaching equilibrium
    double lagDB = static_cast<double>(lag);
    double loopEqDB = static_cast<double >(loopEq);
    int counter=0;

    int remainingDataNum = this->dataNumTotal - static_cast<int>(std::floor(loopEqDB / lagDB * 2 / 3));
    int remainingLoopNum = remainingDataNum * lag;

    if (remainingLoopNum <= 0) {
        return;
    }
    double remainingLoopNumDB = static_cast<double>(remainingLoopNum);
    double NNDB = static_cast<double>(N*N);
    int remainingSweepNum = std::ceil(remainingLoopNumDB / NNDB);
    double remainingSweepNumDB = static_cast<double >(remainingSweepNum);
    double sweepNumInOneFlushDB = static_cast<double >(sweepNumInOneFlush);
    double remainingFlushNumDB = std::ceil(remainingSweepNumDB / sweepNumInOneFlushDB);
    int remainingFlushNum = static_cast<int>(remainingFlushNumDB);

    //init
    std::random_device rd;
    std::uniform_int_distribution<int> indsAll(0, 1);
    std::uniform_int_distribution<int> flipInds(0, N - 1);
    std::vector<std::vector<double>> sCurr;//init s
    for(int i=0;i<N;i++){
        std::vector<double> sOneRow;
        for(int j=0;j<N;j++){
            sOneRow.push_back(0);
        }
        sCurr.push_back(sOneRow);
    }

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            sCurr[i][j]=sRange[indsAll(rd)];
        }
    }
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);
    int flipNum = 0;
    int noFlipNum = 0;

    double ECurr=this->ETot(sCurr);
    double sValCurr=this->sAvg(sCurr);
    //output directory
    std::string outDir = "./part" + std::to_string(this->part)+"data"+"/part"+std::to_string(this->part)+ "/T" + std::to_string(this->T) + "/";
    std::string outEAllSubDir = outDir + "EAll/";

    std::string outSAllSubDir = outDir + "sAll/";
    const auto tMCStart{std::chrono::steady_clock::now()};

    std::cout<<"remaining flush number: "<<remainingFlushNum<<std::endl;

    for (int fls = 0; fls < remainingFlushNum; fls++) {
        std::vector<double>EAll;
        std::vector<double> sAll;
        int loopStart =loopEq+ fls * this->sweepNumInOneFlush * N*N;
        for (int i = 0; i < this->sweepNumInOneFlush * N*N; i++) {
            //perform a flip
            auto sNext=std::vector<std::vector<double>>(sCurr);
            int flip_i=flipInds(rd);
            int flip_j=flipInds(rd);
            sNext[flip_i][flip_j]*=-1;
            double ENext=this->ETot(sNext);
            double DeltaE=(ENext-ECurr);
            //decide if flip is accepted
            if(DeltaE<=0){
                sCurr=std::vector<std::vector<double>>(sNext);
                ECurr=ENext;
                flipNum++;
            }else{
                double r = distUnif01(e2);
                if(r <= std::exp(-this->beta * DeltaE)){
                    sCurr=std::vector<std::vector<double>>(sNext);
                    ECurr=ENext;
                    flipNum++;
                }else{
                    noFlipNum++;
                }

            }
            sAll.push_back(sValCurr);
            EAll.push_back(ECurr);
            counter += 1;
        }//end of sweeps in 1 flush
        int loopEnd = loopStart + this->sweepNumInOneFlush * N*N - 1;
        std::string filenameMiddle = "loopStart" + std::to_string(loopStart) +
                                     "loopEnd" + std::to_string(loopEnd) + "T"  + "J" + std::to_string(this->J);
        std::string outEFileTmp = outEAllSubDir + filenameMiddle + ".EAll.xml";
        this->saveVecToXML(outEFileTmp,EAll);
        std::string outSFileTmp = outSAllSubDir + filenameMiddle + ".sAll.xml";
        this->saveVecToXML(outSFileTmp,sAll);


        const auto tflushEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_seconds{tflushEnd - tMCStart};
        std::cout << "flush " << fls << std::endl;
        std::cout << "time elapsed: " << elapsed_seconds.count() / 3600.0 << " h" << std::endl;


    }//end flush loop
    std::ofstream outSummary(outDir + "summaryAfterEq.txt");
    int loopTotal = flipNum + noFlipNum;


    const auto tMCEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
    outSummary << "total mc time: " << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    outSummary << "total sweep number: " << static_cast<int>(loopTotal / (N*N)) << std::endl;
    outSummary << "total loop number: " << loopTotal << std::endl;

    outSummary << "flip number: " << flipNum << std::endl;
    outSummary << "no flip number: " << noFlipNum << std::endl;


    outSummary.close();



}