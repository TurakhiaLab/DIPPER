#ifndef MASHDISTANCE_HPP
#include "mashDistance.hpp"
#endif

float jaccardEstimate(std::priority_queue<hash_t> A, std::priority_queue<hash_t> B, size_t sketchSize)
{
    // std::cout << A.size() << "--" << B.size() << "\n";

    std::vector<hash_t> ADash, BDash, ADashUnionBDash;
    while (!A.empty()) 
    {
        ADash.push_back(A.top());
        A.pop();
    }

    while (!B.empty()) 
    {
        BDash.push_back(B.top());
        B.pop();
    }

    float inter = 0;
    float uni = 0;
    int ADashPointer=ADash.size()-1;
    int BDashPointer=BDash.size()-1;

    while (true)
    {
        if (ADashPointer==-1) break;
        if (BDashPointer==-1) break;
        if (uni>sketchSize) break;

        if (ADash[ADashPointer] == BDash[BDashPointer]) { inter++;uni++;ADashPointer--;BDashPointer--; }
        else if (ADash[ADashPointer] < BDash[BDashPointer]) { uni++; ADashPointer--; }
        else { uni++; BDashPointer--; }

    }

    while (uni<sketchSize && (ADashPointer>=0 || BDashPointer>=0))
    {
        if(ADashPointer>=0) {uni++; ADashPointer--;}
        if(BDashPointer>=0) {uni++; BDashPointer--;}
    }

    // std::cout << inter << "\t" << uni << "\t";
    return (inter/uni);
}

float mashDistance(std::priority_queue<hash_t>& A, std::priority_queue<hash_t>& B, int k, size_t sketchSize)
{
    float j = 1 - jaccardEstimate(A, B, sketchSize);
    if (j==0) j=0.000001;
    if (j==1) j=0.999999;
    // std::cout << j << "--" << ((1/(float)k)*log(2*j/(1+j))) << std::endl;
    return (-(1/(float)k)*log(2*j/(1+j)));
}