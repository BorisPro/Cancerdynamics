#include "randomdice.h"
#include <time.h>

RandomDice::RandomDice():generator(time(0)){
    generator.seed(time(0));
}

double RandomDice::rollExpDist(){
    std::exponential_distribution<double> distribution(1.);
    return distribution(generator);
}
double RandomDice::rollExpDist(double Lambda){
    std::exponential_distribution<double> distribution(Lambda);
    return distribution(generator);
}

double RandomDice::rollContUnifDist(double tmp){
    std::uniform_real_distribution<double> distribution(0,tmp);
//    std::uniform_real_distribution<double> UnifDist();
//    std::uniform_real_distribution<double>::param_type parameter(0,tmp);
//    UnifDist.param(parameter);
    return distribution(generator);
}
double RandomDice::rollDiscrUnifDist(int min, int max)
{
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
//    if(max > 0){
//        std::uniform_int_distribution<int> distribution(min, max);
//        return distribution(generator);
//    }
//    else return 0;
}
