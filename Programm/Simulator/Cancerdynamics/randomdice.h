#ifndef RANDOMDICE_H
#define RANDOMDICE_H

#include <random>
#include <mutex>


class RandomDice
{
public:
    // Methods
    RandomDice();
    double rollExpDist(double Lambda);
    double rollDiscrUnifDist(int min, int max);
    double rollContUnifDist(double);

    // TESTING
    double rollExpDist();

private:

    std::default_random_engine generator;
};

#endif // RANDOMDICE_H
