// charVariety/utils/helpersTests.cpp

#define BOOST_TEST_MODULE MyTest
#include <boost/test/included/unit_test.hpp> 
#include "../utils/helpers.h" 
namespace utf = boost::unit_test;

#define DELTA 1e-6

BOOST_AUTO_TEST_CASE(testTangentSpaceAction)
{   
    Pair p(1, 2, 3, 4, 5, 6);

    Pair expectedX(1, 3, 1, 4, 6, 13);
    Pair foundX = DTX(p);

    Pair expectedY(3, 2, 5, 6, 5, 23);
    Pair foundY = DTY(p);

    Pair expectedX1(1, -1, 2, 4, 7, 5);
    Pair foundX1 = DTX1(p);
 
    Pair expectedY1(-1, 2, 1, 7, 5, 4);
    Pair foundY1 = DTY1(p); 

    BOOST_CHECK(expectedX == foundX);
    BOOST_CHECK(expectedY == foundY);
    BOOST_CHECK(expectedX1 == foundX1);
    BOOST_CHECK(expectedY1 == foundY1);

    Pair expectedXX1 = p;
    Pair foundXX1 = DTX(DTX1(p));
    Pair foundX1X = DTX1(DTX(p));

    Pair expectedYY1 = p;
    Pair foundYY1 = DTY(DTY1(p));
    Pair foundY1Y = DTY1(DTY(p)); 

    BOOST_CHECK(expectedXX1 == foundXX1);
    BOOST_CHECK(expectedXX1 == foundX1X);
    BOOST_CHECK(expectedYY1 == foundYY1);
    BOOST_CHECK(expectedYY1 == foundY1Y);
}

BOOST_AUTO_TEST_CASE(testAction)
{
    Pair p(1, 2, 3, 4, 5, 6);
    int maps[][2] = {
        {0},
        {1},
        {2},
        {3},
        {0, 2},
        {1, 3},
        {0, 1},
        {1, 0},
        {2, 0},
        {3, 1},
        {2, 3}
    };
    int nIter[] = {1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2};
    int nTests = sizeof(nIter) / sizeof(nIter[0]);

    Pair expectedImages[] = {
        {1, 3, 1, 4, 6, 13},
        {3, 2, 5, 6, 5, 23},
        {1, -1, 2, 4, 7, 5},
        {-1, 2, 1, 7, 5, 4},
        {1, 2, 3, 4, 5, 6},
        {1, 2, 3, 4, 5, 6},
        {3, 5, 13, 6, 23, 94},
        {1, 3, 2, 13, 6, 41},
        {1, 2, 3, 4, 5, 6},
        {1, 2, 3, 4, 5, 6},
        {-1, -3, 2, 7, 5, 5}
    };

    Pair expectedImageCals[] = {
        DTX(p),
        DTY(p),
        DTX1(p),
        DTY1(p),
        DTX(DTX1(p)),
        DTY(DTY1(p)),
        DTX(DTY(p)),
        DTY(DTX(p)),
        DTX1(DTX(p)),
        DTY1(DTY(p)),
        DTX1(DTY1(p))
    };

    for (int i = 0; i < nTests; i++){
        Pair foundImage = action(p, maps[i], nIter[i]);
        BOOST_CHECK(expectedImages[i] == foundImage);
        BOOST_CHECK(expectedImageCals[i] == foundImage);
    }
}

BOOST_AUTO_TEST_CASE(testAvgExpansion, * utf::tolerance(DELTA))
{
    Matrix a1(1, 2, 3, 4), a2(5, -6, -7, 8);
    Matrix A[] = {a1, a2};
    int nMatrix = sizeof(A) / sizeof(A[0]);
    double angles[] = {0.0, PI/2, 3*PI, 3*PI/2, PI/4};
    const int nAngle = sizeof(angles) / sizeof(angles[0]);

    double expectedAvgExps[] = {
        (0.5 * log(10) + 0.5 * log(74)) / double(nMatrix),
        (0.5 * log(20) + 0.5 * log(100)) / double(nMatrix),
        (0.5 * log(10) + 0.5 * log(74)) / double(nMatrix),
        (0.5 * log(20) + 0.5 * log(100)) / double(nMatrix),
        (0.5 * log(29) + 0.5 * log(1)) / double(nMatrix)
    };

    double expectedAvgExpCals[nAngle];
    for (int i = 0; i < nAngle; i++){
        expectedAvgExpCals[i] = 0.5 * log(A[0].expansionSq(angles[i])) + 0.5 * log(A[1].expansionSq(angles[i]));
        expectedAvgExpCals[i] /= double(nMatrix);
    }

    double foundAvgExps[nAngle];
    for (int i = 0; i < nAngle; i++){
        foundAvgExps[i] = avgLogExpansion(A, nMatrix, angles[i]);
        BOOST_TEST(expectedAvgExps[i] == foundAvgExps[i]);
        BOOST_TEST(expectedAvgExpCals[i] == foundAvgExps[i]);
    }
}

BOOST_AUTO_TEST_CASE(testMinAvgExpansion, * utf::tolerance(DELTA))
{
    DiagonalMatrix m1(2.0), m2(0.5);
    RotationMatrix m3(PI/2);
    Matrix m4(2.0, 1.0, 1.0, 2.0), 
           m5(2.0, -1.0, -1.0, 2.0);
    Matrix ms[][3] = {{m1}, {m1, m2}, {m3}, {m4}, {m4, m5}, {m1, m3}, {m1, m4, m5}};
    int nMatrices = 7, nMatrix[] = {1, 2, 1, 1, 2, 2, 3};
    double expectedMinAvgExps[] = {log(0.5), 0.0, 0.0, 0.0, 0.5*log(3.0), 0.5*log(0.5), log(2.5)/3};
    for (int i = 0; i < nMatrices; i++){
        double foundMinAvgExp = minAvgExpansionNewton(ms[i], nMatrix[i]);
        BOOST_TEST(expectedMinAvgExps[i] == foundMinAvgExp);
    }
}

// ---------template-----------
// BOOST_AUTO_TEST_CASE(my_boost_test4)
// {
//     double expected = 3.5;
//     double found = abb(-3.50001);
//     BOOST_CHECK(expected == found);
// }