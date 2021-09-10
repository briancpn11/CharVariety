// charVariety/dataStructureTests.cpp

#define BOOST_TEST_MODULE MyTest
#include <boost/test/included/unit_test.hpp> 
#include "../utils/dataStructure.h" 
namespace utf = boost::unit_test;

#define DELTA 1e-6

BOOST_AUTO_TEST_CASE(testSquare, * utf::tolerance(DELTA))
{
    // test for sqr
    double param = -1.2;
    double foundSquare = sqr(param);
    double expectedSquare = 1.44;
    BOOST_TEST(expectedSquare == foundSquare);
}

BOOST_AUTO_TEST_CASE(testCompMax, * utf::tolerance(DELTA))
{
    // test for sqr
    double a = 1e-5, b = -1e-4, c = 1e-3;
    bool foundABC = compMax(a, b, c);
    bool foundBCA = compMax(b, c, a);
    bool foundCAB = compMax(c, a, b);
    bool expectedABC = 0, expectedBCA = 0, expectedCAB = 1;
    BOOST_TEST(expectedABC == foundABC);
    BOOST_TEST(expectedBCA == foundBCA);
    BOOST_TEST(expectedCAB == foundCAB);
}

BOOST_AUTO_TEST_CASE(testMatrixMethods, * utf::tolerance(DELTA))
{
    // test for Matrix product
    Matrix a(1, 2, 3, 4), b(5, -6, -7, 8);
    Matrix foundLeftProd = a*b;
    Matrix expectedLeftProd(-9, 10, -13, 14);
    BOOST_CHECK(expectedLeftProd == foundLeftProd); 

    // test for Matrix product
    Matrix foundRightProd = b*a;
    Matrix expectedRightProd(-13, -14, 17, 18);
    BOOST_CHECK(expectedRightProd == foundRightProd);

    // test for coeff0
    double foundCoeff0 = a.coeff0();
    double expectedCoeff0 = 15.0;
    BOOST_TEST(expectedCoeff0 == foundCoeff0);

    // test for coeff1
    double foundCoeff1 = a.coeff1();
    double expectedCoeff1 = -5.0;
    BOOST_TEST(expectedCoeff1 == foundCoeff1);

    // test for coeff2
    double foundCoeff2 = a.coeff2();
    double expectedCoeff2 = 14.0;
    BOOST_TEST(expectedCoeff2 == foundCoeff2);

    // tests for expansionSq
    double angles[] = {0.0, PI/2, 3*PI, 3*PI/2, PI/4};
    const int nAngle = sizeof(angles) / sizeof(angles[0]);
    double expectedExpSqs[] = {10.0, 20.0, 10.0, 20.0, 29.0};
    double foundExpSqs[nAngle];
    for (int i = 0; i < nAngle; i++){
        foundExpSqs[i] = a.expansionSq(angles[i]);
        BOOST_TEST(expectedExpSqs[i] == foundExpSqs[i]);
    }
    
    // test that expansionSq(x) = coeff0 + coeff1 * cos(2*x) + coeff2 * sin(2*x) by trig identities.
    double foundIDs[nAngle];
    for (int i = 0; i < nAngle; i++){
        foundIDs[i] = a.coeff0() + a.coeff1() * cos(2*angles[i]) + a.coeff2() * sin(2*angles[i]);
        BOOST_TEST(foundExpSqs[i] == foundIDs[i]);
    }
    
    // test for contractDir
    DiagonalMatrix m1(2.0), m2(0.4);
    RotationMatrix m3(PI/2);
    Matrix m4(2.0, 1.0, 1.0, 2.0), m5(3.0, -1.0, -1.0, 3.0);
    Matrix ms[] = {m1, m2, m3, m4, m5};
    const int nMatrix = sizeof(ms) / sizeof(ms[0]);
    double expectedContractDirs[] = {PI/2, 0.0, PI/4, 3*PI/4, PI/4};
    double foundContractDirs[nMatrix];

    for (int i = 0; i < nMatrix; i++){
        foundContractDirs[i] = ms[i].contractDir();
        BOOST_TEST(expectedContractDirs[i] == foundContractDirs[i]);
    }
}
