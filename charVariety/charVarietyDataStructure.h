// charVarietyDataStructure.h

#ifndef CHARVARIETY_DATA_STRUCTURE
#define CHARVARIETY_DATA_STRUCTURE

#include <boost/numeric/interval.hpp>
#include <iostream>

#define EPS 1e-6
#define PI 3.1415926535

// define Interval

struct full_rounding:
  	boost::numeric::interval_lib::rounded_arith_opp<double>
{
	public:
	# define GENR_FUNC(name) \
		double name##_down(double x) { return name(x); } \
		double name##_up  (double x) { return name(x); }
		GENR_FUNC(exp)
		GENR_FUNC(log)
		GENR_FUNC(sin)
		GENR_FUNC(cos)
};

namespace IntervalNS {
	using namespace boost;
	using namespace numeric;
	using namespace interval_lib;
	typedef save_state<full_rounding> R;
	typedef checking_strict<double> P;
	typedef interval<double, policies<R, P> > Interval;
};

typedef IntervalNS::Interval Interval;

// junk
// template<class os_t>
// os_t& operator<<(os_t &os, const Interval &a) {
//   os << '[' << a.lower() << ',' << a.upper() << ']';
//   return os;
// }

// std::ostream& operator<<(std::ostream& os, const Interval& r) {
//   os << "[" << r.lower() << "," << r.upper() << "]";
//   return os;
// }
// end junk

// double sqr(double);

inline double sqr(double a){
  	return a*a;
}

// define a pair of point (x, y, z) in R^3 and a tangent vector (v1, v2, v3)
class Pair{
    public:
        double x; double y; double z;  // a point (x, y, z) in R^3
        double v1; double v2; double v3; // a tangent vector (v1, v2, v3)

        Pair(double x_=0, double y_=0, double z_=0, double v1_=0, double v2_=0, double v3_=0){
            x = x_; y = y_; z = z_; 
            v1 = v1_; v2 = v2_; v3 = v3_;
        }

        void printPair(){
            std::cout << "(x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
            std::cout << "        v = (" << v1 << ", " << v2 << ", " << v3 << ")" << std::endl;
        }

        bool operator==(Pair const &B){
            bool pointEqual = (x == B.x) && (y == B.y) && (z == B.z);
            bool vectorEqual = (v1 == B.v1) && (v2 == B.v2) && (v3 == B.v3);
            return pointEqual && vectorEqual;
        }
};


// define a 2x2 matrix
class Matrix{
    public:
        // a 2x2 matrix [[p11, p12], [p21, p22]] 
        double p11, p12, 
               p21, p22;

        Matrix(double p11_=1, double p12_=0, double p21_=0, double p22_=1){
			p11 = p11_; p12 = p12_;
			p21 = p21_; p22 = p22_;
        }

        void printMatrix(){
			std::cout << "M = [[" << p11 << ", " << p12 << "]," << std::endl;
			std::cout << "     [" << p21 << ", " << p22 << "]]" << std::endl;
        }

        bool operator==(Matrix const &B){
          	return (p11 == B.p11) && (p12 == B.p12) && (p21 == B.p21) && (p22 == B.p22);
        }

        // overload multiplication for matrices
        Matrix operator*(Matrix const &B){
			Matrix C(p11 * B.p11 + p12 * B.p21, 
					 p11 * B.p12 + p12 * B.p22, 
					 p21 * B.p11 + p22 * B.p21, 
					 p21 * B.p12 + p22 * B.p22);
			return C;
        }

        double coeff0(){ return 0.5 * (sqr(p11) + sqr(p21) + sqr(p12) + sqr(p22)); }
        double coeff1(){ return 0.5 * (sqr(p11) + sqr(p21) - sqr(p12) - sqr(p22)); }
        double coeff2(){ return (p11 * p12 + p21 * p22); }

        // treating the 2x2 matrix as a linear transformation, how much does the direction (angle) get expanded?
        double expansionSq(double angle){
			// expansionSq(angle) = coeff0 + coeff1 * cos(2*angle) + coeff2 * sin(2*angle)
			return sqr(p11 * cos(angle) + p12 * sin(angle)) + sqr(p21 * cos(angle) + p22 * sin(angle));
        }

        // return the contracting direction (angle) of A in the range [0, pi]
        double contractDir(){
			double c1 = coeff1();
			double c2 = coeff2();
			
			double angle1 = abs(c1)<EPS ? PI/2 : atan(c2/c1);
			double angle2 = angle1 + PI;
			angle1 = angle1<0 ? angle1+2*PI : angle1;
			double f1 = c1*cos(angle1) + c2*sin(angle1);
			double f2 = c1*cos(angle2) + c2*sin(angle2);

			double angle = f1>f2 ? angle2 : angle1;
			return 0.5 * angle;
        }
};

class RotationMatrix: public Matrix{
  	public:
    	RotationMatrix(double t): Matrix(cos(t), -sin(t), sin(t), cos(t)){}
};

class DiagonalMatrix: public Matrix{
  	public:
    	DiagonalMatrix(double lambda): Matrix(lambda, 0, 0, 1/lambda){}
};




#endif