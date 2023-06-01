#include <armadillo>
#include <iostream>

int main(){
	arma::Mat<double> a= arma::randu(2,2);
	std::cout<<a;
}
