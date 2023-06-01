// Compile with:
// $ g++ -std=c++11 main.cpp -o file_name -O2 -larmadillo
//c++ Armadillo_prueba.cpp -larmadillo
#include <iostream>
#include <armadillo>
#include <cmath>

int main()
{
                                                //    ^
  // Position of a particle                     //    |
  arma::vec Pos = {{0},                         //    | (0,1)
                   {1}};                        //    +---x-->

  // Rotation matrix 
  double phi = -3.1416/2; 
  arma::mat RotM = {{+cos(phi), -sin(phi)},
                    {+sin(phi), +cos(phi)}};

  Pos.print("Current position of the particle:");
  std::cout << "Rotating the point " << phi*180/3.1416 << " deg" << std::endl;

  Pos = RotM*Pos;

  Pos.print("New position of the particle:");   //    ^
                                                //    x (1,0)
                                                //    | 
                                                //    +------>

  return 0;
}

