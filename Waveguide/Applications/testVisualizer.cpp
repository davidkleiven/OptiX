#include "visualizer.hpp"
#include <armadillo>

int main( int argc, char** argv )
{
  arma::mat matrix(10000,10000);
  matrix.randu();
  Visualizer vis;
  vis.init();

  while (vis.isOpen())
  {
    sf::Event event;
    while ( vis.pollEvent(event) )
    {
      if ( event.type == sf::Event::Closed )
      {
        vis.close();
      }
    }
    vis.clear();
    vis.fillVertexArray( matrix );
    vis.display();
  }
  return 0;
}
