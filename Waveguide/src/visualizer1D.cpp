#include "visualizer1D.hpp"
#include <cstdlib>
#include <iostream>

using namespace std;
void Visualizer1D::setLimits( double minVal, double maxVal )
{
  min = minVal;
  max = maxVal;
}

unsigned int Visualizer1D::getY( double val ) const
{
  int indx = ( (val-min)/(max-min) )*height;
  indx = indx > height ? height-1:indx;
  indx = indx < 0 ? 0:indx;
  return height-indx;
}

void Visualizer1D::fillVertexArray( const arma::vec &vec )
{
  if ( vArray == NULL )
  {
    vArray = new sf::VertexArray( sf::Lines, 2 );
  }

  unsigned int color = rand()%256;
  sf::Color col;
  col.r = rand()%256;
  col.g = rand()%256;
  col.b = rand()%256;

  if ( vec.n_elem > width )
  {
    //unsigned int step = vec.n_elem/width;
    double step = static_cast<double>(vec.n_elem)/static_cast<double>(width);
    for ( unsigned int i=0;i<width-1;i++ )
    {
      unsigned int y1 = getY( vec( static_cast<unsigned int>(i*step)) );
      unsigned int y2 = getY( vec( static_cast<unsigned int>(i*step+step)) );
      (*vArray)[0].position = sf::Vector2f(i,y1);
      (*vArray)[1].position = sf::Vector2f(i+1,y2);
      (*vArray)[0].color = col;
      (*vArray)[1].color = col;
      window->draw( *vArray );
    }
  }
  else
  {
    unsigned int pixStep = width/vec.n_elem;
    for ( unsigned int i=0;i<vec.n_elem-1;i++ )
    {
      unsigned int y1 = getY( vec(i) );
      unsigned int y2 = getY( vec(i+1) );
      (*vArray)[0].position = sf::Vector2f(i*pixStep,y1);
      (*vArray)[1].position = sf::Vector2f(i*pixStep+pixStep,y2);
      (*vArray)[0].color = col;
      (*vArray)[1].color = col;
      window->draw( *vArray );
    }
  }
}
