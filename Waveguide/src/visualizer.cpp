#include "visualizer.hpp"
#include <iostream>
//#include <SFML/Color.hpp>

using namespace std;
Visualizer::Visualizer(){};

Visualizer::~Visualizer()
{
  if ( view != NULL ) delete view;
  delete window;

  if ( vArray != NULL ) delete vArray;
}

void Visualizer::init()
{
  window = new sf::RenderWindow(sf::VideoMode(width, height), "Waveguide");
  window->setKeyRepeatEnabled(false);
  window->clear( sf::Color::Black );

  // Define view
  //view = new sf::View( sf::FloatRect(width/2, height/2, width, height) );
  //window->setView(*view);
}

void Visualizer::fillVertexArray( const arma::mat &values )
{
  if ( vArray == NULL )
  {
    if ( values.n_rows*values.n_cols > width*height )
    {
      vArrayNrow = height;
      vArrayNcol = width;
    }
    else
    {
      vArrayNrow = values.n_rows;
      vArrayNcol = values.n_cols;
    }
    vArray = new sf::VertexArray(sf::Points, vArrayNrow*vArrayNcol);
  }

  double dr = static_cast<double>(width)/static_cast<double>(values.n_rows);
  double dc = static_cast<double>(height)/static_cast<double>(values.n_cols);
  unsigned int rowStep = values.n_rows/height;
  unsigned int colStep = values.n_cols/width;
  rowStep = rowStep == 0 ? 1:rowStep;
  colStep = colStep == 0 ? 1:colStep;

  for ( unsigned int row=0;row<vArrayNrow;row++ )
  {
    for ( unsigned int col=0;col<vArrayNcol;col++)
    {
      (*vArray)[row*vArrayNrow+col].position = sf::Vector2f(row,col);
      sf::Color color;
      int value = 255.0*(values(row*rowStep,col*colStep)-colorMin)/(colorMax-colorMin);
      if ( value > 255 )
      {
        value = 255;
      }
      else if ( value < 0 )
      {
        value = 0;
      }
      color.r = value;
      color.g = value;
      color.b = value;
      (*vArray)[row*vArrayNrow+col].color = color;
    }
  }

  // Draw onto screen
  window->draw( *vArray );
}

bool Visualizer::isOpen() const
{
  if ( window != NULL )
  {
    return window->isOpen();
  }
  return false;
}

bool Visualizer::pollEvent( sf::Event &event ) const
{
  if ( window != NULL )
  {
    return window->pollEvent( event );
  }
  return false;
}
