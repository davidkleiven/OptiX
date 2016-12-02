#include "visualizer.hpp"
#include <iostream>
#include <stdexcept>
//#include <SFML/Color.hpp>

using namespace std;

Visualizer::Visualizer(){};

Visualizer::~Visualizer()
{
  if ( view != NULL ) delete view;
  if ( window != NULL ) delete window;

  if ( vArray != NULL ) delete vArray;
  if ( tx != NULL ) delete tx;
}

void Visualizer::init( const char *windowName )
{
  window = new sf::RenderWindow(sf::VideoMode(width, height), windowName);
  window->setKeyRepeatEnabled(false);
  window->clear( sf::Color::Black );
  window->setVerticalSyncEnabled(true);

  tx = new sf::RenderTexture();

  if ( !tx->create(width, height) )
  {
    throw (runtime_error("Could not create texture!"));
  }

  // Define view
  //view = new sf::View( sf::FloatRect(width/2, height/2, width, height) );
  //window->setView(*view);
}

void Visualizer::init()
{
  init("Default Window");
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
      return; // This case is not yet implemented
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
      unsigned int rStart = row*rowStep;
      unsigned int colStart = col*colStep;
      double avgValue = average( values.submat(rStart, colStart, rStart+rowStep, colStart+colStep) );
      setColor( avgValue, color);
      (*vArray)[row*vArrayNrow+col].color = color;
    }
  }
  
  tx->draw(*vArray);
  tx->display();
  sf::Sprite sprite( tx->getTexture());
  // Draw onto screen
  window->draw( sprite );
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

double Visualizer::average( const arma::mat &mat )
{
  double avg = 0.0;
  for ( unsigned int i=0;i<mat.n_rows;i++ )
  {
    for ( unsigned int j=0;j<mat.n_cols;j++ )
    {
      avg += mat(i,j);
    }
  }
  return avg/(mat.n_rows*mat.n_cols);
}

void Visualizer::setColor( double value, sf::Color &color ) const
{
  int indx = 255*(value-colorMin)/(colorMax-colorMin);
  indx  = indx > 255 ? 255:indx;
  indx = indx < 0 ? 0:indx;

  switch ( cmap )
  {
    case Colormap_t::VIRIDIS:
      color.r = 255.0*cmaps.viridis[indx][0];
      color.g = 255.0*cmaps.viridis[indx][1];
      color.b = 255.0*cmaps.viridis[indx][2];
      break;
    case Colormap_t::GREYSCALE:
      color.r = indx;
      color.g = indx;
      color.b = indx;
      break;
  }
}
