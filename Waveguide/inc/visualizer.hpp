#ifndef VISUALIZER_H
#define VISUALIZER_H
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <armadillo>
#include "colormaps.hpp"

/** Matrix visualizer using the SFML libray */
class Visualizer
{
public:
  /** Available colomaps */
  enum class Colormap_t {GREYSCALE, VIRIDIS};
  Visualizer();
  ~Visualizer();

  /** Initialize the window handler */
  void init( const char *windowName );

  /** Initialize the window handler */
  void init();

  /** Set values to visualize */
  void fillVertexArray( const arma::mat &values );

  /** Checks if the window is still open */
  bool isOpen() const;

  /** Checks for events */
  bool pollEvent( sf::Event &event ) const;

  /** Close the window */
  void close(){ window->close(); };

  /** Dump image to screen */
  void display(){ window->display(); };

  /** Fill the window with a black background color */
  void clear(){ window->clear(sf::Color::Black); };

  /** Set upper limit of the colorscale */
  void setColorMax( double max ){ colorMax = max; };

  /** Set color map */
  void setCmap( Colormap_t cm ){ cmap = cm; };

  /** Down sample by simple averaging */
  static double average( const arma::mat &mat );

  /** Save current scene */
  sf::Image capture(){ return window->capture(); };
protected:
  sf::RenderWindow *window{NULL};
  sf::View *view{NULL};
  sf::VertexArray *vArray{NULL};
  sf::RenderTexture *tx{NULL};

  double width{640};
  double height{480};
  double colorMax{1.0};
  double colorMin{0.0};
  unsigned int vArrayNrow{0};
  unsigned int vArrayNcol{0};
  Colormap_t cmap{Colormap_t::VIRIDIS};
  Colormaps cmaps;

  /** Set color corresponding to value */
  void setColor( double value, sf::Color &color ) const;
};
#endif
