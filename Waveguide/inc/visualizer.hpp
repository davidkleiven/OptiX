#ifndef VISUALIZER_H
#define VISUALIZER_H
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <armadillo>

class Visualizer
{
public:
  Visualizer();
  ~Visualizer();
  void init();
  void fillVertexArray( const arma::mat &values );
  bool isOpen() const { return window->isOpen(); };
  bool pollEvent( sf::Event &event ) const { return window->pollEvent(event); };
  void close(){ window->close(); };
  void display(){ window->display(); };
  void clear(){ window->clear(sf::Color::Black); };
private:
  sf::RenderWindow *window{NULL};
  sf::View *view{NULL};
  sf::VertexArray *vArray{NULL};

  double width{640};
  double height{480};
  double colorMax{1.0};
  double colorMin{0.0};
  unsigned int vArrayNrow{0};
  unsigned int vArrayNcol{0};
};
#endif
