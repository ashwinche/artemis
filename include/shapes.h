#ifndef SHAPES_H
#define SHAPES_H

#include <vector>
#include <glm/glm.hpp>
#include <SFML/Graphics.hpp>

struct ArGeometry3D{
  std::vector<glm::vec3> lines;
  std::vector<sf::Color> lines_c;
  std::vector<glm::vec3> triangles;
};

#endif