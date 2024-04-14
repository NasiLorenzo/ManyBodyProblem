#include <SFML/Graphics.hpp>
// #include <SFML/OpenGL.hpp>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>
// #include "imgui.h"
// #include "imgui-SFML.h"

struct par {
  double sigma{0.01};
  static const unsigned int dim{2};
  double reproduction{20};
  double deltaT{0.01};
  static const unsigned int n = 2;
  unsigned int size{100};
  unsigned int rate{
      1};  // rapporto tra la dimensione dello schermo e della generazione
  unsigned int rate2{100};
  double vel_factor{100000};
  std::array<unsigned int, 2> pixel{1010 * rate, 710 * rate};
  double pi = 3.141592;
  double theta{pi / 12};
  double G{6.7 * pow(10, 4)};
  double radius{2};
};
par para;

struct boidstate {
  std::array<double, para.dim> pos;
  std::array<double, para.dim> vel;
  std::array<double, para.dim> acc;
  double massa;
  boidstate(std::array<double, para.dim> Pos, std::array<double, para.dim> Vel,
      std::array<double, para.dim> Acc, double Massa)
      : pos(Pos), vel(Vel), acc(Acc), massa(Massa) {}
  // boidstate() : pos(), vel(), acc() {}
};

auto generate(std::default_random_engine
        eng) {  // genera pos e vel di un boid distribuiti secondo
  // una gauss centrata in 0
  boidstate boid({}, {}, {}, 0);
  std::normal_distribution<double> dist(0.0, para.sigma);
  for (auto it = boid.pos.begin(); it != boid.pos.end(); ++it) {
    *it = dist(eng);
  }
  for (auto it = boid.vel.begin(), last = boid.vel.end(); it != last; ++it) {
    *it = (para.vel_factor * dist(eng));
  }
  return boid;
}

double distance(const boidstate& a, const boidstate& b) {  // sqrt dispendiosa
  double s{};
  for (auto it = a.pos.begin(), index = b.pos.begin(); it != a.pos.end();
       ++it, ++index) {
    s += pow((*it) - (*index), 2);
  }
  return s;
}

using stormo = std::vector<boidstate>;

auto generator(std::default_random_engine eng) {
  stormo set;
  for (int i = 0; i < para.size; i++) {
    auto pix = para.pixel.begin();
    boidstate boidprova{generate(eng)};
    for (auto it = boidprova.pos.begin(); it != boidprova.pos.end();
         ++it, ++pix) {
      std::uniform_real_distribution<double> dis(0, *pix);
      *it += dis(eng);
    }
    set.push_back(boidprova);
  }
  return set;
}

auto neighbors(std::vector<boidstate>& set, boidstate& boid, double d) {
  stormo neighbors{};
  for (auto index = set.begin(); index != set.end(); ++index) {
    if (distance(boid, *index) < pow(d, 2) && distance(boid, *index) != 0)
      neighbors.push_back(*index);
  }
  return neighbors;
}

auto meanvel(stormo& set) {
  double s{};
  for (auto it = set.begin(); it != set.end(); ++it) {
    s += pow((*it).vel[0], 2) + pow((*it).vel[1], 2);
  }
  return sqrt(s) / set.size();
}
auto compx(stormo& set) {
  double s{};
  for (auto it = set.begin(); it != set.end(); ++it) {
    s += (*it).vel[0];
  }

  return s / set.size();
}
auto compy(stormo& set) {
  double s{};
  for (auto it = set.begin(); it != set.end(); ++it) {
    s += (*it).vel[1];
  }

  return s / set.size();
}
auto mod_vel(const boidstate& boid) {
  double sum{};
  for (auto it = boid.vel.begin(); it != boid.vel.end(); ++it) {
    sum += pow(*it, 2);
  }
  return sum;
}
auto mod_pos(const boidstate& boid) {
  double sum{};
  for (auto it = boid.pos.begin(); it != boid.pos.end(); ++it) {
    sum += pow(*it, 2);
  }
  return sqrt(sum);
}

auto rotate(boidstate& boid, const double angle) {
  auto mod{mod_vel(boid)};
  auto alpha = acos(boid.vel[0] / mod);
  std::cout << "alpha " << boid.vel[0] << "\n";
  boid.vel[0] = mod * cos(alpha - angle);
  boid.vel[1] = mod * sin(alpha - angle);

  return boid;
}

auto gravity(stormo& set, boidstate& boid) {
  boid.acc = {0, 0};
  std::cout << "Massa " << boid.massa << "\n";
  for (auto jt = set.begin(); jt != set.end(); ++jt) {
    if (distance(*jt, boid) != 0) {
      std::array<double, para.dim> verij;
      for (auto index = verij.begin(), j = (*jt).pos.begin(),
                i = boid.pos.begin();
           index != verij.end(); ++index, ++j, ++i) {
        *index = *j - *i;
      }
      auto factor = para.G * (jt->massa) / pow(distance(boid, *jt), (3.0 / 2));
      assert(distance(boid, *jt) != 0);
      for (auto index = boid.acc.begin(), i = verij.begin();
           index != boid.acc.end(); ++index, ++i) {
        if (distance(*jt, boid) < pow(para.radius, 2)) factor= -factor; 
        *index += factor * (*i);
      }
    }
  }
  return boid;
}
auto potenziale(const stormo& set) {
  double U{};
  for (auto it = set.begin(); it != set.end() - 1; ++it) {
    for (auto jt = it + 1; jt != set.end(); ++jt) {
      U += -para.G * (it->massa) * (jt->massa) / sqrt(distance(*jt, *it));
    }
  }
  return U;
}
auto cinetica(const stormo& set) {
  double K{};
  for (auto it = set.begin(); it != set.end(); ++it) {
    K += 1. / 2 * ((*it).massa) * mod_vel(*it);
  }
  return K;
}
class ensemble {
  stormo set;
  stormo newset{set};
  double zoom{1.};

 public:
  ensemble(stormo& old) : set{old} {}
  auto set_() { return set; }
  auto newset_() { return newset; }
  auto size_() { return set.size(); }
  void setZoom(double NewZoom) { zoom *= NewZoom; }
  auto getZoom() { return zoom; }
  void update() {
    newset = set;
    for (auto it = set.begin(), jt = newset.begin(); it != set.end();
         ++it, ++jt) {
      *jt = gravity(set, *jt);
      auto pix = para.pixel.begin();
      // std::cout<<"moduli vel "<<mod_vel(*jt)<<"\n";
      for (auto index = (*jt).pos.begin(), accind = (*jt).acc.begin(),
                velind = (*jt).vel.begin();
           index != (*jt).pos.end(); ++index, ++accind, ++velind, ++pix) {
        (*velind) += (*accind) * para.deltaT;
        (*index) += (*velind) * para.deltaT;
        /* (*index) = fmod(*index, *pix);
         if (*index <= 0) *index += *pix;
         assert(*index <= *pix);*/
      }
    }
    std::cout << "Energia totale " << potenziale(newset) + cinetica(newset)
              << "\n";
    std::cout << "Livello di zoom " << getZoom() << "\n";
    set = newset;
  }
};

int main() {
  std::random_device r;
  std::default_random_engine eng(r());
  boidstate pianeta1({500, 400}, {-100, 0}, {0, 0}, 100);
  boidstate pianeta2({400, 300}, {100, 50}, {0, 0}, 100);
  boidstate pianeta3({300, 200}, {0, 100}, {0, 0}, 100);
  stormo flock = {pianeta1, pianeta2, pianeta3};
  ensemble prova(flock);
  std::cout << "Dimensione generazione" << prova.size_() << "\n";
  auto y = prova.set_().size();
  /*while(true){
  prova.update();
  }*/
  // std::cout<<"Dimensione set "<<y<<"\n";
  sf::RenderWindow window(
      sf::VideoMode(para.pixel[0] / para.rate, para.pixel[1] / para.rate),
      "Boids Simulation");

  // Desired frame rate
  const sf::Time frameTime = sf::seconds(1.f / 20.f);

  sf::Clock clock;
  sf::Time accumulator = sf::Time::Zero;

  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) window.close();
    }
    if (event.type == sf::Event::MouseWheelScrolled) {
      if (event.mouseWheelScroll.delta > 0) {
        prova.setZoom(1.1);
      }
      if (event.mouseWheelScroll.delta < 0) {
        prova.setZoom(0.9);
      }
    }
    // Calculate elapsed time for this frame
    sf::Time elapsedTime = clock.restart();
    accumulator += elapsedTime;

    // Update the simulation while we have enough time accumulated
    while (accumulator >= frameTime) {
      prova.update();
      accumulator -= frameTime;
    }

    window.clear(sf::Color::White);

    // Draw boids
    for (auto& boid : prova.set_()) {
      sf::CircleShape circle(5);
      // std::cout<<prova.set_().size()<<"\n";
      circle.setFillColor(sf::Color::Black);
      circle.setPosition(boid.pos[0] / para.rate * prova.getZoom(),
          boid.pos[1] / para.rate *
              prova.getZoom());  // Assuming x and y are in pos[0]
                                 // and pos[1] respectively
      window.draw(circle);
    }

    window.display();

    // Delay to achieve desired frame rate
    sf::sleep(frameTime - clock.getElapsedTime());
  }
}
