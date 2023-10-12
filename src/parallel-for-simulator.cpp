#include "quad-tree.h"
#include "world.h"
#include <algorithm>
#include <iostream>
#include <omp.h>

// TASK 2.2

// NOTE: You may modify this class definition as you see fit, as long as the
// class name, and type of simulateStep and buildAccelerationStructure remain
// the same. You may modify any code outside this class unless otherwise
// specified.

const int QuadTreeLeafSize = 128;

class ParallelForNBodySimulator : public INBodySimulator {
public:
  // TODO: implement a function that builds and returns a quadtree containing
  // particles. You do not have to preserve this function type.
  std::unique_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> &particles,
                                              Vec2 bmin, Vec2 bmax) {

    std::unique_ptr<QuadTreeNode> quadTree = std::make_unique<QuadTreeNode>();
    if (particles.size() <= QuadTreeLeafSize) {  
        quadTree->isLeaf = true;
        quadTree->particles = particles;
    } else {
        quadTree->isLeaf = false;
        Vec2 mid = (bmin + bmax) * 0.5f;
        std::vector<std::vector<Particle>> children(4);
        for (int i = 0; i < particles.size(); i++){
          auto p = particles[i];
          int posx, posy;
          if (p.position.x < mid.x) 
            posx = 0;
          else posx = 1;
          if (p.position.y < mid.y) 
            posy = 0;
          else posy = 1;
          children[posx + 2 * posy].push_back(p);
        }
        Vec2 b_mins[] = {bmin, Vec2(mid.x, bmin.y), Vec2(bmin.x, mid.y), mid};
        Vec2 b_maxs[] = {mid, Vec2(bmax.x, mid.y), Vec2(mid.x, bmax.y), bmax};
        #pragma omp parallel for schedule(static) if (particles.size() > 400)
        for (int i = 0; i < 4; i++){
            quadTree->children[i] = buildQuadTree(children[i], b_mins[i], b_maxs[i]);
        }
    }
    return quadTree;
  }

  // Do not modify this function type.
  virtual std::unique_ptr<AccelerationStructure>
  buildAccelerationStructure(std::vector<Particle> &particles) {
    // build quad-tree
    auto quadTree = std::make_unique<QuadTree>();

    // find bounds
    Vec2 bmin(1e30f, 1e30f);
    Vec2 bmax(-1e30f, -1e30f);

    for (auto &p : particles) {
      bmin.x = fminf(bmin.x, p.position.x);
      bmin.y = fminf(bmin.y, p.position.y);
      bmax.x = fmaxf(bmax.x, p.position.x);
      bmax.y = fmaxf(bmax.y, p.position.y);
    }

    quadTree->bmin = bmin;
    quadTree->bmax = bmax;

    // build nodes
    quadTree->root = buildQuadTree(particles, bmin, bmax);
    if (!quadTree->checkTree()) {
      std::cout << "Your Tree has Error!" << std::endl;
    }

    return quadTree;
  }

  // Do not modify this function type.
  virtual void simulateStep(AccelerationStructure *accel,
                            std::vector<Particle> &particles,
                            std::vector<Particle> &newParticles,
                            StepParameters params) override {
    // TODO: implement parallel version of quad-tree accelerated n-body
    // simulation here, using quadTree as acceleration structure
    #pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < (int)particles.size(); i++) {
        auto pi = particles[i];
        Vec2 force = Vec2(0.0f, 0.0f);
        std::vector<Particle> nearParticles;
        accel->getParticles(nearParticles, pi.position, params.cullRadius);
        // accumulate attractive forces to apply to particle i
        for (size_t j = 0; j < (int)nearParticles.size(); j++) {
          if (nearParticles[j].id == pi.id) continue;
          if ((pi.position - nearParticles[j].position).length() < params.cullRadius)
            force += computeForce(pi, nearParticles[j], params.cullRadius);
        }
        // update particle state using the computed force
        newParticles[i] = updateParticle(pi, force, params.deltaTime);
      }
  }
};

// Do not modify this function type.
std::unique_ptr<INBodySimulator> createParallelForNBodySimulator() {
  return std::make_unique<ParallelForNBodySimulator>();
}
