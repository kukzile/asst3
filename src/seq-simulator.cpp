#include "quad-tree.h"
#include "world.h"
#include <algorithm>
#include <iostream>

// TASK 1

// NOTE: You may modify any of the contents of this file, but preserve all
// function types and names. You may add new functions if you believe they will
// be helpful.
const int QuadTreeLeafSize = 8;
class SequentialNBodySimulator : public INBodySimulator {
  std::unique_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> &particles,
                                              Vec2 bmin, Vec2 bmax) {
    // TODO: implement a function that builds and returns a quadtree containing
    // particles
    std::unique_ptr<QuadTreeNode> quadTree = std::make_unique<QuadTreeNode>();
        if (particles.size() <= QuadTreeLeafSize) {  
        quadTree->isLeaf = true;
        quadTree->particles = particles;
    } else {
        quadTree->isLeaf = false;
        Vec2 mid = (bmin + bmax) * 0.5f;
        std::vector<Particle> children[4];
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
        quadTree->children[0] = buildQuadTree(children[0], bmin, mid);
        quadTree->children[1] = buildQuadTree(children[1], Vec2(mid.x, bmin.y), Vec2(bmax.x, mid.y));
        quadTree->children[2] = buildQuadTree(children[2], Vec2(bmin.x, mid.y), Vec2(mid.x, bmax.y));
        quadTree->children[3] = buildQuadTree(children[3], mid, bmax);
    }
    return quadTree;
  }
  
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
    printf("bmin: (%f, %f)\n", bmin.x, bmin.y);
    printf("bmax: (%f, %f)\n", bmax.x, bmax.y);

    // build nodes
    quadTree->root = buildQuadTree(particles, bmin, bmax);
    if (!quadTree->checkTree()) {
      std::cout << "Your Tree has Error!" << std::endl;
    }

    return quadTree;
  }
  virtual void simulateStep(AccelerationStructure *accel,
                            std::vector<Particle> &particles,
                            std::vector<Particle> &newParticles,
                            StepParameters params) override {
    // TODO: implement sequential version of quad-tree accelerated n-body
    // simulation here, using quadTree as acceleration structure
    #pragma omp parallel for
    for (int i = 0; i < (int)particles.size(); i++) {
        auto pi = particles[i];
        Vec2 force = Vec2(0.0f, 0.0f);
        std::vector<Particle> nearParticles;
        accel->getParticles(nearParticles, pi.position, params.cullRadius);
        // accumulate attractive forces to apply to particle i
        for (size_t j = 0; j < (int)nearParticles.size(); j++) {
          if ((pi.position - nearParticles[j].position).length() < params.cullRadius)
            force += computeForce(pi, nearParticles[j], params.cullRadius);
        }
        // update particle state using the computed force
        newParticles[i] = updateParticle(pi, force, params.deltaTime);
      }
  }
};

std::unique_ptr<INBodySimulator> createSequentialNBodySimulator() {
  return std::make_unique<SequentialNBodySimulator>();
}
