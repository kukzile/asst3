#include "quad-tree.h"
#include "world.h"
#include <algorithm>
#include <iostream>
#include <omp.h>

// TASK 2.1

// NOTE: You may modify this class definition as you see fit, as long as the
// class name, and type of simulateStep and buildAccelerationStructure remain
// the same. You may modify any code outside this class unless otherwise
// specified.

const int QuadTreeLeafSize = 64;
class ParallelTaskNBodySimulator : public INBodySimulator {
public:
  std::vector<std::vector<Particle>>
  particles_get(std::vector<Particle> particles, Vec2 center) {
    std::vector<std::vector<Particle>> new_particles(4);
    for (auto &p : particles) {
      if (p.position.x < center.x && p.position.y < center.y) {
        new_particles[0].push_back(p);
      } else if (p.position.x >= center.x && p.position.y < center.y) {
        new_particles[1].push_back(p);
      } else if (p.position.x < center.x && p.position.y >= center.y) {
        new_particles[2].push_back(p);
      } else if (p.position.x >= center.x && p.position.y >= center.y) {
        new_particles[3].push_back(p);
      }
    }
    return new_particles;
  }

  std::vector<std::vector<Particle>>
  parallel_particles_get(std::vector<Particle> particles, Vec2 center) {
    std::vector<std::vector<Particle>> new_particles(4);
    new_particles[0].reserve(particles.size());
    new_particles[1].reserve(particles.size());
    new_particles[2].reserve(particles.size());
    new_particles[3].reserve(particles.size());

#pragma omp parallel shared(new_particles)
    {

      int num_threads = omp_get_num_threads();
      int thread_id = omp_get_thread_num();

      std::vector<std::vector<std::vector<Particle>>> particles_private(
          4, std::vector<std::vector<Particle>>(num_threads));

      // do blocked partitioning
      int start_num = thread_id * (particles.size() / num_threads);
      int end_num = (thread_id + 1) * (particles.size() / num_threads);
      end_num = std::min(end_num, (int)particles.size());

      for (int i = start_num; i < end_num; i++) {
        auto &p = particles[i];
        if (p.position.x < center.x && p.position.y < center.y) {
          particles_private[0][thread_id].push_back(p);
        } else if (p.position.x >= center.x && p.position.y < center.y) {
          particles_private[1][thread_id].push_back(p);
        } else if (p.position.x < center.x && p.position.y >= center.y) {
          particles_private[2][thread_id].push_back(p);
        } else if (p.position.x >= center.x && p.position.y >= center.y) {
          particles_private[3][thread_id].push_back(p);
        }
      }
#pragma omp barrier
      // if (thread_id == 0) {

      for (int i = 0; i < 4; i += 1) {
        int j = thread_id;
#pragma omp critical
        new_particles[i].insert(new_particles[i].end(),
                                particles_private[i][j].begin(),
                                particles_private[i][j].end());
      }
    }
    return new_particles;
  }

  // pointer to the root of the quad-tree
  std::unique_ptr<QuadTreeNode>
  recurse(std::vector<std::vector<Particle>> &particles, Vec2 bmin, Vec2 bmax,
          Vec2 center) {

    auto quadTree = std::make_unique<QuadTreeNode>();

    quadTree->isLeaf = false;

    quadTree->children[0] = buildQuadTree(particles[0], bmin, center);
    quadTree->children[1] = buildQuadTree(particles[1], Vec2(center.x, bmin.y),
                                          Vec2(bmax.x, center.y));
    quadTree->children[2] = buildQuadTree(particles[2], Vec2(bmin.x, center.y),
                                          Vec2(center.x, bmax.y));
    quadTree->children[3] = buildQuadTree(particles[3], center, bmax);

    return quadTree;
  }

  std::unique_ptr<QuadTreeNode>
  task_recurse(std::vector<std::vector<Particle>> &particles, Vec2 bmin,
               Vec2 bmax, Vec2 center) {

    auto quadTree = std::make_unique<QuadTreeNode>();

    quadTree->isLeaf = false;

#pragma omp parallel
#pragma omp single
#pragma omp taskgroup
    {
#pragma omp task
      quadTree->children[0] = buildQuadTree(particles[0], bmin, center);
#pragma omp task
      quadTree->children[1] = buildQuadTree(
          particles[1], Vec2(center.x, bmin.y), Vec2(bmax.x, center.y));
#pragma omp task
      quadTree->children[2] = buildQuadTree(
          particles[2], Vec2(bmin.x, center.y), Vec2(center.x, bmax.y));
#pragma omp task
      quadTree->children[3] = buildQuadTree(particles[3], center, bmax);
    }
    return quadTree;
  }

  // TODO: implement a function that builds and returns a quadtree containing
  // particles. You do not have to preserve this function type.
  std::unique_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> &particles,
                                              Vec2 bmin, Vec2 bmax) {

    // if the number of particles is less than QuadTreeLeafSize, then it is a
    // leaf node
    if (particles.size() <= QuadTreeLeafSize) {
      auto quadTree = std::make_unique<QuadTreeNode>();
      quadTree->isLeaf = true;
      quadTree->particles = particles;
      return quadTree;
    } else {
      // find the center of the box
      Vec2 center;
      center.x = (bmin.x + bmax.x) / 2.0f;
      center.y = (bmin.y + bmax.y) / 2.0f;

      std::vector<std::vector<Particle>> new_particles(4);
      if (particles.size() < 300) {
        new_particles = particles_get(particles, center);
      } else {
        // std::vector<std::vector<Particle>> new_particles(4
        new_particles = parallel_particles_get(particles, center);
      }
      std::unique_ptr<QuadTreeNode> quadTree;
      if (particles.size() < 00) {
        quadTree = recurse(new_particles, bmin, bmax, center);
      } else {
        quadTree = task_recurse(new_particles, bmin, bmax, center);
      }
      return quadTree;
    }
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
        if (nearParticles[j].id == pi.id)
          continue;
        if ((pi.position - nearParticles[j].position).length() <
            params.cullRadius)
          force += computeForce(pi, nearParticles[j], params.cullRadius);
      }
      // update particle state using the computed force
      newParticles[i] = updateParticle(pi, force, params.deltaTime);
    }
  }
};

// Do not modify this function type.
std::unique_ptr<INBodySimulator> createParallelTaskNBodySimulator() {
  return std::make_unique<ParallelTaskNBodySimulator>();
}
