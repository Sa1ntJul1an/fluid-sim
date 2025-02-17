#ifndef PARTICLE_SIM_HEADER
#define PARTICLE_SIM_HEADER

#include "particle.h"
#include <vector>
#include <algorithm>
#include <set>

class ParticleSim {
    public: 
        ParticleSim(float gravity, float interaction_radius, std::vector<Particle*> particles, int width, int height);

        ~ParticleSim();

        void addParticle(Particle*);

        void updateParticles(float);

        void setViscosity(const double);

        const double getViscosity();

        void setGravity(const double);

        const double getGravity();

        std::vector<Particle*> getParticles() const;

        void reset();

    private:
        float _gravity;
        float _interaction_radius;

        int _width;
        int _height;

        std::vector<Particle*> _particles;

        float _delta_x(const Particle*, const Particle*);
        float _delta_y(const Particle*, const Particle*);
        float _euclideanDistance(const Particle*, const Particle*);
        float _theta(const Particle*, const Particle*);
        float _xComponent(const float, const float);
        float _yComponent(const float, const float);
        std::vector<float> _calculateSurfaceTension(Particle*, Particle*);
        void _collision(Particle *, Particle *);
        std::vector<float> _calculateBuoyantFoce(const Particle*);

        // custom pair comparator:
        struct _PointerPairComparator {
            bool operator()(const std::pair<Particle*, Particle*>& a, const std::pair<Particle*, Particle*>& b) const {
                // Canonicalize the pairs by sorting the pointers
                auto canonical_a = std::minmax(a.first, a.second);
                auto canonical_b = std::minmax(b.first, b.second);
                // Compare the canonicalized pairs  
                return canonical_a < canonical_b;
            }
        };

        std::set<std::pair<Particle*, Particle*>, _PointerPairComparator> _collidingPairs;
};


#endif /* !PARTICLE_SIM_HEADER */
