#ifndef PARTICLE_HEADER
#define PARTICLE_HEADER

#include <vector>

class Particle{
    public:
        Particle();    // default (random) constructor

        Particle(float, float, float, float, std::vector<int>, std::vector<float>, std::vector<float>, std::vector<float>);

        Particle(Particle&);

        std::vector<float> getPosition() const;

        void setPosition(std::vector<float>);

        float getRadius() const;

        void setRadius(float);

        float getMass() const;

        void setMass(float);

        float getSurfaceTensionStrength() const;

        float getPreviousTime();

        void setPreviousTime(float);
        
        std::vector<int> getColor() const;

        void setColor(std::vector<int>);

        std::vector<float> getVelocity() const;

        void setVelocity(std::vector<float>);

        std::vector<float> getAcceleration();

        void setAcceleration(std::vector<float>);

        bool getIsColliding() const;

        void setIsColliding(bool);

        void updateParticle(float);

    private:
        float _radius;
        float _mass;
        float _surfaceTensionStrength;
        float _previous_time;
        std::vector<int> _rgb;
        std::vector<float> _position;
        std::vector<float> _velocity;
        std::vector<float> _acceleration;
        
        bool _isColliding;
};

#endif /* !PARTICLE_HEADER */
