#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <sstream>
#include <random>
#include <functional>

#include "configuration_menu.h"
#include "collisionModels.h"
#include "particle.h"
#include "particle_sim.h"
#include "render_config.h"

using namespace sf;
using namespace std;

const int WIDTH = 1000;
const int HEIGHT = 1000;

const int MENU_WIDTH = 400;
const int MENU_HEIGHT = 800;

const Color posVectorColor = Color(255, 255, 0);
const Color velVectorColor = Color(255, 0, 255);
const Color accVectorColor = Color(0, 255, 255);


mt19937 mt(time(nullptr));

struct testParticle{
    float radius = 4.0;
    float mass = 50000000000000.0;
    vector<int> rgb = {255, 0, 0};
    vector<float> position = {int(WIDTH / 2), int(HEIGHT / 2)};
    vector<float> velocity = {0.0, 0.0};
    vector<float> acceleration = {0.0, 0.0};
};

struct TracePoint{
    Vector2f position; 
    vector<int> color;
    int spawnIteration; 
    bool markedForDelete = false;
};

// convert coordinate from between window and sim coordinates (invert y axis)
vector<float> convertCoords(vector<float> coords){
    return {coords[0], HEIGHT - coords[1]};
}

Vector2f convertCoords(Vector2f coords){
    return {coords.x, HEIGHT - coords.y};
}


void testToggleCallback(bool value) {
  cout << value << "\n";
}


int main(){

    float G = 0.0000000000667;
    long viscosityOfMedium = 100000000000;

    vector<Particle*> particles;
    vector<TracePoint> traces;

    bool drawVelocityVectors = true;
    bool drawAccelerationVectors = true;
    bool drawTraces = true;
    int traceLifeTime = 120;

    float vectorScaleFactor = 60.0;

    bool displayValues = true;

    bool collideWithWalls = true;
    bool isFrictionEnabled = true;

    CollisionModels collisionModel = CollisionModels::Elastic;

    // mouse position in config menu coordinates
    Vector2i mousePositionMenu; 
    // mouse position in particle menu coordinates
    vector<float> mousePosition = {-10.0, -10.0};
    vector<float> fixed_particle_position = {0.0, 0.0};
    bool spawning_particle = false;
    bool setting_velocity = false;

    testParticle particle_struct;

    Font font;
    FileInputStream fontIn;
    fontIn.open("slkscr.ttf");
    font.loadFromStream(fontIn);

    // RENDER WINDOWS 
    // =======================================================================
    RenderWindow particleWindow(VideoMode(WIDTH, HEIGHT), "Particle Sim");
    particleWindow.setFramerateLimit(60);
    RenderWindow menuWindow(VideoMode(MENU_WIDTH, MENU_HEIGHT), "Settings");
    menuWindow.setFramerateLimit(60);
    // =======================================================================

    RenderConfig renderConfig;

    bool sim_running = true;
    bool space_pressed = false;
    bool mouse_held = false;

    Text pausedIndicator; 
    pausedIndicator.setFont(font);
    pausedIndicator.setFillColor(Color::Red);
    pausedIndicator.setPosition(Vector2f(5, 5));
    pausedIndicator.setCharacterSize(30);
    pausedIndicator.setString("Paused.");

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> rgb_uniform_dist(0, 255);

    Clock renderTime;

    vector<CircleShape> particle_shapes;

    ParticleSim particleSim(G, viscosityOfMedium, collisionModel, particles, WIDTH, HEIGHT, collideWithWalls, isFrictionEnabled);
    
    // define callbacks that configuration menu will call to alter particle sim parameters
    // set gravitational constant callback
    function<void(double)> setG = particleSim.setG;
    // set viscosity callback
    function<void(double)> setViscosity = particleSim.setViscosity;

    function<void(bool)> setShowValues = renderConfig.setShowValues;
    function<void(bool)> setShowVectors = renderConfig.setShowVectors;
    function<void(bool)> setShowTraces = renderConfig.setShowTraces;

    ConfigurationMenu configurationMenu(menuWindow, font);
    configurationMenu.addTrackbar(setG, "Gravitational\nConstant", 0.0, 0.0000000005, G);
    configurationMenu.addTrackbar(setViscosity, "Viscosity", 10000000000, 1000000000000, viscosityOfMedium);

    configurationMenu.addToggle(setShowValues, "Show Values", renderConfig.getShowValues());
    configurationMenu.addToggle(setShowVectors, "Show Vectors", renderConfig.getShowVectors());
    configurationMenu.addToggle(setShowTraces, "Show Particle Traces", renderConfig.getShowTraces());

    // spawn a particle with a fixed position at first, until mouse released then unfix position 
    Particle particle = Particle(particle_struct.radius, particle_struct.mass, 0, particle_struct.rgb, convertCoords(mousePosition), particle_struct.velocity, particle_struct.acceleration);

    Particle * particle_pointer = nullptr;
    
    CircleShape particle_shape;
    CircleShape tracePointCircle;
    Text position_text;
    Text velocity_text;
    Text acceleration_text;
    stringstream info_stream;

    int renderIteration = 0;
    while(particleWindow.isOpen() && menuWindow.isOpen()){

        Time time = renderTime.getElapsedTime();
        float time_seconds = time.asSeconds();

        mousePosition = {float(Mouse::getPosition(particleWindow).x), float(Mouse::getPosition(particleWindow).y)};
        mousePositionMenu = Mouse::getPosition(menuWindow);

        // if mouse pressed within bounds of render window and render window has OS focus
        if (Mouse::isButtonPressed(Mouse::Left) && mousePosition[0] < WIDTH && mousePosition[0] >= 0 && mousePosition[1] < HEIGHT && mousePosition[1] >= 0 && particleWindow.hasFocus()){

            if (!mouse_held) {
                mouse_held = true;
                
                if (spawning_particle) {
                    spawning_particle = false;
                    setting_velocity = true;
                    fixed_particle_position = convertCoords(particle_pointer->getPosition());
                } else if (setting_velocity) {
                    setting_velocity = false;
                    particle_pointer = nullptr;
                } else {
                    spawning_particle = true;
                    // generate random particle color, ensuring it is not too dark
                    int r, g, b;
                    do {
                        r = rgb_uniform_dist(gen);
                        g = rgb_uniform_dist(gen);
                        b = rgb_uniform_dist(gen);
                    } while (r + b + g < 200);
                
                    particle.setColor({r, g, b});
                    particle_pointer = new Particle(particle);
                    particleSim.addParticle(particle_pointer);
                    particle_pointer->setPosition(convertCoords(mousePosition));
                }
            }

        } else {        // mouse is not pressed, reset mouse held bool
            mouse_held = false;

            if(setting_velocity) {
                setting_velocity = false;
            }
        }

        // if mouse pressed within bounds of menu window and menu window has OS focus
        if (Mouse::isButtonPressed(Mouse::Left) && mousePositionMenu.x < MENU_WIDTH && mousePositionMenu.x >= 0 && mousePositionMenu.y < MENU_HEIGHT && mousePositionMenu.y >= 0 && menuWindow.hasFocus()) {
            configurationMenu.evaluateMouseClick(mousePositionMenu);
        } else {
          configurationMenu.mouseReleased();
        }


        // CLOSE WINDOWS IF X PRESSED
        // ==========================================================
        Event particleWindowEvent;
        Event menuWindowEvent;

        while(particleWindow.pollEvent(particleWindowEvent)){
            
            if(particleWindowEvent.type == Event::Closed){
                particleWindow.close();
                menuWindow.close();
            }
        }

        while(menuWindow.pollEvent(menuWindowEvent)){
            
            if(menuWindowEvent.type == Event::Closed){
                particleWindow.close();
                menuWindow.close();
            }
        }
        // ==========================================================

        // UPDATE PARTICLES AND DISPLAY
        // ==========================================================

        // KEYBOARD EVENTS =========================================
        if (Keyboard::isKeyPressed(Keyboard::Space)){   // space to pause/unpause

            if (!space_pressed) {
                sim_running = !sim_running;
                space_pressed = true;
            }
        } else {
            space_pressed = false;
        }

        if (Keyboard::isKeyPressed(Keyboard::R)) {
            particleSim.reset();
            traces.clear();
        }

        particleWindow.clear();
        menuWindow.clear(Color(100, 100, 100));


        if (sim_running) {
            particleSim.updateParticles(time_seconds);
        } else {
            particleWindow.draw(pausedIndicator);
        }

        if (spawning_particle) {
            particle_pointer->setPosition(convertCoords(mousePosition));
            particle_pointer->setVelocity({0.0, 0.0});
        } else if (setting_velocity) {
            particle_pointer->setPosition(convertCoords(fixed_particle_position));
            particle_pointer->setPreviousTime(time_seconds);

            float delta_x = particle_pointer->getPosition()[0] - convertCoords(mousePosition)[0];
            float delta_y = particle_pointer->getPosition()[1] - convertCoords(mousePosition)[1];

            particle_pointer->setVelocity({delta_x / vectorScaleFactor, delta_y / vectorScaleFactor});
        } else {
            fixed_particle_position = {0.0, 0.0};
            particle_pointer = nullptr;
        }

        particles = particleSim.getParticles();
        for (int i = 0; i < particles.size(); i ++){
            Particle* particle = particles[i];

            Vector2f particle_pos_sim = Vector2f(particle->getPosition()[0], particle->getPosition()[1]);
            Vector2f particle_pos_window = convertCoords(particle_pos_sim);

            particle_pos_window = Vector2f(particle_pos_window.x - particle->getRadius(), particle_pos_window.y - particle->getRadius());

            particle_shape.setRadius(particle->getRadius());
            particle_shape.setFillColor(Color(particle->getColor()[0], particle->getColor()[1], particle->getColor()[2]));
            particle_shape.setPosition(particle_pos_window);

            particleWindow.draw(particle_shape);

            if (renderConfig.getShowValues()){
                int char_size = 15;
                int y_offset_multiplier = 2;

                info_stream.str(string());
                info_stream.precision(2);


                position_text.setFillColor(posVectorColor);
                position_text.setCharacterSize(char_size);
                position_text.setFont(font);

                info_stream << fixed << "pos: (" <<  particle->getPosition()[0] << ", " <<  particle->getPosition()[1] << ")";

                position_text.setString(info_stream.str());

                float x_offset = position_text.getLocalBounds().width / 2.0;
                position_text.setPosition(Vector2f(particle_pos_window.x - x_offset, particle_pos_window.y + y_offset_multiplier * particle->getRadius()));
                particleWindow.draw(position_text);


                velocity_text.setFillColor(velVectorColor);
                velocity_text.setCharacterSize(char_size);
                velocity_text.setFont(font);

                info_stream.str(string());
                info_stream << fixed << "vel: (" <<  particle->getVelocity()[0] << ", " <<  particle->getVelocity()[1] << ")";

                velocity_text.setString(info_stream.str());

                x_offset = velocity_text.getLocalBounds().width / 2.0;
                velocity_text.setPosition(Vector2f(particle_pos_window.x - x_offset, particle_pos_window.y + y_offset_multiplier * particle->getRadius() + char_size));
                particleWindow.draw(velocity_text);


                acceleration_text.setFillColor(accVectorColor);
                acceleration_text.setCharacterSize(15);
                acceleration_text.setFont(font);

                info_stream.str(string());
                info_stream << fixed << "acc: (" << particle->getAcceleration()[0] << ", " << particle->getAcceleration()[1] << ")";

                acceleration_text.setString(info_stream.str());

                x_offset = acceleration_text.getLocalBounds().width / 2.0;
                acceleration_text.setPosition(Vector2f(particle_pos_window.x - x_offset, particle_pos_window.y + y_offset_multiplier * particle->getRadius() + char_size * 2));
                particleWindow.draw(acceleration_text);
            }

            if (renderConfig.getShowVectors()){
                
                Vector2f velocity_vector;

                velocity_vector.x = particle_pos_sim.x + particle->getVelocity()[0] * vectorScaleFactor;
                velocity_vector.y = particle_pos_sim.y + particle->getVelocity()[1] * vectorScaleFactor;

                velocity_vector = convertCoords(velocity_vector);                
                
                Vertex line[] =
                {
                    Vertex(convertCoords(particle_pos_sim), velVectorColor),
                    Vertex(velocity_vector, velVectorColor)
                };

                particleWindow.draw(line, 2, sf::Lines);

            }
            
            if (renderConfig.getShowVectors()){
                
                Vector2f acceleration_vector;

                acceleration_vector.x = particle_pos_sim.x + particle->getAcceleration()[0] * vectorScaleFactor;
                acceleration_vector.y = particle_pos_sim.y + particle->getAcceleration()[1] * vectorScaleFactor;

                acceleration_vector = convertCoords(acceleration_vector);

                Vertex line[] =
                {
                    Vertex(convertCoords(particle_pos_sim), accVectorColor),
                    Vertex(acceleration_vector, accVectorColor)
                };

                particleWindow.draw(line, 2, sf::Lines);
            }

            TracePoint tracePoint; 
            tracePoint.spawnIteration = renderIteration;
            tracePoint.color = particle->getColor();
            tracePoint.position = Vector2f(convertCoords(particle->getPosition())[0], convertCoords(particle->getPosition())[1]);
            traces.push_back(tracePoint);
            
            auto it = traces.begin();
            while (it != traces.end()) {
                TracePoint& point = *it;

                tracePointCircle.setRadius(1);
                int alpha = (1.0 - (float(renderIteration - point.spawnIteration) / float(traceLifeTime))) * 255.0;
                if (alpha <= 1) {
                    alpha = 0;
                }
                tracePointCircle.setFillColor(Color(point.color[0], point.color[1], point.color[2], alpha));
                tracePointCircle.setPosition(point.position);

                if (renderConfig.getShowTraces()) {
                  particleWindow.draw(tracePointCircle);
                }

                if (renderIteration - point.spawnIteration >= traceLifeTime) {
                    it = traces.erase(it);
                } else {
                    ++it;
                }
            }
        }

        configurationMenu.drawUI();
        menuWindow.display();

        particleWindow.display();
        // ==========================================================


        renderIteration ++;
    }

    return 0;
}
