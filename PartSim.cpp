#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <conio.h>
#include <ctime>
#include <string>
#include <SFML/Graphics.hpp>
using namespace std;

const double G = 5.0;
const double dt = 0.001;
const int steps = 400000;
const double box_size = 300.0;
const double radius = 1.5;

struct Vec2 {
    double x, y;

    Vec2 operator+(const Vec2& o) const { return { x + o.x, y + o.y }; }
    Vec2 operator-(const Vec2& o) const { return { x - o.x, y - o.y }; }
    Vec2 operator*(double s) const { return { x * s, y * s }; }
    Vec2& operator+=(const Vec2& o) { x += o.x; y += o.y; return *this; }
    Vec2& operator-=(const Vec2& o) { x -= o.x; y -= o.y; return *this; }
};

double norm(const Vec2& v) {
    return sqrt(v.x * v.x + v.y * v.y) + 1e-5;
}

struct Particle {
    Vec2 pos;
    Vec2 vel;
    double mass;
    bool active = true;
};

void compute_gravity(vector<Particle>& p, vector<Vec2>& forces) {
    int N = p.size();
    forces.assign(N, { 0.0, 0.0 });
    for (int i = 0; i < N; ++i) {
        if (!p[i].active) continue;
        for (int j = i + 1; j < N; ++j) {
            if (!p[j].active) continue;
            Vec2 r = p[j].pos - p[i].pos;
            double d = norm(r);
            Vec2 dir = r * (1.0 / d);

            Vec2 f = { 0.0, 0.0 };

            if (d > radius) {
                double Fg = G / (d * d);
                f += dir * Fg;
            }
            if (d <= 1.9 * radius) {
                double k = 1.0;
                double Fp = -k / (d * d);
                f += dir * Fp;
            }
            forces[i] += f;
            forces[j] -= f;
        }
    }
}

void odbicia(vector<Particle>& p) {
    int N = p.size();
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            Vec2 r = p[j].pos - p[i].pos;
            double d = norm(r);
            if (d > 1.0 * radius && d < 2 * radius) {
                Vec2 n = r * (1.0 / d);
                double v1n = p[i].vel.x * n.x + p[i].vel.y * n.y;
                double v2n = p[j].vel.x * n.x + p[j].vel.y * n.y;
                if (v1n - v2n < 0) continue;
                double a = 0.0;
                if (v1n - v2n > 10) a = 0.1;
                double m1 = p[i].mass;
                double m2 = p[j].mass;
                double v1n_new = (v1n * (m1 - m2) + 2 * m2 * v2n) / (m1 + m2);
                double v2n_new = (v2n * (m2 - m1) + 2 * m1 * v1n) / (m1 + m2);
                Vec2 dv1 = n * (v1n_new - v1n) * (1 - a);
                Vec2 dv2 = n * (v2n_new - v2n) * (1 - a);
                p[i].vel += dv1;
                p[j].vel += dv2;
            }
        }
    }
}

int main() {
    const int N_start = 300;
	const int sps = 1000; // steps per second
    vector<Particle> particles;
    mt19937 rng(time(nullptr));
    uniform_real_distribution<> pos_dist(-box_size, box_size);
    uniform_real_distribution<> vel_dist(0.0, 0.4);

    sf::Clock fpsClock;
    int frameCount = 0;

    for (int i = 0; i < N_start; ++i) {
        particles.push_back({
            {pos_dist(rng), pos_dist(rng)},
            {vel_dist(rng), vel_dist(rng)},
            1.0
            });
    }

    vector<Vec2> forces;

    unsigned int width = 1600, height = 1000;
    sf::RenderWindow* window = new sf::RenderWindow(sf::VideoMode({ width, height }), "Particle Simulation");
	window->setFramerateLimit(sps);


    sf::CircleShape circle(4.0f);
	circle.setOrigin(circle.getGeometricCenter());
	circle.setPosition({ width / 2.0f, height / 2.0f });
	circle.setFillColor(sf::Color::Blue);
    
    int step = 0;

    while (window->isOpen()) {
        while (const optional event = window->pollEvent()) {
            if (event->is<sf::Event::Closed>()) {
                window->close();
            }
        }

        for (int i = 0; i < 10; ++i) {
            odbicia(particles);
            compute_gravity(particles, forces);
            for (size_t i = 0; i < particles.size(); ++i) {
                Particle& pt = particles[i];
                pt.vel += forces[i] * (dt / pt.mass);
                pt.pos += pt.vel * dt;
            }
            step++;
        }

        /*double kinetic = 0.0;
        for (const auto& pt : particles) {
            kinetic += 0.5 * pt.mass * (pt.vel.x * pt.vel.x + pt.vel.y * pt.vel.y);
        }

        double potential = 0.0;
        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                Vec2 r = particles[j].pos - particles[i].pos;
                double d = norm(r);
                potential -= G * particles[j].mass * particles[i].mass / d;
            }
        }*/
                
        window->clear(sf::Color::Black);
                
        for (size_t i = 0; i < particles.size(); ++i) {
            circle.setPosition({
                static_cast<float>((particles[i].pos.x + box_size) * (width / (2.0 * box_size))),
                static_cast<float>((particles[i].pos.y + box_size) * (height / (2.0 * box_size)))
                });
            window->draw(circle);
        }

		frameCount++;
        if (fpsClock.getElapsedTime().asSeconds() >= 1.0f) {
            cout << "FPS: " << frameCount << endl;
            frameCount = 0;
            fpsClock.restart();
        }
        window->display();
            
    }

    
	delete window;

    return 0;
}
