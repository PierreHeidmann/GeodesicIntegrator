#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

#include "DimensionDefinitions.hpp"
#include "geodesic_shooter.hpp"
#include "render.hpp"
#include "rk4.hpp"
#include "Minkowski.hpp"
#include "tensor.hpp"

int main(void)
{
    // ==========================================
    // ========== Shooting some test geod========
    // ==========================================

    // Setting up inital data

    const double center_x = 15;  // x-distance to the horizon starting point
    // const double center_y = -5.1;  // upper geodesics distance for highlightening the photon shell
    const double center_y = -12.5;  // upper geodesics distance
    const double center_z = 0.0;  // z-position
    const double start_time = 0.0;
    const double velocity_x = -1.0; // speed just the angle is important
    const double velocity_y = -1.0;
    const double velocity_z = 0.0;

    // const double shift_y = 0.005; // for highlightening the photon shell
    // const int numberofgeodesics = 10; // number of the geodesics for highlightening the photon shell

    const double shift_y = 0.25; 
    const int numberofgeodesics = 200; // number of the geodesics 

    const double lapse = -1.0;
    const bool null_geodesic = true;
    const double end_time = 100;  // total time of the geodesics
    const double dt = 0.05; // time steps 

    const Vec3 initial_data(center_x, center_y, center_z, start_time,
                            velocity_x, velocity_y, velocity_z, lapse);

    geodesic_shooter<Minkowski> pewpew;

    pewpew.shoot(initial_data, shift_y, numberofgeodesics, null_geodesic,
                 end_time, start_time, dt);

    return 0;
}
