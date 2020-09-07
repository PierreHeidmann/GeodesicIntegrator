#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>


#include "geodesic_shooter.hpp"
#include "render.hpp"
#include "rk4.hpp"
#include "schwarzschild.hpp"
#include "tensor.hpp"
#include "DimensionDefinitions.hpp"
#include "Oscilloton.hpp"


using namespace std;

int main(void)
{
    // ==========================================
    // ========== Shooting some test geod========
    // ==========================================

    // Setting up inital data

    double center_x = 15;
    double center_y = -12.5;
    double center_z = 0.0;
    double start_time = 0.0;
    double velocity_x = -1.0;
    double velocity_y = 0.0;
    double velocity_z = 0.0;
    double lapse = -1.0;
    bool null_geodesic = true;

    const Vec3 initial_data(center_x, center_y, center_z, start_time,
                            velocity_x, velocity_y, velocity_z, lapse);

    geodesic_shooter<Black_Hole> pewpew;

    pewpew.shoot(initial_data, 0.25, 100, null_geodesic);

    return 0;
}
