#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include <queue>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stack>

using namespace std;
// Header file for ray tracer
#include "ray_tracer.hpp"

stack<pair<double, pair<Ray<double>, int>>> pq;

vector<pair<double, Vector3D<double>>> detected;

double OIL_THICKNESS = .05;
const double DETECTOR_DISTANCE = 2.500;

const double PTFE_THICKNESS = 1.0 / 3;
const double FLOW_PIXELS = 1 + OIL_THICKNESS;

// SETUP GEOMETRICO GOCCIOLINA
double theta_0 = M_PI / 2, theta_1 = M_PI / 2, r = FLOW_PIXELS - OIL_THICKNESS, w = 2.2;
// SETUP GEOMETRICO FLUSSO e tubicino
double r_flow = FLOW_PIXELS;
double r_ptfe = r_flow + OIL_THICKNESS + PTFE_THICKNESS;
// SETUP GEOMETRICO DETECTOR
const double DETECTOR_SIZE = 1;
const int DATA_POINTS = 300;
const double DETECTION_WINDOW = (w + r) * 1.05;

const double START_Z = -DETECTION_WINDOW - DETECTOR_SIZE, END_Z = DETECTION_WINDOW;
// SETUP REFRACTION INDECES
const double n_0 = 1.0, n_1 = 1.38, n_2 = 1.467, n_3 = 1.333;

// SETUP SIMULATION PARAMETERS
const double INTENSITY_THRESHOLD = 0.005;
const double MAX_DIST = 1e9;
const double MAX_CYCLES = 700;
const double RAY_RESOLUTION = 0.2;
const double PIPE_LENGTH = 500;
const double DEBUG_NO_DROP = false;

// CONSTRUCTION OF GEOMETRIES
double r_0 = r / sin(theta_0), r_1 = r / sin(theta_1);
auto drop_cyl = Cylinder<double>(r, w, Vector3D<double>(0, 0, 0), Vector3D<double>(0, 0, 1));
auto drop_l = SphereSector<double>(r_0, theta_0, Vector3D<double>(0, 0, -(w / 2 - r_0 * cos(theta_0))), Vector3D<double>(0, 0, -1));
auto drop_r = SphereSector<double>(r_1, theta_1, Vector3D<double>(0, 0, w / 2 - r_1 * cos(theta_1)), Vector3D<double>(0, 0, 1));
auto drop = Drop<double>(drop_cyl, drop_l, drop_r);
auto flow = Cylinder<double>(r_flow, PIPE_LENGTH *w, Vector3D<double>(0, 0, 0), Vector3D<double>(0, 0, 1));
auto ptfe = Cylinder<double>(r_ptfe, PIPE_LENGTH *w, Vector3D<double>(0, 0, 0), Vector3D<double>(0, 0, 1));
auto det = Plane<double>(Vector3D<double>(DETECTOR_DISTANCE, 0, 0), Vector3D<double>(-1, 0, 0));

// AXES STRUCTURE
//  ========================== (ptfe)  y
//       __________________  flow     O ------> z
//      /                  \          |
//     /                    \        x|
//    (         drop         )        |
//     \                    /         v
//      \__________________/
//                           flow
//  ========================== (ptfe)

// RUNS A SINGLE SIMULATION -- CAN BE PARALLELIZED
double simulate(Vector3D<double> const &DETECTOR_CORNER)
{
    // INITIALIZE PRIORITY QUEUE
    size_t total_emittance = 0;
    double steps = DETECTOR_SIZE / RAY_RESOLUTION;
    for (size_t i = 0; i < steps; i++)
    {
        for (size_t j = 0; j < steps; j++)
        {
            Vector3D<double> p = DETECTOR_CORNER + Vector3D<double>(0, i * RAY_RESOLUTION, j * RAY_RESOLUTION);
            Vector3D<double> v = Vector3D<double>(1, 0, 0);
            if (in_shape(p, DETECTOR_CORNER, DETECTOR_SIZE))
            {
                total_emittance += 1;
                pq.push({1, {Ray<double>(p, v), 0}});
            }
        }
    }
    size_t RAYS_NUMBER = pq.size();

    // ITERATIVELY SIMULATE RAYS
    int j = 0;
    while (!pq.empty() && j < MAX_CYCLES * RAYS_NUMBER)
    {
        j++;
        auto [t, r_p] = pq.top();
        pq.pop();
        auto [r, medium] = r_p;

        // DISCARD LOW INTENSITY RAYS
        if (abs(t) < INTENSITY_THRESHOLD)
            continue;
        // DISCARD STRAY RAYS
        if (norm(r.p) > MAX_DIST)
            continue;

        // INTERSECT RAY WITH GEOMETRIES
        switch (medium)
        {
        case 0:
        {
            if (r < ptfe)
            {
                auto inters = r & ptfe;
                double R = reflectivity(n_0, n_1, acos(abs(normalize(r.v) * ptfe.normal(inters))));
                auto refl = reflect(r.v, ptfe.normal(inters)) * R;
                auto refr = refract(r.v, ptfe.normal(inters), n_0, n_1) * (1. - R);
                pq.push({norm(refl), {Ray<double>(inters, refl), 0}});
                pq.push({norm(refr), {Ray<double>(inters, refr), 1}});
                break;
            }
            if (!(r < det)) // HANDLE STRAY RAYS
                break;

            // DETECT RAYS
            Vector3D<double> inters = r & det;
            detected.push_back({t, inters});
            break;
        }
        case 1:
        {
            if (r < flow)
            {
                auto inters = r & flow;
                double R = reflectivity(n_1, n_2, acos(abs(normalize(r.v) * flow.normal(inters))));
                auto refl = reflect(r.v, flow.normal(inters)) * R;
                auto refr = refract(r.v, flow.normal(inters), n_1, n_2) * (1. - R);
                pq.push({norm(refl), {Ray<double>(inters, refl), 1}});
                pq.push({norm(refr), {Ray<double>(inters, refr), 2}});
                break;
            }
            if (r < ptfe)
            {
                auto inters = r & ptfe;
                double R = reflectivity(n_1, n_0, acos(abs(normalize(r.v) * ptfe.normal(inters))));
                auto refl = reflect(r.v, ptfe.normal(inters)) * R;
                auto refr = refract(r.v, ptfe.normal(inters), n_1, n_0) * (1. - R);
                pq.push({norm(refl), {Ray<double>(inters, refl), 1}});
                pq.push({norm(refr), {Ray<double>(inters, refr), 0}});
                break;
            }
            if (!(r < det)) // HANDLE STRAY RAYS
                break;
            break;
        }
        case 2:
        {
            if (r < drop && !DEBUG_NO_DROP)
            {
                auto inters = r & drop;
                double R = reflectivity(n_2, n_3, acos(abs(normalize(r.v) * drop.normal(inters))));
                auto refl = reflect(r.v, drop.normal(inters)) * R;
                auto refr = refract(r.v, drop.normal(inters), n_2, n_3) * (1. - R);
                pq.push({norm(refl), {Ray<double>(inters, refl), 2}});
                pq.push({norm(refr), {Ray<double>(inters, refr), 3}});
                break;
            }
            if (r < flow)
            {
                auto inters = r & flow;
                double R = reflectivity(n_2, n_1, acos(abs(normalize(r.v) * flow.normal(inters))));
                auto refl = reflect(r.v, flow.normal(inters)) * R;
                auto refr = refract(r.v, flow.normal(inters), n_2, n_1) * (1. - R);
                pq.push({norm(refl), {Ray<double>(inters, refl), 2}});
                pq.push({norm(refr), {Ray<double>(inters, refr), 1}});
                break;
            }
            if (!(r < det)) // HANDLE STRAY RAYS
                break;
            break;
        }
        case 3:
        {
            if (r < drop)
            {
                auto inters = r & drop;
                double R = reflectivity(n_3, n_2, acos(abs(normalize(r.v) * drop.normal(inters))));
                auto refl = reflect(r.v, drop.normal(inters)) * R;
                auto refr = refract(r.v, drop.normal(inters), n_3, n_2) * (1. - R);
                pq.push({norm(refl), {Ray<double>(inters, refl), 3}});
                pq.push({norm(refr), {Ray<double>(inters, refr), 2}});
                break;
            }
            // THERE ARE NO STRAY RAYS TO HANDLE
        }
        }
    }

    // COMPUTE DETECTED INTENSITY (ONLY RAYS IN THE DETECTOR SHAPE ARE CONSIDERED)
    double sum = 0;
    for (auto const &[t, p] : detected)
    {
        sum += in_shape(p, DETECTOR_CORNER, DETECTOR_SIZE) ? t : 0;
    }

    pq = stack<pair<double, pair<Ray<double>, int>>>();
    detected.clear();

    return sum / total_emittance; // RETURN NORMALIZED INTENSITY
}

int main(int argc, char const *argv[])
{

    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " OIL_THICKNESS filename\n";
        return 1;
    }

    OIL_THICKNESS = atof(argv[1]);
    const char *filename = argv[2];
    cerr << "OIL_THICKNESS: " << OIL_THICKNESS << endl;

    // RESET CONSTRUCTION PARAMETERS ACCORDING TO INPUT
    // ----------------------------
    r = FLOW_PIXELS - OIL_THICKNESS;
    r_flow = FLOW_PIXELS;
    r_ptfe = r_flow + OIL_THICKNESS + PTFE_THICKNESS;
    r_0 = r / sin(theta_0), r_1 = r / sin(theta_1);
    drop_cyl = Cylinder<double>(r, w, Vector3D<double>(0, 0, 0), Vector3D<double>(0, 0, 1));
    drop_l = SphereSector<double>(r_0, theta_0, Vector3D<double>(0, 0, -(w / 2 - r_0 * cos(theta_0))), Vector3D<double>(0, 0, -1));
    drop_r = SphereSector<double>(r_1, theta_1, Vector3D<double>(0, 0, w / 2 - r_1 * cos(theta_1)), Vector3D<double>(0, 0, 1));
    drop = Drop<double>(drop_cyl, drop_l, drop_r);
    flow = Cylinder<double>(r_flow, PIPE_LENGTH * w, Vector3D<double>(0, 0, 0), Vector3D<double>(0, 0, 1));
    ptfe = Cylinder<double>(r_ptfe, PIPE_LENGTH * w, Vector3D<double>(0, 0, 0), Vector3D<double>(0, 0, 1));
    det = Plane<double>(Vector3D<double>(DETECTOR_DISTANCE, 0, 0), Vector3D<double>(-1, 0, 0));
    // ----------------------------

    // FASTER IO
    ios_base::sync_with_stdio(false);
    auto f = freopen(filename, "w", stdout);
    cerr << setprecision(2) << endl;

    for (size_t k = 0; k < DATA_POINTS / 2; ++k)
    {
        double z = START_Z + (END_Z - START_Z) / DATA_POINTS * k;
        Vector3D<double> DETECTOR_CORNER = Vector3D<double>(-DETECTOR_DISTANCE, -DETECTOR_SIZE / 2, z);

        // PRINT RESULTS
        cout << (z + DETECTOR_SIZE / 2) << " " << simulate(DETECTOR_CORNER) << "\n";

        // LOADING BAR (WORKS ON LINUX)
        cerr << "\033c";
        cerr << "simulating: " << (z - START_Z) / (END_Z - START_Z) * 100 << "%\tz: " << (z + DETECTOR_SIZE / 2) << "\n";
        cerr << "<";
        for (int i = 0; i < 100; i++)
            cerr << (i < (z - START_Z) / (END_Z - START_Z) * 100 ? "=" : ".");
        cerr << ">\n";
    }
}