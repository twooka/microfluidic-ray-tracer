# Microfluidic Ray Tracer
Simple ray-tracing simulation to calculate intensity response on a microfluidic setup.

## Hardware and Software Requirements and Setup
This software is written in C++17, whereas the Plotter and Automator are written in Python3.9.
Any computer capable of running g++, python3, and pip can run this software. No specific Hardware configuration is required.
As physics simulation is a high-performance task, using an older or less capable setup can lead to prolonged computation times.
The code runs on a single core, thus it is CPU-heavy and runs optimally on a single-threaded High-Frequency CPU setup.


## Usage
Run the executable and specify an OIL_THICKNESS parameter and a FILENAME to write the results.
```Plotter.py``` can be used to plot the results, and ```automator.py``` is a simple automation tool to run simulations varying one parameter.

## The Model
The represented model is a drop of fluid (dispersed phase) immersed in a second fluid (continuous phase).
The drop is modeled as a cylinder with two (possibly different) spherical caps.
Both fluids are contained in a tube, modeled as two concentric cylinders (also concentric with the cylindrical body of the drop).

The rays are collimated (orthogonally to the tube axis) and shot off a circular guide.
A detector the same size as the light guide is placed symmetrically.

The resulting intensity is calculated as a percentage of the original intensity.

The physical model is based on Geometrical Optics, where the effects considered are the following:
- Reflection
- Refraction

Whereas are ignored:
- Any interference effect
- Any polarization effect

Refraction angles are computed through Snell's Law [https://en.wikipedia.org/wiki/Snell%27s_law], and reflection coefficients through Fresnel's Equations [https://en.wikipedia.org/wiki/Fresnel_equations]

## The simulation
The simulation is iterative and computes for each ray a pair of rays (if this ray encounters a surface) or the detected intensity (if this ray encounters the detector).
A few parameters determine the accuracy of the simulation:
```RAY_RESOLUTION```: Determines the distance between rays in the original beam (lower is more precise)
```INTENSITY_THRESHOLD```: Determines the threshold intensity below which a ray is discarded (lower is more precise)
```MAX_DIST```: Determines the threshold distance from the detector above which a ray is discarded (higher is more precise)
```MAX_CYCLES```: Determines the maximum number of simulator cycles per ray after which the simulation is cut (higher is more precise)

### Usage notes:
```RAY_RESOLUTION``` is the parameter upon which the simulation mostly depends on. No good results have been achieved with values greater than ```0.1```.

```INTENSITY_THRESHOLD``` also influences the simulation significantly, but no measurable improvements seem to occur below ```0.001```

