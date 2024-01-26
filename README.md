# Moradeke Onamusi Fluid Simulator

## Interactive Fluid Simulator

This C++ program simulates a fluid system using a particle-based approach. It utilizes the NGL library and Qt for user interface integration.
![FluidSimulator](https://github.com/NCCA/ase-assignment-Radekeo/assets/66294398/da5485fa-7494-43b7-ad1d-e67ef20e8472)

## Features

- Particle-based fluid simulation.
- Real-time visualization of fluid dynamics.
- Adjustable parameters for fluid behavior.

## Prerequisites

- C++ compiler supporting C++11 standard.
- NGL (NCCA Graphics Library)
- Qt 5 0r 6

## Installation

1. Clone or download the repository to your local machine.
2. Build the source files ( e.g using CMAKE).
3. Ensure that NGL and Qt libraries are properly linked during compilation.
4. Run the executable file generated after compilation.

## Usage

1. Upon running the program, you will see the fluid simulation window.
2. Adjust the viscosity parameter for fluid behavior as needed.
3. Interact with the simulation using the provided user interface elements.
4. Observe the real-time visualization of fluid dynamics, including particle movement and collisions.

## Contributing
Contributions are welcome! If you have suggestions for improvements, feel free to open an issue or submit a pull request.

### Approach
- Modified particle system from ase-labs
- Applied equation for Smoothed-particle hydrodynamics for making particles act more like fluids
    - Reference Article: https://andrew.gibiansky.com/blog/physics/computational-fluid-dynamics/ 
- Implemented a QT gui that takes input for different fluid properties which would affect fluid behaviours in the simulation
