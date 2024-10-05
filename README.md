# Wave-Propagator

A shallow-water propagator program written in C++ that provides visual showcase of wave propagation in a Nx by Ny grid for 4 different test cases. Files can be read automatically and the grid size, step size and time step can be changed. Program is object-oriented and parallelised with OpenMP. Uses RK4 for time integration of PDEs (eqns of motion of waves). 

## Table of Contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [Features](#features)
4. [Contributing](#contributing)
5. [License](#license)
6. [Contact](#contact)

## Installation

1. Clone the repository:
    ```bash
    https://github.com/Ave04/Wave-Propagator.git
    ```

2. Navigate to the project directory:
    ```bash
    cd your-repo-name
    ```


## Usage

Enter values of dt, T, Nx, Ny and ic (test case to be used) in the commandline.
Run main.cpp
Using the Output.txt generated plot the graphs in MATLAB or any other graphing tool

