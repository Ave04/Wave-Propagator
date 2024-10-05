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

1. Enter values of dt, T, Nx, Ny and ic (test case to be used) in the commandline.
2. Run main.cpp.
3. Using the Output.txt generated plot the graphs in MATLAB or any other graphing tool.

## Features

1. Object-oriented: Separate class created for shallowwater which includes member functions for discretisation and time integration of PDEs as well as for writing output to txt file.
2. Parallelised: OpenMP included to optimise runtime by parallelising bulk of operations.
3. Blas implementation: BLAS version of code included but not used as it does not perform as well as the parallelised loop-based approach.

## Contributions

Contributions are welcome! Please follow these steps:

1. Fork the repository.
2. Create a new branch:
   ```bash
   git checkout -b feature/your-feature
   ```
3. Commit your changes:
   ```bash
   git commit -m 'Add your feature'
   ```
4. Push to the branch
   ```bash
   git push origin feature/your-feature
   ```
6. Open a pull request.

## License
This project was licensed by Imperial College London

## Contact
1. Name: Ave
2. Email: avish.madhow@gmail.com
3. LinkedIn: https://www.linkedin.com/in/avish01/
