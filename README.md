# PiCLES ![Build Status](https://github.com/mochell/PiCLES.jl/actions/workflows/CI.yml/badge.svg?branch=main)
PiCLES is a fast and efficient wave model for Earth System Models, using Particle-in-Cell methods for better performance.

## Quick Start
A brief guide on how to use PiCLES.

## How to Install

### Step 1: Install Julia
1. Download Julia v1.10.2 or higher from the [official Julia website](https://julialang.org/downloads/).
2. Follow the installation instructions for your operating system.

### Step 2: Download PiCLES Repository
1. Open a terminal or command prompt.
2. Clone the PiCLES repository using Git:
   ```
   git clone https://github.com/mochell/PiCLES.git
   ```
3. Navigate into the cloned repository directory:
   ```
   cd PiCLES
   ```

### Step 3: Install Dependencies and Activate Environment
1. Start Julia by typing `julia` in your terminal within the PiCLES.jl directory.
2. Activate the project environment:
   ```julia
   using Pkg
   Pkg.activate(".")
   ```
3. Install the required dependencies:
   ```julia
   Pkg.instantiate()
   ```

You are now ready to use PiCLES for your simulations.

## How to Run Tests
To run the `T04_2D_reg_test.jl` file from the command line, follow these steps:

1. Open a terminal or command prompt.
2. Navigate to the directory where the `T04_2D_reg_test.jl` file is located.
3. Start Julia by typing `julia` in your terminal.
4. In the Julia REPL, include the `T04_2D_reg_test.jl` file:
    ```julia
    include("T04_2D_reg_test.jl")
    ```
5. The test will run and display the results in the terminal.

## How to Cite

