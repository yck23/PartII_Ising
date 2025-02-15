## Ising Model Simulation

This project implements a Monte Carlo simulation of the Ising Model on a 2D square lattice using the Metropolis algorithm. The model is designed to study magnetization, energy fluctuations, heat capacity, and hysteresis in ferromagnetic systems.

## Features

* Simulates a 16×16 spin lattice with periodic boundary conditions.
* Uses the Planck unit system where the Boltzmann constant ($k_B$) is set to 1.
* Implements the Metropolis algorithm to determine spin flips.
* Provides multiple experimental tasks, including:
    * **Task 1:** Visualization of lattice spin configurations over time. Outputs spin configurations to `drawRand.dat` for visualization (e.g., using a plotting tool).
    * **Task 2:** Measurement of equilibrium magnetization. Calculates and outputs the average magnetization at equilibrium.
    * **Task 3:** Temperature-dependent magnetization analysis. Calculates and outputs magnetization as a function of temperature to `Task3N128.dat`.
    * **Task 4:** Calculation of heat capacity from energy fluctuations. Calculates and outputs heat capacity as a function of temperature to `Heatcapacity8.dat`.
    * **Task 6:** Hysteresis loop analysis. Simulates the system under a varying external magnetic field to generate a hysteresis loop, outputting the data to `hysteresis.dat`.
    * **Speed Test:** Performance evaluation of the simulation. Measures the execution time to assess the efficiency of the code.
    * **Susceptibility:** Computation of magnetic susceptibility as a function of temperature, outputting to `susceptibility12.dat`.

## Requirements

* C++ Compiler (e.g., g++, clang++)
* Standard C++ libraries (iostream, fstream, cmath, ctime, chrono, cstdlib)

## Output Files

The following output files will be generated:

* **drawRand.dat**: Contains spin configurations for visualization.
* **Task1eq.dat**: Contains energy evolution data for equilibrium checking.
* **Task3N128.dat**: Contains temperature-dependent magnetization data.
* **Heatcapacity8.dat**: Contains heat capacity data as a function of temperature.
* **hysteresis.dat**: Contains data for the hysteresis loop analysis.
* **speedTest.dat**: Contains performance testing data.
* **susceptibility12.dat**: Contains magnetic susceptibility data.

## Theory

### The Ising Model

The Ising model describes a lattice of spins, where each spin $s_i$ can take a value of either $+1$ (spin up) or $-1$ (spin down). The spins interact with their nearest neighbors on a 2D lattice, and the system's energy is determined by the following Hamiltonian:

$$
H = -J \sum_{\langle i,j \rangle} s_i s_j - \mu H \sum_i s_i
$$

Where:

- $J$ is the interaction energy between neighboring spins.
- $s_i$ is the spin of site $i$.
- The sum $\langle i,j \rangle$ runs over all pairs of neighboring spins.
- $\mu$ is the magnetic moment of the lattice.
- $H$ is the external magnetic field.

The system is studied using the Metropolis algorithm, which probabilistically accepts or rejects spin flips based on the change in energy. The probability of accepting a flip is given by the Boltzmann factor:

$$
P(\Delta E) = \begin{cases}
1, & \text{if } \Delta E < 0 \\
e^{-\Delta E / k_B T}, & \text{if } \Delta E \geq 0
\end{cases}
$$

Where $\Delta E$ is the change in energy due to the flip, $k_B$ is the Boltzmann constant, and $T$ is the temperature. The simulation evolves the system to equilibrium and measures macroscopic quantities like magnetization, energy, heat capacity, and susceptibility.

### Metropolis Algorithm

The Metropolis algorithm is a Monte Carlo method used to simulate physical systems. It works by iterating over all spins in the lattice and attempting to flip each spin. If the flip decreases the energy, it is always accepted. If the flip increases the energy, it is accepted with a probability $P(\Delta E) = e^{-\Delta E / k_B T}$. This probabilistic acceptance ensures that the system can reach thermal equilibrium at a given temperature.

The algorithm is used to sample configurations of the system according to the Boltzmann distribution, allowing the calculation of macroscopic properties like the average magnetization, energy, and heat capacity.

### Heat Capacity and Susceptibility

The heat capacity $C$ and susceptibility $\chi$ are two important thermodynamic quantities that can be derived from the energy and magnetization fluctuations, respectively.

The heat capacity is related to the fluctuation in energy:

$$
C = \frac{\langle E^2 \rangle - \langle E \rangle^2}{k_B T^2}
$$

The susceptibility is related to the fluctuation in magnetization:

$$
\chi = \frac{\langle M^2 \rangle - \langle M \rangle^2}{k_B T}
$$
