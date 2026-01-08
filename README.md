# Welcome to Lattice Protein Game!
This project explores lattice-based protein folding using a genetic algorithm. It includes both a lattice enumerator in C++ and Python scripts for evaluating and visualizing GA performance.

## Repository Structure
- frontend/        D3.js visualization
- algorithms/      Genetic algorithm implementation
- dev/             Experimental scripts
- frontend/lattice-enum/  C++ lattice enumerator


### Installation Instructions
1. Clone the repository:
   ```bash
   git clone https://github.com/giakwon/protein.git
   cd protein/
   ```
2. Build the lattice enumerator:
   ```bash
   cd frontend/lattice-enum/
   make
   cd ../../
   ```
3. Make Python virtual environment with Pipenv:
   ```bash
   pipenv install
   ```


### How to Run Algorithms
1. Go to ```dev``` directory in a Pipenv shell:
   ```bash
   pipenv shell
   cd dev
   ```
2. To run the GA performance plot (~3 hrs):
   ```bash
   python3 test_ga.py
   ```
3. To run the maximum fitness plot (~3 hrs):
    ```bash
    python3 max_fitness.py
   ```

### How to Run Frontend
  1. From the repository root un:
     ```bash
     python3 -m http.server 8000
     ```
  4. Open your browser and go to:
     ```bash
     http://localhost:8000
     ```
