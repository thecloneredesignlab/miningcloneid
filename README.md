# <p>  <b>OrangeRx Simulation </b> </p>
OrangeRx is a simulation framework designed to model tumor dynamics and predict optimal drug regimens based on the ploidy composition of cell populations. It utilizes a Monte Carlo Tree Search (MCTS) algorithm to identify the most effective treatment sequences over multiple cycles.

Prerequisites
To run this simulation, ensure you have both R and Python installed. The application relies on the following packages:

R Requirements
- shiny
- deSolve
- reticulate

Python Requirements
- numpy
- scipy
- matplotlib

Getting Started
Run the Code: Open your R environment and run the script titled mcts_shiny_app.R.

Access the Interface: Click on the local URL link (e.g., http://127.0.0.1:xxxx) produced in the R console. This will open the OrangeRx simulation in your web browser.

Using the Simulation
1. Set Initial Ploidy Composition

In the application interface, enter your starting ploidy composition fractions as a list. The system interprets the list elements in the following order:

1st element: Fraction of 2.0 (diploid) cells.

2nd element: Fraction of 3.0 (triploid) cells.

3rd element: Fraction of 4.0 (tetraploid) cells.

Subsequent elements: Fractions for higher ploidy levels (5.0, etc.).

2. Generate Optimal Predictions

Click the "Run MCTS Optimal Drug Prediction" button. The simulation will run a Monte Carlo Tree Search to calculate the optimal drug regimen for the next 5 cycles.

3. Run a Specific Drug Cycle

To manually apply a treatment to the current population:

Drug Name: Type the name of the drug you wish to apply into the box labeled "Drug Name" (e.g., gemcitabine, alisertib, bay1895344, or ispinesib).

Cycle Length: Enter the number of days you want the treatment to last in the box labeled "Cycle Length (Days)".

Execute: Click "Run Next Cycle" to update the tumor state and view the results.

Model Details
The simulation uses a Hill-type kill rate function (ϕ) and pharmacokinetic (PK) models for both IV bolus and oral steady-state dosing to predict how different cell ploidies respond to specific treatments. Optimal regimens are determined by minimizing the total tumor burden while considering the varying sensitivities of cell populations of a specific ploidy.
