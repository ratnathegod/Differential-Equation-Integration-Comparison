# Differential-Equation-Integration-Comparison
Solves ODEs in MATLAB using Euler and RK4 integration, compares them with accurate solutions, and analyzes RMSEs step-by-step.

#  Differential Equation Integration Comparison

---

## Overview  
This project compares **Euler integration** and **Runge-Kutta 4th Order (RK4)** methods for solving ordinary differential equations in MATLAB. It benchmarks both methods against accurate reference data and an exact solution, then calculates and visualizes RMSEs to show accuracy differences.

---

## Save Results  
- Final RMSEs and computed results are stored in a structured variable `answers` for easy reference and reporting.

---

## Key Features  
- Implements Euler and RK4 numerical integration  
- Compares numerical results with exact and reference solutions  
- Interpolates accurate data to align with computed grids  
- Calculates and prints RMSEs for each method and step size  
- Generates clear comparison plots for Euler, RK4, and exact solutions

---

## How to Run  
- Clone this repository  
- Place `AccurateDataSP25.mat` in the project folder  
- Open MATLAB and run the main script  
- View generated plots and check the `answers` structure for final results
