# Microscopic Dynamics and Wave Equation Models

This repository contains numerical implementations of microscopic models for social and physical systems, alongside solvers for wave propagation equations using Finite Volume methods.

## 1. Microscopic Models
Implementation of agent-based dynamics focusing on social interactions and transportation:
* **Opinion Dynamics**: Modeling interaction and compromise between agents to study consensus and cluster formation.
* **Vehicular Traffic**: Microscopic simulation of vehicle flow and interaction patterns, analyzing speed and density relationships.

## 2. Wave Equation
Numerical solution of the wave equation focusing on conservation laws and flux approximations.
* **Numerical Schemes**: Implementation of Finite Volume methods (FVM).
* **Solvers**: Upwind and Rusanov (Local Lax-Friedrichs) schemes for the approximation of numerical fluxes.
* **Objective**: Analysis of wave propagation, stability, and capture of discontinuities in the solution.

## Repository Structure
* **/microscopic_models**: MATLAB scripts for opinion dynamics and traffic flow simulations.
* **/wave_equation**: Source code for Finite Volume solvers (Upwind and Rusanov).
* **/docs**: Technical documentation and analysis of the mathematical models.
