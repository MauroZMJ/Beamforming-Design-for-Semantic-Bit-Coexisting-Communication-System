# Beamforming-Design-for-Semantic-Bit-Coexisting-Communication-System

This repository contains the code for the paper **"Beamforming Design for Semantic-Bit Coexisting Communication System"**. Semantic communication is naturally task oriented, which leads to a distinct performance metric with respect to traditional metric such as Shannon Channel Capacity. In this work, we take beamforming design as an example to explore the physical layer design for the semantic-bit coexisting systems. 

## Key Components

### 1. Core Algorithm
*   **[`MM_FP`](./MM_FP.m)**: This is the main contribution. It implements a **Minorization-Maximization Fractional Programming (MM-FP)** algorithm to solve the beamforming optimization problem. It iteratively updates precoding matrices to maximize semantic performance while satisfying QoS constraints for bit users.

*   **[`LP-MM-FP`](./MM_FP_naive_v1.m)**: This is a low-complexity version of MM-FP, which is motivated by the derived semi closed-form solution. 

### 2. Baseline Methods
The repository includes several standard beamforming techniques for comparison:
*   **[`MRT_beamforming.m`](./MRT_beamforming.m)**: Maximum Ratio Transmission.
*   **[`ZF_beamforming.m`](./ZF_beamforming.m)**: Zero Forcing.
*   **[`WMMSE_beamforming.m`](./WMMSE_beamforming.m)**: Weighted Minimum Mean Square Error.

### 3. Simulation & Evaluation
*   **[`beamforming_comparison_qos.m`](./beamforming_comparison_qos.m)**: A main script to run simulations. It compares the performance of the different beamforming methods (MRT, ZF, WMMSE, MM-FP) under varying QoS constraints (controlled by the `beta` parameter).

*   **[`imagenet_performance.m`](./imagenet_performance.m)**: Performance evaluation on ImageNet dataset, including the original semantic metric and the one after data fitting.

*   **[`surrogate_function_evaluation.m`](./surrogate_function_evaluation.m)**: This script is for evaluating the proposed surrogate function as in Fig. 5.


### 4. Data Generation
*   **[`generate_channel.m`](./generate_channel.m)**: Generates random channel matrices for the simulations (MIMO channel setup).

## How to Run

To reproduce the results or run the simulation, you would typically start with **[`beamforming_comparison_qos.m`](./beamforming_comparison_qos.m)**. This script sets up the system parameters (antennas `Nt`, users `B` & `T`, noise `sigma`, etc.), generates channels, and loops through the algorithms to compute and compare their performance.
