# README: Medical Diagnosis Bayesian Network

## Project Overview
This project focuses on using Bayesian Networks for medical diagnosis. The primary task is to compute the parameters of a provided Bayesian Network based on healthcare datasets, some of which include missing values. The learned parameters can then be used for diagnostic purposes.

## Features
- **Bayesian Network Learning**: Parameter estimation for a given network structure.
- **Handling Missing Data**: Implements strategies to infer missing values.
- **Efficiency**: Processes datasets of over 10,000 records within 2 minutes.

## File Structure
- **Input Files**:
  - `alarm.bif`: The Bayesian Network structure with missing probability values (`-1`).
  - `records.dat`: Dataset with patient records, including potential missing values.
- **Output File**:
  - `solved_alarm.bif`: The learned Bayesian Network with computed probability values.
- **Support Files**:
  - `Format_Checker.cpp`: Validates the format of the output.
  - `startup_code.cpp`: Includes helper functions for parsing and managing the Bayesian Network.

## Usage Instructions
### Prerequisites
- A machine compatible with the GCL environment (e.g., `todi`).
- Ensure all provided files are in the working directory.

### Compilation and Execution
1. **Compilation**: Use the provided `compile.sh` script to compile the code.
   ```bash
   ./compile.sh
