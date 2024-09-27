# Multivariable Control Systems
This repo is supposed to be about a few projects and functions regarding multivariable control systems: state-space modeling, uncertainty analysis, and robust control design. Here is also included implementations and functions from Khaki Sedigh's book about multivariable control systems.

## 1. State-Space Modeling and Uncertainty Analysis
Model scripts, simulation, and analysis of multivariable systems by state-space methods are dealt with in this section. Testing the robustness of the system under different uncertainty scenarios of both unstructured and structured uncertainties is to be done.

- **`chapter6.m`**: Unstructured uncertainty in a multivariable system is considered. For various levels of input and output uncertainty, the system response is simulated and analyzed via step responses and Bode plots. The script uses `ultidyn` to model uncertain systems and assess how uncertainty impacts system stability.
- **`chapter8.m`**: Expands uncertainty analysis to structured uncertainty. Robustness of a system with structured uncertainties is tested both without and with feedback. Develops robust stability margins and uses step responses, frequency-domain analysis for assessing system performance.

## 2. Fuzzy Logic Control
This section covers the implementation of FLC and state-space modeling. Underlying fuzzy inference systems control the system at different operating conditions, whereas all simulations are performed by using Simulink models.

- **`fuzzy_data.m`**: Define the state-space matrices \\\( A \\\\), \\\( B \\\\), \\\( C \\\\), and \\\( D \\\\) of the system which will be used by the FLC model for state-space simulations.

- **`StateSpace_fuzzy.fis`**: This is the FIS file in the development and simulation of fuzzy logic controllers through Simulink. It includes the rules and membership functions that shall be involved in the regulation of the system by the fuzzy controller.

## 3. Functions from Khaki Sedigh's Book
These functions provide some important functionalities for the analysis and design of control systems according to the book of Khaki Sedigh about multivariable systems; essential functions for robust control design, uncertainty modeling, and system analysis.

- **`MFD_Form.m`**: Computes the Minimal Realization of a system using MFD method.
- **`Realization_Controllable.m`**: Com­putes the transfer function of a system in its controllable canonical form.
- **`Realization_Observable.m`**: From a system transfer function determine the observable canonical form.

## 4. Project: Robustness of a System and Control Design
The script `Paperwork.m` presents in detail the multivariable control design for a reacting system. The script covers methods based on state-space and transferfunction and methods for robust control design.

### Major components of `Paperwork.m`:
• **System Dynamics - State-Space Modeling**: A model of the system is developed in terms of state-space equations, and time-varying behavior is simulated.
 • **Transfer-Function Derivation**: This script derives the transfer function from state-space matrices and plots system responses.
 • **Pole-Zero Analysis**: The script performs pole and zero placement analysis in order to modify system dynamics for stability.
- **Controllability and Observability**: To do controllability and observability checks with appropriate state-space matrices and Grammians to see if the system was fully controllable and observable.
- **Controller Design**: Controller designs were done: diagonal for some subsystems, robust controllers via SVD designed for disturbances rejection.
- **Singular Value and Sensitivity Analysis**: The script determines the singular values of the system's transfer function and sensitivity functions, by which the robustness to disturbances and system uncertainty are determined.
- **Hankel Matrix Realization**: Observable and controllable canonical forms are realized by the Hankel matrices that help perform the model order reduction and simplify the system modeling.

## 5. Additional System Analysis and Control Functions
The following utility functions for the analysis of multi-variable systems, testing of robustness, and the design of controllers are also included in this repository: unbyp4.m: This processes the history of system states in order to prepare the input data for the identification of the system. Order_Reduction.m: Model order reduction is carried out by the application of model reduction techniques such as balanced truncation and residualization that will reduce the order of system models while preserving dominant dynamics.
