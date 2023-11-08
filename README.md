# Newton-Raphson-analysis
Intialization :
Start by initializing the system parameters, including the bus voltage magnitudes and angles
Choose an initial guess for voltage magnitudes and angles for all buses, typically starting with flat voltage profiles 
Formulate the Y bus Matrix:
This matrix represents the network's nodal admittances, which include both the shunt and line admittances between buses
Calculate the power mismatch:
calculate the real and reactive power mismatches
These mismatches are the differences between the actual power injections at each bus and the calculated power injections based on the initial voltage profiles
Formulate the Jacobian matrix:
The Jacobian matrix is a partial derivative of the power flow equations with respect to voltage magnitudes and angles
Calculate the elements of the Jacobian matrix
Solve the linearized equations:
Using the Jacobian matrix, solve a linear system of equations to obtain the voltage increments required to reduce the power mismatches
Update voltage profiles:
Update the voltage magnitudes and angles at each bus based on the calculated voltage increments
Repeat the power flow calculations with the updated voltage profiles
Check convergence:
Check if the system has reached convergence, which typically means that the power mismatches are small enough
If convergence is not achieved, return to step 4 and repeat the process until convergence is reached
