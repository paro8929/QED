# QED
Numerical codes for QED in 2+1 dimensions

by Paul Romatschke, August 2019 

If you are using (parts of) this code, it would be great
if you can reference the corresponding paper(s):

https://arxiv.org/pdf/1908.xxxxx.pdf

File list:

pressure-QED3.cpp [main part of QED3 pressure calculation]
QED3.h [header file]
compile.h [bash shell script for compiling executable]
QED3.exe [executable for QED3 pressure]

Directories:
LegendreData [stencils and weights for Gauss-Legendre quadrature]
Published_Results_Data_Files [Data text files matching manuscript]

---------------------------------------
Setting parameters
---------------------------------------

Modify the file "pressureQED3.cpp" and recompile

Look for "Run parameters which you might want to change" in .cpp file:

* N0, N1: resolution of the thermal sum (N0) and integration (N1). Note that depending on your choice of N1, an appropriate file "Leg-Quad-NXX.dat" must be present in the "LegendreData" subdirectory where xx=N0-1

* generate tabulated output: set "write_to_file" parameter to 1, this will initialize a loop over interaction strengths (alpha=e^2 Nf/4Pi) and output the results to a file'

* single interaction value: set "write_to_file" parameter to 0, only one value of alpha will be run (set by "alpha_default")

* IR sets the infrared regulator for the spatial integration. Results should be independent for sufficiently small IR

--------------------------------------

Compiling the code: use ./compile

Output: QED3.exe

--------------------------------------

Running the code:

Use ./QED3.exe

Depending on "write_to_file" flag, you will either see a list of alpha values that the code is running (and an associated output file name), or a single value of alpha.

---------------------------------------

Output files:

results-NTxxx-Nxxx-IRcutxxx.dat

contains 5 colums: alpha, fA, fB, fV and sum

These are: interaction strength (alpha=e^2 Nf/4/Pi), photon free energy density component A and B, photon vacuum component fV and sum (fa+fb+2*fv).
 
