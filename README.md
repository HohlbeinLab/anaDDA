# anaDDA
README AnaDDA


How to run the software

First generate an input file by running input = Generateinputfile

In the input file all parameters can be adapted to the experimental parameters of your software. (Important is the loc error (input.sigmaerror) and the range of step sizes you want to use in the fitting: framerange (Default [1:8])

Subsequently you can fit experimental data by running Fitdata(input) which will generate a prompt to select a file or if you already have a 
D loaded in MATLAB you can run Fitdata(input,D)


Trial Data
Trialdata is available in the folder Trial data which can be selected to try if software works. 


Toolboxes required:
Statistics and Machine Learning Toolbox --> Contains MLE function
Curve Fitting Toolbox --> required for fitting relationship confined functions
Parallel Computing Toolbox --> parallel for loops can increase performance
