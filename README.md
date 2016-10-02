# IRIS-Segmentation

**Keywords**: MATLAB, C, Daugman Integro-Differential, IRIS Segmentation

# About

**C Based Human Eye IRIS Segmentation Algorithm based on Daugman's Itegro-Differential Operator
[![](https://raw.githubusercontent.com/ghazi94/IRIS-Segmentation/master/sample_output.jpg)](https://github.com/ghazi94/IRIS-Segmentation)

<br>

## How to use/run the codes
Image is being processed directly via MATLAB and a CSV (Matrix) of the BITMAP is generated.
The matlab file "Equivalent_MATLAB_Code/MainApplication.m" automatically creates a CSV file for the image put in this folder as "image1.jpg". So just put any image from CASIA database into this folder and rename it as "image1.jpg" to be processed by the MATLAB code.

All the BITMAP Matrix the values are from 0 to 255.
(The program can run fine on image processed values of 0 to 1, just change the threshold value in the C code (line 500, variable pixintcol) to 0.5)


Verify the ROWS and COLS in the CSV file ("1.csv") and make sure that the ROWS and COLS match those of the CSV file. For the image1 included in the folder the dimensions are 480*640 which is resized to 240*320 and stored n the CSV file. The ROWS and COLS has already been set to these values.

Wait for the program to execute and verify the result from the matlab result. It comes almost or exactly same. The result shows the iris and pupil co-ordinates and radius of both. Matlab even circles the output as white circles on the image.
