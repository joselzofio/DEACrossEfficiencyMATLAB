DEA Cross-Efficiency Toolbox for MATLAB
Version 1.0
Authors: Christian Kaps and José L. Zofío.

To use the Toolbox add the folder called "cetoolbox" to the MATLAB path.

The toolbox contains the following folders:

- cetoolbox: all the functions needed to calculate cross-efficiency are here. This is the folder you must add to the MATLAB path.
	- crossEff: The main function to call to calculate cross-efficiency scores, which uses all other functions in the folder to compute outputs. 
					Pass X,Y and a string. X is a matrix of Inputs. Y a matrix of Outputs and the string is the model you want to run.
	- crossEffClassic: 		The handler function to execute classic cross-efficiency calculations
	- crossEffGameTheory: 		The handler function to execute game theory cross-efficiency calculations
	- crossEffMultiplicativeCRS: 	The handler function to execute the multiplicative cross-efficiency calculations
	- crossEffOrdinal: 		The handler function to execute ordinal cross-efficiency calculations
	- crossEffRatio: 		The handler function to execute ratio cross-efficiency calculations
	- dea				Helper function to calculate DEA scores
	- deaout			Helper function to facilitate DEA score calculation
	- funcRatioScore		Helper function to facilitate ratio score cross-efficiency calculation
	- getDEAoptions			Helper function to facilitate DEA score calculation
	- test				Minimum working example (MWE) of how to use this toolbox, employing the example data and calling all five cross-efficiency methods.

- example data: warehouse data studied in Balk, B.M. & de Koster, M.B.M. & Kaps, C. & Zofío, J.L., 2017. "An Evaluation of Cross-Efficiency Methods, Applied to Measuring Warehouse Performance," ERIM Report Series Research in Management ERS-2017-015-LIS, Erasmus Research Institute of Management (ERIM). https://repub.eur.nl/pub/103185. The dataset, along with accompanying information on the survey methods and variable definitions, are available for downloading in spreadsheet format at https://doi.org/10.25397/eur.8279426.

How to cite:

Christian Kaps and José L. Zofío (2019). DEA Cross Efficiency Toolbox for Matlab. DOI: ??? https://github.com/joselzofio/DEACrossEfficiencyMATLAB.git

License: Code is distributed under the GNU-GPL3
http://www.gnu.org/licenses/gpl-3.0.html

