# TMTtargMatlab
[Matlab](https://www.mathworks.com/products/matlab.html) functions for generating targeting/inclusion lists for tandem mass tag/TMT data for use with [Xcalibur](https://www.thermofisher.com/order/catalog/product/OPTON-30487) or [Maxquant.Live](http://maxquant.live/) from [Maxquant](https://www.maxquant.org/) evidence.txt files. 

* **targXcal** : generates Xcalibur targeting lists for TMT data.
* **targMql**  : generates MaxQuant.Live targeting lists for TMT data.

**Please note:** that these functions are designed for generating targeting lists for **TMT-based** data from maxquant evidence.txt files. The mass and m/z columns in Maxquant evidence files omit the mass of the TMT reporter. These functions correct for this, but as a result are **not suitable for label-free or SILAC-based** inclusion list generation.

These functions have been tested on Windows 10, and Max OS X and Matlab versions R2017b and R2018a. Maxquant v1.6.0.16 was used for generating evidence.txt files. Output files were tested with Xcalibur v4.1.31.9 and Maxquant.live v1.0

*Note that although the **targXcal** function is used for all the examples below, the input and output arguments and use of the **targMql** function is identical, though the files generated are not cross-compatible.*

## Basic use
```matlab
% Basic function call for targXcal:
  targXcal( data , 'filename.csv');
  
% Basic function call for targMql:
  targMql( data , 'filename.csv');
```

**data** can be either a table or a structure. In most cases, a user will have imported the **evidence.txt** file generated by Maxquant into Matlab using the **readtable()** function to generate a table, e.g.
```matlab
% Import using readtable
  data = readtable('evidence.txt');
```
The **filename** can include path details, if the user does not add .csv to the end of the filename this will be done automatically.

## Function default behavior
The function has been written on the basis that after importing the evidence.txt file (using e.g. readtable), the user will further processed the data to select peptides of interest. However the function can apply a number of broad filters to the data depending on the parameters set. By default, the function removes peptides identified as contaminants or reverse database hits ('FilterRevCon', true). It also selects for unique peptide sequence/charge state combinations. If multiple entries for the same sequence/charge state are identified, the function selects the entry with the highest intensity ('Unique', true).

## Optional input arguments
* **rtWidth**    : Retention time window, +/- X minutes. **Default = 3**.
* **mzDecimals** : m/z decimal places. **Default = 5**.
* **rtDecimals** : Retention time decimal places. **Default 2**.

Name/value pair input parameters:
* **'FilterConRev'** : Remove entries identified by maxquant as contaminants, or reverse database hits. **Default: True**.
* **'PEP'**          : Remove entries not meeting this PEP threshold. **Default: 1** (no filter applied).
* **'Score'**        : Remove entries not meeting this Score threshold. **Default: 0** (no filter applied).
* **'Intensity'**    : Remove entries not meeting this Intensity threshold. **Default: 0** (no filter applied).
* **'Unique'**       : Remove entries sharing both sequence and charge state, keeping the most intense. **Default: true**.
* **'Include'**      : Include only entries matching with accession numbers matching those provided in user-provided cell array. **Default: false**.
* **'Exclude'**      : Exclude entries matching with accession numbers matching those provided in user-provided cell array. If both 'Include' and 'Exclude' options are used, 'Include' is applied first **Default: false**.

## Optional output arguments
The function can optionally output the table of targets generated that is used to create the inclusion list.csv file.
```matlab
% Output table
outputTable = targXcal( data , 'filename.csv');
```

## Example function calls for different filtering requirements
#### Applying a PEP filter
```matlab
% Create inclusion list and apply PEP filter of 0.02
  targXcal( data , 'filename.csv' , 'PEP', 0.02);
```
#### Apply a Score threshold, and select only peptides with accession numbers matching those on 'Include'
```matlab
% Create cell array with target accession numbers
  accessions = {'A6NKH3';'P05386';'P05386-2';'P05387'};

% Create inclusion list, apply Score threshold of 50, and select only peptides with accession 
% numbers matching 'accessions'/
  targXcal( data , 'filename.csv','Score',50,'Include',accessions);
```
#### Generate output table, and do not remove Reverse or Contaminant hits.
```matlab
  outputTable = targXcal( data , 'filename.csv' , 'FilterRevCon' , false);
```
#### Adjust rtWidth, mzDecimals or rtDecimals
```matlab
% There are two ways to do this:

% 1. rtWidth, mzDecimals and rtDecimals are positional arguments in that order. 
% To set these to 2,4 and 1 respectively:
  targXcal( data , 'filename.csv' , 2 , 4 , 1);
  
% 2. These can be used as name-value pairs in any order:
  targXcal( data , 'filename.csv' , 'rtWidth' , 2 , 'mzDecimals' , 4 , 'rtWidth' , 1);

```

## Using the targeting lists generated by the function
### xCalibur
The **.csv** inclusion list files generated by this function can be directly loaded into xcalibur. 

##### Step 1, open the inclusion list section
![Open Inclusion List Section](/img/xcal_incl_1.PNG)


##### Step 2 Select File -> Import, and open the .csv file you just generated
![Import inclusion list](/img/xcal_incl_2.PNG)
This should populate the inclusion list. Hit done.


##### Step 3 Turn the inclusion list on in the method properties.
![Turn inclusion list on](/img/xcal_incl_3.PNG)


### MaxQuant.Live
MaxQuant.Live does not provide a means of directly importing the targeting data. Instead, open the generated **filename.csv** in spreadsheet software, and copy/paste the columns into MaxQuant.Live.

##### Step 1, open the inclusion list file in spreadsheet software e.g. Excel, Open Office, and select and copy the columns
![Select and copy columns](/img/mql_incl_1.PNG)


##### Step 2, open the targeting app in MaxQuant.Live, use the 'Paste' button to populate the targeting list
![Paste targeting list into Maxquant.live](/img/mql_incl_2.PNG)



Created: Ed Emmott, Northeastern University, Boston, MA, USA May 2019.

Email: e.emmott@northeastern.edu
WWW: http://edemmott.co.uk
Twitter: https://twitter.com/edemmott
