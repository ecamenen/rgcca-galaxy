# TUTORIAL RGCCA GALAXY-TOOL 

##### Version: 1.0

##### Author: Etienne CAMENEN

##### Key-words: 
omics, RGCCA, multi-block

##### EDAM operation: 
analysis, correlation, visualisation

##### Contact: 
arthur.tenenhaus@l2s.centralesupelec.fr

##### Short description:
Performs multi-variate analysis (PCA, CCA, PLS, RGCCA) and projects the variables and samples into a bi-dimensional space.

---

## Description

For the sake of comprehension of the use of the RGCCA package, the theoretical foundations of RGCCA and variations - 
that were previously published (Tenenhaus and Tenenhaus 2011 ; Tenenhaus and Tenenhaus 2014 ; Tenenhaus, and Frouin 2015 ;
Tenenhaus, Tenenhaus, and Groenen 2017).

We consider J data matrices *X1, ..., XJ*. Each *n × pj* data matrix *Xj = [ xj1, ..., xjpj ]* is called a block and represents 
a set of *pj* variables observed on *n* individuals.
The number and the nature of the variables may differ from one block to another, but the individuals must be the same across blocks.
We assume that all variables are centered. The objective of RGCCA is to find, for each block, a weighted composite of 
variables (called block component) *yj = Xj * aj, j=1, ... , J* (where *aj* is a column-vector with  *pj* elements) 
summarizing the relevant information between and within the blocks.
The block components are obtained such that (i) block components explain well their own block and/or (ii) 
block components that are assumed to be connected are highly correlated. 
In addition, RGCCA integrates a variable selection procedure, called SGCCA, allowing the identification of the most 
relevant features. Finally, as a component-based method, RGCCA/SGCCA can provide users with graphical 
representations to visualize the sources of variability within blocks and the amount of correlation between blocks.

## Load the inputs

In localhost:8080, in the tool-shed (left panel), select the analysis menu and the « multiOmics Toolbox » tool (Fig. 1). 
This tool is a wrapper for RGCCA analysis (for more information, see  [1], [2] or [3])

Download the pre-formatted files [here](https://github.com/BrainAndSpineInstitute/rgcca_Rpackage/tree/master/data). 
This folder includes three blocks with the same individuals (corresponding to the countries here) but different
types of variables (agriculture, industry and politic). According to Russett (1964), a high agriculture inequality
and a low industrial development lead to unstable political regime. 

Then, download them in Galaxy (with the download button in green, **Fig. 1**).

![Fig 1](tools/analysis/img/toolShed.png)

*Fig. 1 : Tool-shed of Galaxy with the “multiOmics Toolbox” emplacement*

The structure of the dataset should be seen in the history panel (**Fig. 2**).

![Fig 2](tools/analysis/img/history.png) 

*Fig. 2. History of Galaxy after downloading the three blocks from Russett data*

The block could be then used in the « multiOmics Toolbox » tool (**Fig. 3**).  
Use ```agriculture.tsv``` as dataset. Click on “Insert New dataset” to make the red panel appear and to add 
```industry.tsv``` and ```politic.tsv``` as new dataset.

![Fig 3](tools/analysis/img/tool.png) 

*Fig. 3. Graphical interface in Galaxy of “multiOmics Toolbox”. By default, an only parameter is required : an input file
to analyze. Another dataset could be added for a multi-bloc analyze (in red).*

## Customize the parsing

All the files should be delimited by the same type of separator (tabulation, by default) and by considering the first 
row as header. These parameters could be customized  by clicking “Yes” in parsing settings and selecting semicolon (**Fig 4**).

![Fig 4](tools/analysis/img/advParse.png) 

*Fig. 4. The panel of parsing settings should appear after clicking on “Yes”*


## Customize the analysis

Two other files are available in the data folder: connection.tsv and response.tsv (**Fig. 5**). 
1. The first one is the matrix representation of the connection/correlation between the blocks. A priori knowledge 
could be used to design the matrix. By default, none of the blocks was connected to each other, but they are all 
connected to a “superblock”, the concatenation matrix of the blocks. This superblock is used to visualize the data 
of all the blocks together in a common space. The design matrix is symmetric with a number of rows/columns equals to
the number of blocks (+1 with the superblock). As the other files, the type of separator should be defined in the 
parsing settings.
2. The response is a one-column file either a qualitative, or a quantitative variable or multiple columns containing
a disjunctive table. It could be used to visualize the group of responses for each sample (**Fig 7B**).

![Fig 5](tools/analysis/img/files.png) 

*Fig. 5. Supplementary files that could be used respectively to customize the connection between the blocks in the 
analysis and to visualize the group of a response in the samples plot.*

After clicking “Yes” to customize the analysis settings, the user could modify the default matrix by loading a 
customized one in the corresponding field (**Fig 6A**). In the same way, a response file could be added for visualization.

![Fig 6A](tools/analysis/img/advAn1.png) 

*Fig. 6A. The panel of analysis settings should appear after clicking on “Yes”. By adding a corresponding file, the user
could modify the connection between the blocks and color the samples according to a group of responses.*

##### For advanced users only :

The RGCCA could be customized with g, the (Scheme) function with *g(cov (Xj * aj, Xk * ak))* with 
four schemes : Horst (*g(x)=x*), centroid (*g(x)=|x|*),  factorial (*g(x)=x^2*) and its extension *g(x)^m* (**Fig 6B**). 
The shrinkage parameter for RGCCA, “tau” set to 1 by default could be varied until 0 to maximize the covariance 
between each variable of the block instead of the correlation between blocks). The algorithm could run randomly instead
of by using Singular Value Decomposition. Various other options could be disabled : blocks scaling, the bias estimator
for the variance and the use of a superblock (to visualize all the blocks together in the same plot). Finally, the 
number of component set two component could be increased until five.

![Fig 6A](tools/analysis/img/advAn2.png) 

*Fig. 6B. By scrolling down the advanced mode for the analysis, the user could disabled the automatic scaling, the use 
of a “bloc concatenation” for a better visualization, the estimator of the variance and the correlation, the shrinkage 
estimator “tau”, the number of components, the scheme function and the mode of initialization of the algorithm.*

## Customize the graphics

By default, the graphical outputs are based on the first component of the analysis on the X-axis and the second one on 
the Y-axis. By choosing more components in the advanced analysis mode, more components could be visualized (**Fig 7**). 
On the fingerprint output, a hundred of biomarkers is shown; this limitation could be customized.

![Fig 6A](tools/analysis/img/advGraph.png) 

*Fig. 7. In the graphical mode, the components visualized and the number of top biomarkers could be customized.*

## Visualize the plots

By executing the analysis (blue button at the bottom), four images should appear in the history panel.
For each axis of the block, the corresponding percent of average explained variance is indicated.
```sample_space.pdf``` is the projection of individuals coordinates in the selected component of the analysis, by default, on the
superblock (a concatenation of all the blocks) (**Fig. 5**). If a ```response``` file is loaded, each samples is colored according to
this group of responses.

![Fig 5](https://raw.githubusercontent.com/BrainAndSpineInstitute/rgcca_Rpackage/master/img/samples_space.png)
*Fig. 5 : Samples coordinates on the two first components for the superblock of the RGCCA*

```corcircle.pdf```  corresponds to the Pearson correlation between the variables of the block and the
selected components in the analysis (by default, the two first one) (**Fig. 6**).
For this plot and the next one, if the superblock is selected for each variable, their belonging to each block is 
illustrated by groups of color.


![Fig 6](https://raw.githubusercontent.com/BrainAndSpineInstitute/rgcca_Rpackage/master/img/variables_space.png)
*Fig. 6 : Correlation between each variable of each block (by using the superblock) and the two first components of RGCCA*

```fingerprint.pdf```  represents the weight of each variable ordered decreasingly on the RGCCA component selected in X-axis (1, by default) 
and for the selected block (superblock, by default) (**Fig. 7**). The top potential biomarkers are in the top. Their number is set by the graphical 
associated parameter.

![Fig 7](https://raw.githubusercontent.com/BrainAndSpineInstitute/rgcca_Rpackage/master/img/best_biomarkers.png)
*Fig. 7 : Top potential biomarkers among all the blocks with higher weight for the first component of RGCCA*

In ```ave.pdf```, the last one, the average variance explained of each variable in the selected block is ordered decreasingly (**Fig. 8**).

![Fig 8](https://raw.githubusercontent.com/BrainAndSpineInstitute/rgcca_Rpackage/master/inst/shiny/img/ave.png)
*Fig. 8 : Average variance explained in the first component of the superblock*

## References

1. Tenenhaus M, Tenenhaus A, Groenen PJF, (2017) Regularized generalized canonical correlation analysis: A framework for sequential multiblock component methods, Psychometrika, vol. 82, no. 3, 737–777
2. Tenenhaus  A. and Guillemot V. (2017): RGCCA Package. https://cran.r-project.org/web/packages/RGCCA/vignettes/vignette_RGCCA.pdf
3. Tenenhaus A, Tenenhaus M (2011) Regularized generalized canonical correlation analysis, vol. 76, pp. 257-284, Psychometrika
4. Russett, B.M. 1964. “Inequality and Instability: The Relation of Land Tenure to Politics.” World Politics 16:3: 442–54