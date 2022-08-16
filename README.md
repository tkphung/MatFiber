# MatFiber

MatFiber is used to analyze orientations of fiber structures in imaging data. This repository includes 

* the original `MatFiber.m` algorithm with `SampleData` analyzing muscle fiber orientations
* a new set of functions (`sfvectors` & `sfhetmetric`) applying the same algorithm with experimental data (`test_images`) analyzing actin stress fibers within cells

- [MatFiber](#matfiber)
  - [History: MatFiber version of Fiber3 algorithm for automated fiber orientation analysis](#history-matfiber-version-of-fiber3-algorithm-for-automated-fiber-orientation-analysis)
    - [`SampleData`](#sampledata)
    - [`MatFiber.m`](#matfiberm)
  - [Function Implementation](#function-implementation)
    - [`sfvectors.m`](#sfvectorsm)
    - [`sfhetmetric.m`](#sfhetmetricm)
    - [Executing program](#executing-program)
  - [Authors](#authors)
  - [Acknowledgments](#acknowledgments)

## History: MatFiber version of Fiber3 algorithm for automated fiber orientation analysis

In 1998, the Cardiac Mechanics Research Group at UCSD (CMRG) published an algorithm for automated analysis of muscle fiber orientations in histologic slides: [Karlon WJ, Covell JW, McCulloch AD, Hunter JJ, Omens JH. Automated measurement of myofiber disarray in transgenic mice with ventricular expression of ras. Anat Rec 252(4):612-625, 1998](http://www.ncbi.nlm.nih.gov/pubmed/9845212). The CMRG named their implementation Fiber3, and a version of Fiber 3 is available as part of their [freely available modeling software Continuity](http://continuity.ucsd.edu/Continuity).

Our group has used the algorithm primarily for analyzing collagen fiber orientations in picrosirius red-stained sections of myocardial scar tissue. We developed a MATLAB implementation of the algorithm named MatFiber, and verified that MatFiber gives the same results as Fiber3 when analyzing the same image using the same settings.

This repository contains the MatFiber code as well as some sample images to practice using the analysis.

Please note that the results of your analysis will depend heavily on your choice of subregion size (the algorithm produces one average orientation vector for each subregion) and your choice of threshold (the threshold excludes subregions with weaker overall image gradient information)! It is very important to test a range of subregion sizes and thresholds to establish parameters that produce an accurate analysis of your images. More information about these parameters is provided in comments within the MatFiber.m file.

If you use our MatFiber code, please cite the paper where we first employed it: [Fomovsky GM, Holmes JW. Evolution of scar structure, mechanics, and ventricular function after myocardial infarction in the rat. Am J Physiol Heart Circ Physiol 298:H221-H228, 2010](http://www.ncbi.nlm.nih.gov/pubmed/19897714).

### `SampleData`

Three sample images are provided (2 simulated fiber fields, 1 experimental data) for analysis using MatFiber. Results (visualization of fiber orientations) for image `CP15S219_SUBTRACTED.tif` using different parameters in MatFiber are included.

### `MatFiber.m`

The original MatFiber code from [(Fomovsky & Holmes 2010)](http://www.ncbi.nlm.nih.gov/pubmed/19897714) set to anlayze sample `CP15S219_SUBTRACTED.tif`.

## Function Implementation

### `sfvectors.m`

This function is a reimplementation of MatFiber. It takes in an image of fibers and determines the orientations across the field of view. The input parameters are indicated by name-value pairs (see `help sfvectors`). The function outputs a vector field based on the fiber orientations.


### `sfhetmetric.m`

This function takes in the fiber orientation vector field and calculates heterogeneity of alignment using the method described by [(Richardson & Holmes 2016)](https://pubmed.ncbi.nlm.nih.gov/27224491/). _Alignment Area Under the Curve_ is single metric describing how much the local alignment of fibers exceeds the global alignment. A higher AUC indicates more local alignment.

### Executing program

The script `Shell_SF_heterogeneity.m` includes blocks of code using the fiber analysis functions and visualizing results.

## Authors

Thien-Khoi N. Phung
[@tkphung](https://twitter.com/tkphung)


## Acknowledgments

Inspiration, code snippets, etc.
* [Fomovsky & Holmes 2010](http://www.ncbi.nlm.nih.gov/pubmed/19897714).
* [Richardson & Holmes 2016](https://pubmed.ncbi.nlm.nih.gov/27224491/)