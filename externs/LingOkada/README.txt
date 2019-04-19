IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING

By downloading, copying, installing or using the code, 
you agree with the following license claims:

COPYRIGHT (C) 2006, 2007, HAIBIN LING AND KAZUNORI OKADA, ALL RIGHTS RESERVED.

REDISTRIBUTION AND USE IN SOURCE AND BINARY FORMS, WITH OR WITHOUT MODIFICATION,
ARE PERMITTED FOR NON-PROFIT RESEARCH USE ONLY, PROVIDED THAT THE FOLLOWING 
CONDITIONS ARE MET:

REDISTRIBUTION'S OF SOURCE CODE MUST RETAIN THE ABOVE COPYRIGHT NOTICE,
THIS LIST OF CONDITIONS AND THE FOLLOWING DISCLAIMER.

REDISTRIBUTION'S IN BINARY FORM MUST REPRODUCE THE ABOVE COPYRIGHT NOTICE,
THIS LIST OF CONDITIONS AND THE FOLLOWING DISCLAIMER IN THE DOCUMENTATION
AND/OR OTHER MATERIALS PROVIDED WITH THE DISTRIBUTION.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

If you do not agree to this license agreement, 
do not download, install, copy or use the software. 

------------------------------------------------------------------------------------------------------
This code implements the EMD-L1 distance algorithm published in the following paper

	H. Ling and K. Okada, 
	An Efficient Earth Mover's Distance Algorithm for Robust Histogram Comparison, 
	IEEE Transaction on Pattern Analysis and Machine Intelligence (PAMI), 
	29(5):840-853, 2007.

This code is a RE-IMPLEMENTATION of the EMD-L1 distance according to the above paper. 
	There is no guarantee of exact reproduction of the experimental results reported in the paper. 
	In addition, the original implementation and the method are licensed by
 
	Siemens Corporate Research, Inc.
	755 College Road East
	Princeton, NJ 08540-6632
	Telephone: (609) 734-6500
	Fax: (609) 734-6565

	The usage of this code is restricted for non-profit research usage only and 
	using of the code is at the user's risk.

------------------------------------------------------------------------------------------------------
NOTES:
1. This implemention borrowed some basic framework from Rabner's original EMD code.

2. Histogram matrices are assumed to be arranged in the "matlab" style, i.e, the [i,j,k]-th element is 
	located in the position
		i*(n2*n3) + j*(n3) + k

3. The example usage of the code is demonstrated in 
		emdL1_test.cpp

4. This is version 3. It has been tested using VC++ 6.0, and VC 2005 (professional version).

Author Contact Information:
	Haibin Ling (hbling at temple dot edu)
	Kazunori Okada (kazokada at sfsu dot edu)
	9/18/08

