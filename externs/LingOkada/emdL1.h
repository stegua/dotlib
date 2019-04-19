/*
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

Author Contact Information:
	Haibin Ling (hbling at cs dot ucla dot edu)
	Kazunori Okada (kazokada at sfsu dot edu)
	09/23/07

*******************************************************************************/
#ifndef	EMDL1_H_
#define	EMDL1_H_

#pragma warning(disable:4786)

#include <vector>

//======================================================================
//  Data structures and definition
//======================================================================

#define	EMDTYPE	double	// data type
#define EMDINF	1e10
#ifndef MAX
#define MAX(a,b)	((a)>(b)?(a):(b))
#endif

typedef struct EMDEdge * PEmdEdge;
typedef struct EMDNode * PEmdNode;

typedef struct EMDNode {
   int		pos[3];				// grid position
   EMDTYPE	d;					// initial value
   int		u;					//

   // tree maintainance
   int			iLevel;			// level in the tree, 0 means root
   PEmdNode	pParent;		// pointer to its parent
   PEmdEdge	pChild;
   PEmdEdge	pPEdge;			// point to the edge coming out from its parent
} EMDNode;


typedef struct EMDEdge {
   EMDTYPE	flow;				// initial value
   int		iDir;				// 1:outward, 0:inward

   // tree maintainance
   PEmdNode	pParent;		// point to its parent
   PEmdNode	pChild;			// the child node
   PEmdEdge	pNxt;			// next child/edge
} EMDEdge;


typedef std::vector<EMDNode>	EMDNodeArray;
typedef std::vector<EMDEdge>	EMDEdgeArray;
typedef std::vector<EMDNodeArray>	EMDNodeArray2D;
typedef std::vector<EMDEdgeArray>	EMDEdgeArray2D;
typedef std::vector<EMDTYPE>		EMDTYPEArray;
typedef std::vector<EMDTYPEArray>	EMDTYPEArray2D;

//======================================================================
//  class declarition
//======================================================================
class EmdL1 {
 public:
   EmdL1();
   ~EmdL1();

   EMDTYPE EmdDist(EMDTYPE *H1, EMDTYPE *H2, int n1, int n2, int n3=0);
   EMDTYPE EmdDist(EMDTYPE *H1, EMDTYPE *H2, int n1);	// 1-D, special and simple case
   int		SetMaxIteration(int nMaxIt)		{
      m_nMaxIt=nMaxIt;
   };

 private:

   //-- SubFunctions called in the EMD algorithm
   bool	InitMemory(int n1=0, int n2=0, int n3=0);
   bool	Initialize(EMDTYPE *H1, EMDTYPE *H2, int n1, int n2, int n3=0);
   bool	GreedySolution();
   bool	GreedySolution2();
   bool	GreedySolution3();
   void	InitBVTree();					// initialize BVTree from the initial BF solution

   int		UpdateU(PEmdNode	pRoot);
   void	UpdateSubtree(PEmdNode pRoot);
   bool 	IsOptimal();
   void	FindNewSolution();
   void	FindLoopFromEnterBV();

   EMDTYPE	CompuTotalFlow();			// Computing the total flow as the final distance

 private:
   int		m_nDim;
   int		m_n1, m_n2, m_n3;			// the histogram contains m_n1 rows and m_n2 columns
   int		m_nNBV;						// number of Non-Basic Variables (NBV)
   int		m_nMaxIt;

   /*/ tree related data structure
   // The removed code here works well for VC 2005 but not for VC6
   std::vector<std::vector<EMDNode>>	m_Nodes;		// all nodes
   std::vector<std::vector<EMDEdge>>	m_EdgesRight;	// all edges to right
   std::vector<std::vector<EMDEdge>>	m_EdgesUp;		// all edges to upward

   std::vector<std::vector<std::vector<EMDNode>>>	m_3dNodes;			// all nodes for 3D
   std::vector<std::vector<std::vector<EMDEdge>>>	m_3dEdgesRight;		// all edges to right, 3D
   std::vector<std::vector<std::vector<EMDEdge>>>	m_3dEdgesUp;		// all edges to upward, 3D
   std::vector<std::vector<std::vector<EMDEdge>>>	m_3dEdgesDeep;		// all edges to deep, 3D
   /*/
   EMDNodeArray2D	m_Nodes;		// all nodes
   EMDEdgeArray2D	m_EdgesRight;	// all edges to right
   EMDEdgeArray2D	m_EdgesUp;		// all edges to upward

   std::vector<EMDNodeArray2D>	m_3dNodes;			// all nodes for 3D
   std::vector<EMDEdgeArray2D>	m_3dEdgesRight;		// all edges to right, 3D
   std::vector<EMDEdgeArray2D>	m_3dEdgesUp;		// all edges to upward, 3D
   std::vector<EMDEdgeArray2D>	m_3dEdgesDeep;		// all edges to deep, 3D

   std::vector<PEmdEdge>	m_NBVEdges;		// pointers to all NON-BV edges
   std::vector<PEmdNode>	m_auxQueue;		// auxiliary node queue

   PEmdNode	m_pRoot;					// root of the BV Tree
   PEmdEdge	m_pEnter;					// Enter BV edge
   int			m_iEnter;					// Enter BV edge, index in m_NBVEdges
   PEmdEdge	m_pLeave;					// Leave BV edge
   int			m_nItr;						// number of iteration

   // auxiliary variables for searching a new loop
   std::vector<PEmdEdge>	m_fromLoop;
   std::vector<PEmdEdge>	m_toLoop;
   int	m_iFrom;
   int m_iTo;
};

#endif
