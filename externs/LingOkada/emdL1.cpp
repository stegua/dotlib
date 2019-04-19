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

#pragma warning(disable:4786)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "emdL1.h"

/******************************************************************************
	Construction and deconstruction
 ******************************************************************************/
EmdL1::EmdL1() {
   m_pRoot	= NULL;
   m_n1	= 0;
   m_n2	= 0;
   m_n3	= 0;
   m_nDim	= 0;
   m_nMaxIt= std::numeric_limits<int>::max();
}

EmdL1::~EmdL1() {
}


/******************************************************************************
EMDTYPE EmdDist(EMDTYPE *H1, EMDTYPE *H2, int n1, int n2, int n3)

	Earth Mover's Distance (EMD) on histograms with L1 ground distance

	Input:		Two histograms, H1, H2, each with dimension n1,n2 (and n3 for 3D)
	Output:		The EMD-L1 distance

	Note: Histgram matrices are assumed to be arranged in the "matlab" style, i.e,
		the [i,j,k]-th element of H1 is located in the position
				H1 [ i*(n2*n3) + j*(n3) + k ]
******************************************************************************/
EMDTYPE EmdL1::EmdDist(EMDTYPE *H1, EMDTYPE *H2, int n1, int n2, int n3/*=0*/) {
   // For 1D case, use simple one
   if(n2==0)
      return EmdDist(H1, H2, n1);

   // Initialization
   if(!InitMemory(n1, n2, n3)) {	// Initialize memory (reallocate if necessary)
      printf("Memory allocation failed!\n");
      return -1;
   }
   Initialize(H1,H2,n1,n2,n3);		// Initialize histgrams
   GreedySolution();				// Construct an initial Basic Feasible solution
   InitBVTree();					// Initialize BVTree

   // Iteration
   bool bOptimal	= false;
   m_nItr	= 0;
   while(!bOptimal && m_nItr<m_nMaxIt) {
      // Derive U=(u_ij) for row i and column j
      if(m_nItr==0)	UpdateSubtree(m_pRoot);
      else			UpdateSubtree(m_pEnter->pChild);

      // Optimality test
      bOptimal = IsOptimal();

      // Find new solution
      if(!bOptimal)
         FindNewSolution();

      ++m_nItr;
   }
   if (m_nItr >= m_nMaxIt)
      fprintf(stdout, "WARNING: number of iters: %d\n", m_nItr);
   // Output the total flow
   return CompuTotalFlow();
}


/******************************************************************************
EMDTYPE EmdDist(EMDTYPE *H1, EMDTYPE *H2, int n1)

	Earth Mover's Distance (EMD) on 1D histograms

	Input:		Two histograms, H1, H2, each with dimension n1
	Output:		The EMD-L1 distance
******************************************************************************/
EMDTYPE EmdL1::EmdDist(EMDTYPE *H1, EMDTYPE *H2, int n1) {
   EMDTYPE	d	= 0.;
   EMDTYPE	dk	= 0.;
   for(int k=0; k<n1; ++k)	{
      dk	+= H1[k]-H2[k];
      d	+= fabs(dk);
   }
   return d;
}


/************************************************************************
	Initialize Memory
 ************************************************************************/
bool EmdL1::InitMemory(int n1/*=0*/, int n2/*=0*/, int n3/*=0*/) {
   if(m_n1==n1 && m_n2==n2 && m_n3==n3)
      return true;

   m_n1	= n1;
   m_n2	= n2;
   m_n3	= n3;
   if(m_n1==0 || m_n2==0)		m_nDim	= 0;
   else						m_nDim	= (m_n3==0) ? 2 : 3;

   if(m_nDim==2) {
      m_Nodes.resize(m_n1);
      m_EdgesUp.resize(m_n1);
      m_EdgesRight.resize(m_n1);
      for(int i1=0; i1<m_n1; ++i1) {
         m_Nodes[i1].resize(m_n2);
         m_EdgesUp[i1].resize(m_n2);
         m_EdgesRight[i1].resize(m_n2);
      }

      m_NBVEdges.resize(m_n1*m_n2*4+2);
      m_auxQueue.resize(m_n1*m_n2+2);
      m_fromLoop.resize(m_n1*m_n2+2);
      m_toLoop.resize(m_n1*m_n2+2);
   } else if(m_nDim==3) {
      m_3dNodes.resize(m_n1);
      m_3dEdgesUp.resize(m_n1);
      m_3dEdgesRight.resize(m_n1);
      m_3dEdgesDeep.resize(m_n1);
      for(int i1=0; i1<m_n1; ++i1) {
         m_3dNodes[i1].resize(m_n2);
         m_3dEdgesUp[i1].resize(m_n2);
         m_3dEdgesRight[i1].resize(m_n2);
         m_3dEdgesDeep[i1].resize(m_n2);
         for(int i2=0; i2<m_n2; ++i2) {
            m_3dNodes[i1][i2].resize(m_n3);
            m_3dEdgesUp[i1][i2].resize(m_n3);
            m_3dEdgesRight[i1][i2].resize(m_n3);
            m_3dEdgesDeep[i1][i2].resize(m_n3);
         }
      }

      m_NBVEdges.resize(m_n1*m_n2*m_n3*6+4);
      m_auxQueue.resize(m_n1*m_n2*m_n3+4);
      m_fromLoop.resize(m_n1*m_n2*m_n3+4);
      m_toLoop.resize(m_n1*m_n2*m_n3+2);
   } else
      return false;

   return true;
}

/************************************************************************
	Initialization
 ************************************************************************/
bool EmdL1::Initialize(EMDTYPE *H1, EMDTYPE *H2, int n1, int n2, int n3/*=0*/) {
   //- Set global counters
   m_pRoot	= NULL;

   // Graph initialization
   EMDTYPE *p1	= H1;
   EMDTYPE *p2	= H2;
   if(m_nDim==2) {
      for(int c=0; c<m_n2; ++c)
         for(int r=0; r<m_n1; ++r) {
            //- initialize nodes and links
            m_Nodes[r][c].pos[0]	= r;
            m_Nodes[r][c].pos[1]	= c;
            m_Nodes[r][c].d			= *(p1++)-*(p2++);
            m_Nodes[r][c].pParent	= NULL;
            m_Nodes[r][c].pChild	= NULL;
            m_Nodes[r][c].iLevel	= -1;

            //- initialize edges
            // to the right
            m_EdgesRight[r][c].pParent	= &(m_Nodes[r][c]);
            m_EdgesRight[r][c].pChild	= &(m_Nodes[r][(c+1)%m_n2]);
            m_EdgesRight[r][c].flow	= 0;
            m_EdgesRight[r][c].iDir	= 1;
            m_EdgesRight[r][c].pNxt	= NULL;

            // to the upward
            m_EdgesUp[r][c].pParent	= &(m_Nodes[r][c]);
            m_EdgesUp[r][c].pChild	= &(m_Nodes[(r+1)%m_n1][c]);
            m_EdgesUp[r][c].flow	= 0;
            m_EdgesUp[r][c].iDir	= 1;
            m_EdgesUp[r][c].pNxt	= NULL;
         }
   } else if(m_nDim==3) {
      for(int z=0; z<m_n3; ++z)
         for(int c=0; c<m_n2; ++c)
            for(int r=0; r<m_n1; ++r) {
               //- initialize nodes and edges
               m_3dNodes[r][c][z].pos[0]	= r;
               m_3dNodes[r][c][z].pos[1]	= c;
               m_3dNodes[r][c][z].pos[2]	= z;
               m_3dNodes[r][c][z].d		= *(p1++)-*(p2++);
               m_3dNodes[r][c][z].pParent	= NULL;
               m_3dNodes[r][c][z].pChild	= NULL;
               m_3dNodes[r][c][z].iLevel	= -1;

               //- initialize edges
               // to the upward
               m_3dEdgesUp[r][c][z].pParent= &(m_3dNodes[r][c][z]);
               m_3dEdgesUp[r][c][z].pChild	= &(m_3dNodes[(r+1)%m_n1][c][z]);
               m_3dEdgesUp[r][c][z].flow	= 0;
               m_3dEdgesUp[r][c][z].iDir	= 1;
               m_3dEdgesUp[r][c][z].pNxt	= NULL;

               // to the right
               m_3dEdgesRight[r][c][z].pParent	= &(m_3dNodes[r][c][z]);
               m_3dEdgesRight[r][c][z].pChild	= &(m_3dNodes[r][(c+1)%m_n2][z]);
               m_3dEdgesRight[r][c][z].flow	= 0;
               m_3dEdgesRight[r][c][z].iDir	= 1;
               m_3dEdgesRight[r][c][z].pNxt	= NULL;

               // to the deep
               m_3dEdgesDeep[r][c][z].pParent	= &(m_3dNodes[r][c][z]);
               m_3dEdgesDeep[r][c][z].pChild	= &(m_3dNodes[r][c])[(z+1)%m_n3];
               m_3dEdgesDeep[r][c][z].flow		= 0;
               m_3dEdgesDeep[r][c][z].iDir		= 1;
               m_3dEdgesDeep[r][c][z].pNxt		= NULL;
            }
   }
   return true;
}


/************************************************************************
  Greedy method to find an initial solution, analogous to Russel method
 ************************************************************************/
bool EmdL1::GreedySolution() {
   return m_nDim==2 ? GreedySolution2() : GreedySolution3();
}

bool EmdL1::GreedySolution2() {
   //- Prepare auxiliary array, D=H1-H2
   int		c,r;
   EMDTYPEArray2D D(m_n1);
   for(r=0; r<m_n1; ++r)	{
      D[r].resize(m_n2);
      for(c=0; c<m_n2; ++c)		D[r][c]	= m_Nodes[r][c].d;
   }


   // compute integrated values along each dimension
   std::vector<EMDTYPE>	d2s(m_n2);
   d2s[0]	= 0;
   for(c=0; c<m_n2-1; ++c)	{
      d2s[c+1]	= d2s[c];
      for(r=0; r<m_n1; ++r)	d2s[c+1]-= D[r][c];
   }

   std::vector<EMDTYPE>	d1s(m_n1);
   d1s[0]	= 0;
   for(r=0; r<m_n1-1; ++r)	{
      d1s[r+1]	= d1s[r];
      for(c=0; c<m_n2; ++c)		d1s[r+1]-= D[r][c];
   }

   //- Greedy algorithm for initial solution
   PEmdEdge	pBV;
   EMDTYPE	dFlow;
   bool	bUpward	= false;
   m_nNBV	= 0;		// number of NON-BV edges

   for(c=0; c<m_n2-1; ++c)
      for(r=0; r<m_n1; ++r) {
         dFlow	= D[r][c];
         bUpward	= (r<m_n1-1) && (fabs(dFlow+d2s[c+1]) > fabs(dFlow+d1s[r+1]));	// Move upward or right

         // modify basic variables, record BV and related values
         if(bUpward)	{
            // move to up
            pBV	= &(m_EdgesUp[r][c]);
            m_NBVEdges[m_nNBV++]	= &(m_EdgesRight[r][c]);

            D[r+1][c]	+= dFlow;		// auxilary matrix maintanence
            d1s[r+1]	+= dFlow;		// auxilary matrix maintanence
         } else {
            // move to right, no other choice
            pBV	= &(m_EdgesRight[r][c]);
            if(r<m_n1-1)
               m_NBVEdges[m_nNBV++]	= &(m_EdgesUp[r][c]);

            D[r][c+1]	+= dFlow;		// auxilary matrix maintanence
            d2s[c+1]	+= dFlow;		// auxilary matrix maintanence
         }

         pBV->pParent->pChild= pBV;
         pBV->flow	= fabs(dFlow);
         pBV->iDir	= dFlow>0;		// 1:outward, 0:inward
      }

   //- rightmost column, no choice but move upward
   c	= m_n2-1;
   for(r=0; r<m_n1-1; ++r) {
      dFlow		= D[r][c];
      pBV			= &(m_EdgesUp[r][c]);

      D[r+1][c]	+= dFlow;		// auxilary matrix maintanence

      pBV->pParent->pChild= pBV;
      pBV->flow	= fabs(dFlow);
      pBV->iDir	= dFlow>0;		// 1:outward, 0:inward
   }

   return true;
}


bool EmdL1::GreedySolution3() {
   //- Prepare auxiliary array, D=H1-H2
   int		i1,i2,i3;
   std::vector<EMDTYPEArray2D> D(m_n1);
   for(i1=0; i1<m_n1; ++i1) {
      D[i1].resize(m_n2);
      for(i2=0; i2<m_n2; ++i2) {
         D[i1][i2].resize(m_n3);
         for(i3=0; i3<m_n3; ++i3)
            D[i1][i2][i3]	= m_3dNodes[i1][i2][i3].d;
      }
   }

   // compute integrated values along each dimension
   std::vector<EMDTYPE>	d1s(m_n1);
   d1s[0]	= 0;
   for(i1=0; i1<m_n1-1; ++i1)	{
      d1s[i1+1]	= d1s[i1];
      for(i2=0; i2<m_n2; ++i2)
         for(i3=0; i3<m_n3; ++i3)
            d1s[i1+1]	-= D[i1][i2][i3];
   }

   std::vector<EMDTYPE>	d2s(m_n2);
   d2s[0]	= 0;
   for(i2=0; i2<m_n2-1; ++i2)	{
      d2s[i2+1]	= d2s[i2];
      for(i1=0; i1<m_n1; ++i1)
         for(i3=0; i3<m_n3; ++i3)
            d2s[i2+1]	-= D[i1][i2][i3];
   }

   std::vector<EMDTYPE>	d3s(m_n3);
   d3s[0]	= 0;
   for(i3=0; i3<m_n3-1; ++i3)	{
      d3s[i3+1]	= d3s[i3];
      for(i1=0; i1<m_n1; ++i1)
         for(i2=0; i2<m_n2; ++i2)
            d3s[i3+1]	-= D[i1][i2][i3];
   }

   //- Greedy algorithm for initial solution
   PEmdEdge	pBV;
   EMDTYPE	dFlow, f1,f2,f3;
   m_nNBV	= 0;		// number of NON-BV edges

   for(i3=0; i3<m_n3; ++i3)
      for(i2=0; i2<m_n2; ++i2)
         for(i1=0; i1<m_n1; ++i1) {
            if(i3==m_n3-1 && i2==m_n2-1 && i1==m_n1-1)	break;

            //- determine which direction to move, either right or upward
            dFlow		= D[i1][i2][i3];
            f1	= i1<m_n1-1 ? fabs(dFlow+d1s[i1+1]) : EMDINF;
            f2	= i2<m_n2-1 ? fabs(dFlow+d2s[i2+1]) : EMDINF;
            f3	= i3<m_n3-1 ? fabs(dFlow+d3s[i3+1]) : EMDINF;

            if(f1<f2 && f1<f3)	{
               pBV	= &(m_3dEdgesUp[i1][i2][i3]);				// up
               if(i2<m_n2-1)	m_NBVEdges[m_nNBV++]	= &(m_3dEdgesRight[i1][i2][i3]);	// right
               if(i3<m_n3-1)	m_NBVEdges[m_nNBV++]	= &(m_3dEdgesDeep[i1][i2][i3]);		// deep
               D[i1+1][i2][i3]	+= dFlow;						// maintain auxilary matrix
               d1s[i1+1]		+= dFlow;
            } else if(f2<f3)		{
               pBV	= &(m_3dEdgesRight[i1][i2][i3]);			// right
               if(i1<m_n1-1)	m_NBVEdges[m_nNBV++]	= &(m_3dEdgesUp[i1][i2][i3]);		// up
               if(i3<m_n3-1)	m_NBVEdges[m_nNBV++]	= &(m_3dEdgesDeep[i1][i2][i3]);		// deep
               D[i1][i2+1][i3]	+= dFlow;						// maintain auxilary matrix
               d2s[i2+1]		+= dFlow;
            } else				{
               pBV	= &(m_3dEdgesDeep[i1][i2][i3]);				// deep
               if(i2<m_n2-1)	m_NBVEdges[m_nNBV++]	= &(m_3dEdgesRight[i1][i2][i3]);	// right
               if(i1<m_n1-1)	m_NBVEdges[m_nNBV++]	= &(m_3dEdgesUp[i1][i2][i3]);		// up
               D[i1][i2][i3+1]	+= dFlow;						// maintain auxilary matrix
               d3s[i3+1]		+= dFlow;
            }

            pBV->flow	= fabs(dFlow);
            pBV->iDir	= dFlow>0;		// 1:outward, 0:inward
            pBV->pParent->pChild= pBV;
         }

   return true;
}


/************************************************************************
	initialize BVTree from the initial BF solution
	1. Determine the root
	2. Recursively construct the whole tree
 ************************************************************************/
void EmdL1::InitBVTree() {
   //- Using the center of the graph as the root
   int r	= (int)(0.5*m_n1-.5);
   int c	= (int)(0.5*m_n2-.5);
   int z	= (int)(0.5*m_n3-.5);
   m_pRoot	= m_nDim==2 ? &(m_Nodes[r][c]) : &(m_3dNodes[r][c][z]);
   m_pRoot->u		= 0;
   m_pRoot->iLevel	= 0;
   m_pRoot->pParent= NULL;
   m_pRoot->pPEdge	= NULL;

   //- Prepare a queue
   m_auxQueue[0]	= m_pRoot;
   int		nQueue	= 1;		// length of queue
   int		iQHead	= 0;		// head of queue

   //- Recursively build subtrees
   PEmdEdge	pCurE=NULL, pNxtE=NULL;
   PEmdNode	pCurN=NULL, pNxtN=NULL;
   int	iCurL	= -1;
   int	nBin	= m_n1*m_n2*MAX(m_n3,1);
   while(iQHead<nQueue && nQueue<nBin) {
      pCurN	= m_auxQueue[iQHead++];	// pop out from queue
      r		= pCurN->pos[0];
      c		= pCurN->pos[1];
      z		= pCurN->pos[2];

      // check connection from itself
      pCurE	= pCurN->pChild;	// the initial child from initial solution
      if(pCurE)	{
         pNxtN			= pCurE->pChild;
         pNxtN->pParent	= pCurN;
         pNxtN->pPEdge	= pCurE;
         m_auxQueue[nQueue++]	= pNxtN;
      }

      // check four neighbor nodes
      int	nNB	= m_nDim==2 ? 4 : 6;
      for(int k=0; k<nNB; ++k) {
         if(m_nDim==2)
            if(k==0 && c>0)				pNxtN	= &(m_Nodes[r][c-1]);		// left
            else if(k==1 && r>0)		pNxtN	= &(m_Nodes[r-1][c]);		// down
            else if(k==2 && c<m_n2-1)	pNxtN	= &(m_Nodes[r][c+1]);		// right
            else if(k==3 && r<m_n1-1)	pNxtN	= &(m_Nodes[r+1][c]);		// up
            else						continue;
         else if(m_nDim==3)
            if(k==0 && c>0)				pNxtN	= &(m_3dNodes[r][c-1][z]);		// left
            else if(k==1 && c<m_n2-1)	pNxtN	= &(m_3dNodes[r][c+1][z]);		// right
            else if(k==2 && r>0)		pNxtN	= &(m_3dNodes[r-1][c][z]);		// down
            else if(k==3 && r<m_n1-1)	pNxtN	= &(m_3dNodes[r+1][c][z]);		// up
            else if(k==4 && z>0)		pNxtN	= &(m_3dNodes[r][c][z-1]);		// shallow
            else if(k==5 && z<m_n3-1)	pNxtN	= &(m_3dNodes[r][c][z+1]);		// deep
            else						continue;

         if(pNxtN != pCurN->pParent) {
            pNxtE	= pNxtN->pChild;
            if(pNxtE && pNxtE->pChild==pCurN) {	// has connection
               pNxtN->pParent	= pCurN;
               pNxtN->pPEdge	= pNxtE;
               pNxtN->pChild	= NULL;
               m_auxQueue[nQueue++]	= pNxtN;

               pNxtE->pParent	= pCurN;			// reverse direction
               pNxtE->pChild	= pNxtN;
               pNxtE->iDir		= !pNxtE->iDir;

               if(pCurE)	pCurE->pNxt		= pNxtE;	// add to edge list
               else		pCurN->pChild	= pNxtE;
               pCurE	= pNxtE;
            }
         }
      }
   }
}



/************************************************************************
	Updating subtree recursively, including U values and level numbers
 ************************************************************************/
void EmdL1::UpdateSubtree(PEmdNode pRoot) {
   // Initialize auxiliary queue
   m_auxQueue[0]	= pRoot;
   int		nQueue	= 1;		// queue length
   int		iQHead	= 0;		// head of queue

   // BFS browing
   PEmdNode	pCurN=NULL,pNxtN=NULL;
   PEmdEdge	pCurE=NULL;
   while(iQHead<nQueue) {
      pCurN	= m_auxQueue[iQHead++];	// pop out from queue
      pCurE	= pCurN->pChild;

      // browsing all children
      while(pCurE) {
         pNxtN			= pCurE->pChild;
         pNxtN->iLevel	= pCurN->iLevel+1;
         pNxtN->u		= pCurE->iDir ? (pCurN->u - 1) : (pCurN->u + 1);
         pCurE			= pCurE->pNxt;
         m_auxQueue[nQueue++]	= pNxtN;
      }
   }
}


/************************************************************************
	Optimality Test
 ************************************************************************/
bool EmdL1::IsOptimal() {
   int		iC, iMinC = 0;
   PEmdEdge	pE;
   m_pEnter	= NULL;
   m_iEnter	= -1;

   // test each NON-BV edges
   for(int k=0; k<m_nNBV; ++k) {
      pE	= m_NBVEdges[k];
      iC	= 1 - pE->pParent->u + pE->pChild->u;
      if(iC<iMinC) 	{
         iMinC	= iC;
         m_iEnter= k;
      } else 	{
         // Try reversing the direction
         iC	= 1 + pE->pParent->u - pE->pChild->u;
         if(iC<iMinC) {
            iMinC	= iC;
            m_iEnter= k;
         }
      }
   }

   if(m_iEnter>=0) {
      m_pEnter = m_NBVEdges[m_iEnter];
      if(iMinC == (1 - m_pEnter->pChild->u + m_pEnter->pParent->u))	{
         // reverse direction
         PEmdNode	pN			= m_pEnter->pParent;
         m_pEnter->pParent	= m_pEnter->pChild;
         m_pEnter->pChild	= pN;
      }

      m_pEnter->iDir		= 1;
   }

   return m_iEnter==-1;
}


/************************************************************************
	Find new solution
 ************************************************************************/
void EmdL1::FindNewSolution() {
   // Find loop formed by adding the Enter BV edge.
   FindLoopFromEnterBV();

   // Modify flow values along the loop
   PEmdEdge	pE = NULL;
   EMDTYPE	minFlow		= m_pLeave->flow;
   int k;
   for(k=0; k<m_iFrom; ++k)  {
      pE	= m_fromLoop[k];
      if(pE->iDir)	pE->flow	+= minFlow;		// outward
      else			pE->flow	-= minFlow;		// inward
   }
   for(k=0; k<m_iTo; ++k)  {
      pE	= m_toLoop[k];
      if(pE->iDir)	pE->flow	-= minFlow;		// outward
      else			pE->flow	+= minFlow;		// inward
   }

   // Update BV Tree, removing the Leaving-BV edge
   PEmdNode	pLParentN	= m_pLeave->pParent;
   PEmdNode	pLChildN	= m_pLeave->pChild;
   PEmdEdge	pPreE		= pLParentN->pChild;
   if(pPreE==m_pLeave)
      pLParentN->pChild	= m_pLeave->pNxt;		// Leaving-BV is the first child
   else {
      while(pPreE->pNxt != m_pLeave)
         pPreE	= pPreE->pNxt;
      pPreE->pNxt	= m_pLeave->pNxt;				// remove Leaving-BV from child list
   }
   pLChildN->pParent	= NULL;
   pLChildN->pPEdge	= NULL;

   m_NBVEdges[m_iEnter]= m_pLeave;					// put the leaving-BV into the NBV array

   // Add the Enter BV edge
   PEmdNode	pEParentN	= m_pEnter->pParent;
   PEmdNode	pEChildN	= m_pEnter->pChild;

   m_pEnter->flow		= minFlow;
   m_pEnter->pNxt		= pEParentN->pChild;		// insert the Enter BV as the first child
   pEParentN->pChild	= m_pEnter;					//		of its parent

   // Recursively update the tree start from pEChildN
   PEmdNode	pPreN	= pEParentN;
   PEmdNode	pCurN	= pEChildN;
   PEmdNode	pNxtN;
   PEmdEdge	pNxtE, pPreE0;
   pPreE	= m_pEnter;
   while(pCurN) {
      pNxtN		= pCurN->pParent;
      pNxtE		= pCurN->pPEdge;
      pCurN->pParent	= pPreN;
      pCurN->pPEdge	= pPreE;

      if(pNxtN)	{
         // remove the edge from pNxtN's child list
         if(pNxtN->pChild==pNxtE)
            pNxtN->pChild	= pNxtE->pNxt;			// first child
         else {
            pPreE0	= pNxtN->pChild;
            while(pPreE0->pNxt != pNxtE)
               pPreE0	= pPreE0->pNxt;
            pPreE0->pNxt	= pNxtE->pNxt;			// remove Leaving-BV from child list
         }

         // reverse the parent-child direction
         pNxtE->pParent	= pCurN;
         pNxtE->pChild	= pNxtN;
         pNxtE->iDir		= !pNxtE->iDir;
         pNxtE->pNxt		= pCurN->pChild;
         pCurN->pChild	= pNxtE;

         pPreE		= pNxtE;
         pPreN		= pCurN;
      }

      pCurN		= pNxtN;
   }

   // Update U at the child of the Enter BV
   pEChildN->u			= m_pEnter->iDir ? (pEParentN->u - 1) : (pEParentN->u + 1);
   pEChildN->iLevel	= pEParentN->iLevel+1;
}


/************************************************************************
 Find the loop from the enteringBV edge and the Leaving-BV edge simultaneously.
 ************************************************************************/
void EmdL1::FindLoopFromEnterBV() {
   // Initialize Leaving-BV edge
   EMDTYPE		minFlow	= EMDINF;
   PEmdEdge	pE		= NULL;
   int			iLFlag	= 0;	// 0: in the FROM list, 1: in the TO list

   // Using two loop list to store the loop nodes
   PEmdNode	pFrom	= m_pEnter->pParent;
   PEmdNode	pTo		= m_pEnter->pChild;
   m_iFrom	= 0;
   m_iTo	= 0;
   m_pLeave= NULL;

   // Trace back to make pFrom and pTo at the same level
   while(pFrom->iLevel > pTo->iLevel)	{
      pE	= pFrom->pPEdge;
      m_fromLoop[m_iFrom++]	= pE;
      if(!pE->iDir && pE->flow<minFlow) {
         minFlow		= pE->flow;
         m_pLeave	= pE;
         iLFlag		= 0;	// 0: in the FROM list
      }
      pFrom	= pFrom->pParent;
   }

   while(pTo->iLevel > pFrom->iLevel)	{
      pE	= pTo->pPEdge;
      m_toLoop[m_iTo++]	= pE;
      if(pE->iDir && pE->flow<minFlow) {
         minFlow		= pE->flow;
         m_pLeave	= pE;
         iLFlag		= 1;	// 1: in the TO list
      }
      pTo	= pTo->pParent;
   }

   // Trace pTo and pFrom simultaneously till find their common ancester
   while(pTo!=pFrom) {
      pE	= pFrom->pPEdge;
      m_fromLoop[m_iFrom++]	= pE;
      if(!pE->iDir && pE->flow<minFlow)	{
         minFlow		= pE->flow;
         m_pLeave	= pE;
         iLFlag		= 0;	// 0: in the FROM list, 1: in the TO list
      }
      pFrom	= pFrom->pParent;

      pE	= pTo->pPEdge;
      m_toLoop[m_iTo++]	= pE;
      if(pE->iDir && pE->flow<minFlow) {
         minFlow		= pE->flow;
         m_pLeave	= pE;
         iLFlag		= 1;	// 0: in the FROM list, 1: in the TO list
      }
      pTo	= pTo->pParent;
   }

   // Reverse the direction of the Enter BV edge if necessary
   if(iLFlag==0)	{	// 0: in the FROM list
      PEmdNode	pN		= m_pEnter->pParent;
      m_pEnter->pParent	= m_pEnter->pChild;
      m_pEnter->pChild	= pN;
      m_pEnter->iDir		= !m_pEnter->iDir;
   }
}


/************************************************************************
	Computing the total flow as the final distance
 ************************************************************************/
EMDTYPE EmdL1::CompuTotalFlow() {
   EMDTYPE	f = 0;

   // Initialize auxiliary queue
   m_auxQueue[0]	= m_pRoot;
   int		nQueue	= 1;		// length of queue
   int		iQHead	= 0;		// head of queue

   // BFS browing the tree
   PEmdNode	pCurN=NULL,pNxtN=NULL;
   PEmdEdge	pCurE=NULL;
   while(iQHead<nQueue) {
      pCurN	= m_auxQueue[iQHead++];	// pop out from queue
      pCurE	= pCurN->pChild;

      // browsing all children
      while(pCurE) {
         f	+= pCurE->flow;
         pNxtN			= pCurE->pChild;
         pCurE			= pCurE->pNxt;
         m_auxQueue[nQueue++]	= pNxtN;
      }
   }

   return f;
}

