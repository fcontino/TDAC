/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

    Data storage of the chemistryOnLineLibrary according to a binary 
    tree structure.
    
          0 (root node)
        /   \
      0      0
    /   \  /   \
   L    R L     0
              /   \
             L     R
             
L: elementLeft_
R: elementRight_

	The binaryTree class manages two DynamicList: nodeList_ and chemPointList_
	which store pointer to the nodes and leafs (chemPoint) in the binary tree.
	Therefore, there is only one instance of the node and leaf and the class
	using them (e.g. ISAT, chemistryModel,...) are only
	dealing with pointers
              
\*---------------------------------------------------------------------------*/

#ifndef binary_Tree_H
#define binary_Tree_H

#include "binaryNode.H"      
#include "chemPointISAT.H"
#include "scalarField.H"
#include "List.H"


namespace Foam
{

template<class CompType, class ThermoType>
class TDACChemistryModel;

template<class CompType, class ThermoType>
class binaryTree
{

private:
	//Private Data
	
	//- Reference to the chemistryModel
	TDACChemistryModel<CompType, ThermoType>& chemistry_;	
        
	//- Root node of the binary tree
        binaryNode<CompType, ThermoType> *root_;

        //- Maximum number of elements in the binary tree
        label maxElements_;
        
        //- Depth of the binary tree
        label depth_;

        //- Nodes of the tree are also stored in a list
        DynamicList<binaryNode<CompType, ThermoType>*> nodeList_;

        //- Chempoints of the tree are also stored in a list, used to find the less used element, to replace it when
        //  the library is full 
        DynamicList<chemPointISAT<CompType, ThermoType>*> chemPointISATList_;

    // Private Member Functions

        //- Clean the binary tree
        void clean(chemPointISAT<CompType, ThermoType>& oldPoint, binaryNode<CompType, ThermoType> *oldNode, const chemPointISAT<CompType, ThermoType>& newPoint, binaryNode<CompType, ThermoType> *parentNode);

        //- insert the node newNode on the position specified of the parent binaryNode
        void insertNode
	(
	    binaryNode<CompType, ThermoType>*& parent,
	    binaryNode<CompType, ThermoType>*& newNode,
	    label position
	);
		
	//- remove the link between the chemPoint phi0 and its parent node
	label detachLeaf(chemPointISAT<CompType, ThermoType>*& phi0);
		
		
		

   
public:

       //- Constructors
        
        //- Construct from dictionary and chemistryOnLineLibrary
        binaryTree( TDACChemistryModel<CompType, ThermoType>& chemistry,label maxElements);
        //binaryTree(label maxElements);
       //- Access
        
        inline DynamicList<chemPointISAT<CompType, ThermoType>*>& chemPointISATList()
        {
	        return chemPointISATList_;	
        }
        
        inline label size()
        {
        	return chemPointISATList_.size();
        } 
        
        inline label depth()
        {
        	return depth_;
        }  
        
        inline binaryNode<CompType, ThermoType>* root()
        {
        	return root_;
        }
                
       //- Edit

        //- Insert leaf by starting at the root node
		void insertNewLeaf
		(
			const scalarField& phiq,
			const scalarField& Rphiq, 
			const List<List<scalar> >& A, 
			const scalarField& scaleFactor, 
			const scalar& epsTol, 
			const label nCols
		);
		
		//Insert a new leaf starting from the parent node of phi0
		void insertNewLeaf
		(
			chemPointISAT<CompType, ThermoType>*& phi0, 
			const scalarField& phiq,
			const scalarField& Rphiq, 
			const List<List<scalar> >& A, 
			const scalarField& scaleFactor, 
			const scalar& epsTol,
			const label nCols
		);
		
		//Insert a new leaf starting at the root node with a specified index
		//in the chemPointList_ and the nodeList_
		void insertNewLeaf
		(
			const scalarField& phiq,
			const scalarField& Rphiq, 
			const List<List<scalar> >& A, 
			const scalarField& scaleFactor, 
			const scalar& epsTol,
			const label nCols,
			const label nodeIndex,
			const label chemPointISATIndex
		);
		
		//substract 1 to the depth of the binaryNode given in argument
		//and all the children node
		//return true if the branch had a depth equal to depth_
		bool decreaseDepth(binaryNode<CompType, ThermoType>*& bn);
		
		//loop over all the node of the nodeList (except nodeToAvoid)
		//and set depth_ to the deeper leaf
		void correctDepth(label nodeToAvoid);
        
       //-Search

		//Search the binaryTree until the nearest leaf of a specified
		//leaf is found. 
		void binaryTreeSearch(const scalarField& phiq, binaryNode<CompType, ThermoType>* node, chemPointBase*& nearest);

       //- Check
        
        //Report informations on the tree structure : 
        //	1) The number of leaves
        //	2) The maximum depth of the tree (root being level 0)
	//	3) The ratio nbLeaf/depth which is equal to the depth for a perfectly
	//		balanced tree and equal to (1+1/depth) when completly non balanced
	//Print those informations on Info
		void report();

	   //- Delete a leaf from the binary tree and reshape the binary tree for the
	   //  following binary tree search
	   //  Return the index in the nodeList of the removed node (-1 when no node)
	    label deleteLeaf(label chemPointISATIndex);

	//- deleteLeafCompact
	// CAUTION : MRU list is not updated after delete
	void deleteLeafCompact(label chemPointISATIndex);
	

       //- CleanAll
        void clear();

       //- ListFull
        bool isFull();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "binaryTree.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
