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

#ifndef binaryTree_H
#define binaryTree_H

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
    
    public:
        typedef binaryNode<CompType, ThermoType> bn;
        typedef chemPointISAT<CompType, ThermoType> chP;
        
    private:
	//Private Data
	
	//- Reference to the chemistryModel
	TDACChemistryModel<CompType, ThermoType>& chemistry_;	
        
	//- Root node of the binary tree
        bn *root_;
        
        //- Maximum number of elements in the binary tree
        label maxElements_;
        
        //- Size of the BST (= number of chemPoint stored)
        label size_;
        
        //- Secondary retrieve search variables
        label n2ndSearch_;
        label max2ndSearch_;
        
        label minBalanceThreshold_;
        label maxNbBalanceTest_;
        scalar balanceProp_;
        
        
        //- Insert the node newNode on the position specified of the parent binaryNode
        void insertNode
	(
         chP*& phi0,
         bn*& newNode
        );
        
        bool inSubTree
        (
            const scalarField& phiq, 
            bn* y,
            chP*& x
        );
        
        void deleteSubTree(binaryNode<CompType, ThermoType>* subTreeRoot);
                
        inline void deleteSubTree()
        {
            deleteSubTree(root_);
        }
        void transplant(bn* u, bn* v);
        
        chP* chemPSibling(bn* y);
        chP* chemPSibling(chP* x);
        
        bn* nodeSibling(bn* y);
        bn* nodeSibling(chP* x);
        
        void deleteAllNode(bn* subTreeRoot);

    public:
        
        
        //- Constructors
        
        //- Construct from dictionary and chemistryOnLineLibrary
        binaryTree (TDACChemistryModel<CompType, ThermoType>& chemistry,dictionary coeffsDict);
        
        
        inline label size()
        {
            return size_;
        } 
        
        label depth(bn* subTreeRoot);
        
        inline label depth()
        {
            return depth(root_);
        }
        
        inline bn* root()
        {
            return root_;
        }
        
        inline label maxElements()
        {
            return maxElements_;
        }
        
        //Insert a new leaf starting from the parent node of phi0
        //phi0 can be NULL
        void insertNewLeaf
        (
         const scalarField& phiq,
         const scalarField& Rphiq, 
         const List<List<scalar> >& A, 
         const scalarField& scaleFactor, 
         const scalar& epsTol,
         const label nCols,
               chP*& phi0
        );
        
        
        
        //Search the binaryTree until the nearest leaf of a specified
        //leaf is found. 
        void binaryTreeSearch(const scalarField& phiq, bn* node, chemPointBase*& nearest);

        //Perform a secondary binary tree search starting from a failed chemPoint x
        //If another candidate is found return true and x points to the chemPoint
        //The nUsed features is handled at ISAT level
        bool secondaryBTSearch(const scalarField& phiq,chP*& x);


        //- Delete a leaf from the binary tree and reshape the binary tree for the
        //  following binary tree search
        //  Return the index in the nodeList of the removed node (-1 when no node)
        void deleteLeaf(chP*& phi0);
        
        //- Cheap balance function
        //  This function just roughly separate the space in two parts 
        //  with a hyperplane perpendicular to the variance vector.
        //  Then, it repopulate the tree with this hyperplane stored at the root
        //  and by inserting the chemPoint randomly which should give a depth O(log_2(n))
        bool balance();
        
        inline void deleteAllNode()
        {
            deleteAllNode(root_);
        }
        
        chP* treeMin(bn* subTreeRoot);
        
        inline chP* treeMin()
        {
            return treeMin(root_);
        }
        
        chP* treeSuccessor(chP* x);
        
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
