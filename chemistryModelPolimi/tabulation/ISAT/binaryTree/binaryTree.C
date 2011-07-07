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

\*---------------------------------------------------------------------------*/

#include "binaryTree.H"
#include "binaryNode.H"
#include "demandDrivenData.H"
#include "clockTime.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
binaryTree<CompType, ThermoType>::binaryTree
(
 const TDACChemistryModel<CompType, ThermoType>& chemistry,
 label maxElements
	
) 
:
   chemistry_(chemistry),
   root_(NULL),
   maxElements_(maxElements),
   depth_(0)
{


}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * Insert New Leaf * * * * * * //
/*---------------------------------------------------------------------------*\
	Insert a new leaf starting from the parent node of phi0
	Input : phi0 the leaf to replace by a node
			phiq the new composition to store
			Rphiq the mapping of the new composition point
			A the mapping gradient matrix
			B the matrix used to initialize the EOA
			nCols the size of the matrix
	Output: void
	Description : 
			1) Create a new leaf with the data to initialize the EOA and to
				retrieve the mapping by linear interpolation (the EOA is 
				initialize in the chemPoint constructor)
			2) Get the parent node of phi0 and connect a new node in place of the 
				leaf of phi0. This new node is constructed with phi0 on the left
				and phiq on the right (the hyperplane is computed inside the 
				binaryNode constructor)
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
void binaryTree<CompType, ThermoType>::insertNewLeaf
(
	chemPointISAT<CompType, ThermoType>*& phi0, 
	const scalarField& phiq,
	const scalarField& Rphiq, 
	const List<List<scalar> >& A, 
	const scalarField& scaleFactor, 
	const scalar& epsTol,
	const label nCols
)
{
	//create the new chemPoint which holds the composition point
	//phiq and the data to initialize the EOA
	label n = chemPointISATList_.size();
    	//chemPointList_.setSize(n+1);
    	chemPointISATList_.append(new chemPointISAT<CompType, ThermoType>(chemistry_,phiq, Rphiq, A, scaleFactor, epsTol, nCols));  //DynamicList check the size before appending

	//remove leaf from parent node and return its previous position
	//on this node (left=0 or right=1)
	label position = detachLeaf(phi0);

	//insert new node on the parent node in the position of the
	//previously stored leaf (phi0)
	//the new node contains phi0 on the left and phiq on the right
	//the hyper plane is computed in the binaryNode constructor
	label m = nodeList_.size();
	//nodeList_.setSize(m+1);
	//the new binaryNode is initialized with two NULL pointer for the left and right
	//node and with phi0 on the left and phiq on the right
	//the depth of the newly created node is its parent depth + 1
	nodeList_.append(new binaryNode<CompType, ThermoType>(phi0, chemPointISATList_[n], phi0->node(), phi0->node()->depth()+1)); //DynamicList check the size before appending
	phi0->position() = binaryNode<CompType, ThermoType>::LEFT; //phi0 on the left
	chemPointISATList_[n]->position() = binaryNode<CompType, ThermoType>::RIGHT; //new leaf on the right
	chemPointISATList_[n]->listIndex() = n;
	insertNode(phi0->node(), nodeList_[m], position);
	phi0->node()=nodeList_[m];
	nodeList_[m]->listIndex()=m;
	chemPointISATList_[n]->node()=nodeList_[m];
}
        

//Insert a new leaf starting from the root of the binary tree
//(used after cleaning up the binary tree)
template<class CompType, class ThermoType>
void binaryTree<CompType, ThermoType>::insertNewLeaf
(
	const scalarField& phiq,
	const scalarField& Rphiq, 
	const List<List<scalar> >& A, 
	const scalarField& scaleFactor, 
	const scalar& epsTol, 
	const label nCols
)
{
clockTime cpuTimeBST=clockTime();	
cpuTimeBST.timeIncrement();//reset
	if (chemPointISATList_.size()==0) //no chemPoint stored : the tree is empty
	{
        
		//chemPointISATList_.setSize(1);
		chemPointISATList_.append(new chemPointISAT<CompType, ThermoType>(chemistry_,phiq, Rphiq, A, scaleFactor, epsTol, nCols));
		chemPointISATList_[0]->listIndex() = 0;
		//note : depth of the tree remains 0
		//		 position of the leaf is left by default
	}
	else if (chemPointISATList_.size()!=0 && nodeList_.size()==0) //no node stored : the tree contains only one leaf
	{
		//create the second chemPoint of the tree
		//chemPointISATList_.setSize(2);
		chemPointISATList_.append(new chemPointISAT<CompType, ThermoType>(chemistry_,phiq, Rphiq, A, scaleFactor, epsTol, nCols));
		
		//create the first node with the two first chemPoint
		//the last chemPoint added, on the right
		//nodeList_.setSize(1);
		//depth of the root node is 0
		nodeList_.append(new binaryNode<CompType, ThermoType>(chemPointISATList_[0], chemPointISATList_[1], 0));
		chemPointISATList_[0]->position() = binaryNode<CompType, ThermoType>::LEFT;
		chemPointISATList_[1]->position() = binaryNode<CompType, ThermoType>::RIGHT;
		chemPointISATList_[1]->listIndex() = 1;
		//depth of the tree is 1
		depth_ = 1;
		//associate root_ pointer with this node
		root_ = nodeList_[0];
		root_->listIndex()=0;
		chemPointISATList_[0]->node() = root_;
		chemPointISATList_[1]->node() = root_;		
	}
	//Travel until we reach the point where (according to a binary tree search)
	//the leaf is replaced by a node separating the composition space between
	//the new and previously stored leaf
	else
	{
		chemPointBase* phi0Base;
		binaryTreeSearch(phiq, root_,phi0Base);
		chemPointISAT<CompType, ThermoType>* phi0 = dynamic_cast<chemPointISAT<CompType, ThermoType>*>(phi0Base);
		insertNewLeaf(phi0, phiq, Rphiq, A, scaleFactor, epsTol, nCols);
	}
}

//Insert a new leaf starting at the root node with a specified index
//in the chemPointList_ and the nodeList_
template<class CompType, class ThermoType>
void binaryTree<CompType, ThermoType>::insertNewLeaf
(
	const scalarField& phiq,
	const scalarField& Rphiq, 
	const List<List<scalar> >& A, 
	const scalarField& scaleFactor, 
	const scalar& epsTol,
	const label nCols,
	const label nodeIndex,
	const label chemPointISATIndex
)
{
    //when node Index is -1, no node has been deleted 
    //and then no node is added
    if (nodeIndex==-1)
    {
        chemPointISATList_[chemPointISATIndex] = new chemPointISAT<CompType, ThermoType>(chemistry_,phiq, Rphiq, A, scaleFactor, epsTol, nCols);
        chemPointISATList_[chemPointISATIndex]->listIndex()=chemPointISATIndex;
    }
    //otherwise a binaryTreeSearch is performed to find the nearest 
    //leaf stored and a new chemPoint is added to the tree as usual 
    //but the chemPointList_ and nodeList_ are not expanded, the 
    //provided index are used instead
    else if (root_ == NULL) //only one leaf
    {
        //if root_ ==NULL after a deleteLeaf, it is either because chemPointList_[0] or [1] has been deleted
        //then, chemPointList_[1-chemPointIndex] is the other one
        chemPointISATList_[chemPointISATIndex] = new chemPointISAT<CompType, ThermoType>(chemistry_,phiq, Rphiq, A, scaleFactor, epsTol, nCols);
        nodeList_[nodeIndex] = new binaryNode<CompType, ThermoType>(chemPointISATList_[1-chemPointISATIndex], chemPointISATList_[chemPointISATIndex],  0);
        root_ = nodeList_[nodeIndex];
        root_->listIndex()=nodeIndex;
        chemPointISATList_[chemPointISATIndex]->listIndex()=chemPointISATIndex;
        chemPointISATList_[chemPointISATIndex]->node()=root_;
        chemPointISATList_[1-chemPointISATIndex]->node()=root_;
        chemPointISATList_[chemPointISATIndex]->position()=binaryNode<CompType, ThermoType>::RIGHT;
        chemPointISATList_[1-chemPointISATIndex]->position()=binaryNode<CompType, ThermoType>::LEFT;
        depth_ = 1;
    }
    else
    {
        chemPointBase* phi0Base;
        binaryTreeSearch(phiq, root_,phi0Base);
        chemPointISAT<CompType, ThermoType>* phi0 = dynamic_cast<chemPointISAT<CompType, ThermoType>*>(phi0Base);
        
        //a new node is created in the chemPointIndex
        chemPointISATList_[chemPointISATIndex] = new chemPointISAT<CompType, ThermoType>(chemistry_,phiq, Rphiq, A, scaleFactor, epsTol, nCols);
        
        //remove leaf from parent node and return its previous position
        //on this node (left=0 or right=1)
        label position = detachLeaf(phi0);
        
        //insert new node on the parent node in the position of the
        //previously stored leaf (phi0)
        //the new node contains phi0 on the left and phiq on the right
        //the hyper plane is computed in the binaryNode constructor
        //the new binaryNode is initialized with two NULL pointer for the left and right
        //node and with phi0 on the left and phiq on the right
        //the depth of the newly created node is its parent depth + 1
        nodeList_[nodeIndex] = new binaryNode<CompType, ThermoType>(phi0, chemPointISATList_[chemPointISATIndex], phi0->node(), phi0->node()->depth()+1);
        phi0->position() = binaryNode<CompType, ThermoType>::LEFT; //phi0 on the left
        chemPointISATList_[chemPointISATIndex]->position() = binaryNode<CompType, ThermoType>::RIGHT; //new leaf on the right
        chemPointISATList_[chemPointISATIndex]->listIndex() = chemPointISATIndex;
        insertNode(phi0->node(), nodeList_[nodeIndex], position);
        phi0->node()=nodeList_[nodeIndex];
        nodeList_[nodeIndex]->listIndex()=nodeIndex;
        chemPointISATList_[chemPointISATIndex]->node()=nodeList_[nodeIndex];
    }
}		

      
        
// * * * * * * end of insertNewLeaf functions * * * * * * //



//remove the link to the leaf phi0 in the parent node of phi0
//phi0.position() is =0 if phi0 is on the left and =1 if on the right
template<class CompType, class ThermoType>
label binaryTree<CompType, ThermoType>::detachLeaf(chemPointISAT<CompType, ThermoType>*& phi0)
{
//	binaryNode*& parent = phi0->node();
//	label position = phi0->position();
	if (phi0->position())//position !=0 => on the right
	{
		phi0->node()->elementRight() = NULL;
	}
	else //position == 0 => on the left
	{
		phi0->node()->elementLeft() = NULL;
	}
	return phi0->position();
}

//insert new node on the position specified of the parent node
template<class CompType, class ThermoType>
void binaryTree<CompType, ThermoType>::insertNode
(
	binaryNode<CompType, ThermoType>*& parent,
	binaryNode<CompType, ThermoType>*& newNode,
	label position
)
{
	if (position)//position !=0 => on the right
	{
		parent->right() = newNode;
		newNode->position()=binaryNode<CompType, ThermoType>::RIGHT;
	}
	else //position == 0 => on the left
	//position of a new node is left by default
	{
		parent->left() = newNode;
		
	}
	//depth_ is the maximum depth of the leaves (depth of the node +1)
	depth_ = max(depth_,newNode->depth()+1);
}

//Search the binaryTree until the nearest leaf of a specified
//leaf is found. 
template<class CompType, class ThermoType>
void binaryTree<CompType, ThermoType>::binaryTreeSearch(const scalarField& phiq, binaryNode<CompType, ThermoType>* node, chemPointBase*& nearest)
{
    if (nodeList_.size()!=0)
    {
        scalar vPhi=0.0;
        const scalarField v = node->v();
        const scalar a = node->a();
        for (label i=0; i<phiq.size(); i++) vPhi += phiq[i]*v[i];
        if(vPhi > a) //on rigth
        {
            //	Info << "vPhiRight = " << vPhi << endl;
            if (node->right()!=NULL)
            {
                binaryTreeSearch(phiq, node->right(), nearest);
            }
            else //the terminal node is reached, return element on right
            {
                nearest = node->elementRight_;
            }
        }
        else //on left
        {
            //	Info << "vPhiLeft = " << vPhi << endl;
            if (node->left()!=NULL)
            {
                binaryTreeSearch(phiq, node->left(), nearest);
            }
            else //the terminal node is reached, return element on right
            {
                nearest = node->elementLeft_;
            }
        }
    }
    // nodeList_ = 0, if chemPointList !=0, there is only one leaf
    else if (chemPointISATList_.size()!=0)
    {
        nearest = chemPointISATList_[0];
    }
    else
    {
        nearest = NULL;
    }
}//end binaryTreeSearch


//- Delete a leaf from the binary tree and reshape the binary tree for the
//  following binary tree search
//  Return the index in the nodeList of the removed node (-1 when no node)
//Note on the correction of the depth_ variable: 2 possibilities	
//	1) the other branch connected to the parent node of the removed leaf 
//	   has a depth that is lower than the depth_ variable
//	   => the depth of the tree is defined elsewhere, there is no need
//	   to change it
//  2) the other branch's depth is equal to depth_
//	   => the tree may not have another branch with the same depth
//	   it should then been check and correct depth_ if needed  
template<class CompType, class ThermoType>
label binaryTree<CompType, ThermoType>::deleteLeaf(label chemPointISATIndex)
{
label parentNodeIndex = -1;

if (nodeList_.size()==0) //no node stored : the tree contains only one leaf
{
    //if no node, only delete the chemPoint
    deleteDemandDrivenData(chemPointISATList_[chemPointISATIndex]) ;
    
}
else
{
    //find the index of the parent node
    parentNodeIndex = chemPointISATList_[chemPointISATIndex]->node()->listIndex();
    //find the other element (node or leaf) connected to the node 
    //(p: the position of the leaf to delete -> pos. of the other elem. = (1-p))
    label otherElementPosition = 1-chemPointISATList_[chemPointISATIndex]->position();
    //delete the chemPoint
    deleteDemandDrivenData(chemPointISATList_[chemPointISATIndex]) ;
    
    //get the depth of the parent node
    label parentNodeDepth = nodeList_[parentNodeIndex]->depth();
    
    //get the position of the parent node 
    label parentNodePosition = nodeList_[parentNodeIndex]->position();
    //delete the parent node replace it by the other element
    //to know the other element type just see if the parent node 
    //hold a pointer to a node in the position of the other element
    //if not, then it is a leaf
    //change the position of the other element
    
    bool doCorrectDepth = false;
    
    if (otherElementPosition)//!=0 => the other element is on the right
    {
        
        //the other element is on the right and is a node, it will then take the place
        //of the parent node and its depth and the depth of its children are decreased	
        if (nodeList_[parentNodeIndex]->right()!=NULL) //the other element is a node
        {	
            
            //the parent node is the root node, the other element which is a node
            //take the place of root
            if (parentNodeDepth==0)
            {
                root_= nodeList_[parentNodeIndex]->right();
                nodeList_[parentNodeIndex]->right()->parent() = NULL;
                //put root back to default left position
                root_->position()=binaryNode<CompType, ThermoType>::LEFT;
            }
            else if (parentNodePosition)//position !=0 => on the right
            {
                nodeList_[parentNodeIndex]->parent()->right()=nodeList_[parentNodeIndex]->right();
                nodeList_[parentNodeIndex]->right()->parent()=nodeList_[parentNodeIndex]->parent();
                //the moved node was already on the right of the deleted node
            }
            else //position == 0 => on the left
            {
                nodeList_[parentNodeIndex]->parent()->left()=nodeList_[parentNodeIndex]->right();
                nodeList_[parentNodeIndex]->right()->parent()=nodeList_[parentNodeIndex]->parent();
                //set position to left
                nodeList_[parentNodeIndex]->right()->position()=binaryNode<CompType, ThermoType>::LEFT;
            }
            //correct the depth of the moved node and all the children nodes
            doCorrectDepth = decreaseDepth(nodeList_[parentNodeIndex]->right());
        }
        
        //the other element is on the right and is a leaf, the parent node is removed and the link
        //with the great parent is created for the remaining leaf	
        else //the other element is a leaf
        {
            if (parentNodeDepth==0)
            {
                root_ = NULL;
                //set position to the default left position
                nodeList_[parentNodeIndex]->elementRight()->position() = binaryNode<CompType, ThermoType>::LEFT;
                //reset the pointer to the parent node
                nodeList_[parentNodeIndex]->elementRight()->node()=NULL;
            }
            else if (parentNodePosition)//position !=0 => on the right
            {
                nodeList_[parentNodeIndex]->parent()->elementRight()=nodeList_[parentNodeIndex]->elementRight();
                //the pointer to the node is set back to NULL
                nodeList_[parentNodeIndex]->parent()->right()=NULL;
                nodeList_[parentNodeIndex]->elementRight()->node()=nodeList_[parentNodeIndex]->parent();
                //the moved leaf was already on the right of the deleted node
            }
            else //position == 0 => on the left
            {
                nodeList_[parentNodeIndex]->parent()->elementLeft()=nodeList_[parentNodeIndex]->elementRight();
                //the pointer to the node is set back to NULL					
                nodeList_[parentNodeIndex]->parent()->left()=NULL;
                nodeList_[parentNodeIndex]->elementRight()->node()=nodeList_[parentNodeIndex]->parent();
                nodeList_[parentNodeIndex]->elementRight()->position()=binaryNode<CompType, ThermoType>::LEFT;
            }			
            doCorrectDepth = (nodeList_[parentNodeIndex]->depth()+1 == depth_);
        }
    }//end if on the otherElementPosition
    
    else //position == 0 => on the left
    {
        if (nodeList_[parentNodeIndex]->left()!=NULL) //the other element is a node
        {
            if (parentNodeDepth==0)
            {
                root_= nodeList_[parentNodeIndex]->left();
                nodeList_[parentNodeIndex]->left()->parent() = root_;
            }
            else if (parentNodePosition)//position !=0 => on the right
            {
                nodeList_[parentNodeIndex]->parent()->right()=nodeList_[parentNodeIndex]->left();
                nodeList_[parentNodeIndex]->left()->parent()=nodeList_[parentNodeIndex]->parent();
                nodeList_[parentNodeIndex]->left()->position()=binaryNode<CompType, ThermoType>::RIGHT;
            }
            else //position == 0 => on the left
            {
                nodeList_[parentNodeIndex]->parent()->left()=nodeList_[parentNodeIndex]->left();
                nodeList_[parentNodeIndex]->left()->parent()=nodeList_[parentNodeIndex]->parent();
                //the moved node was already on the left of the deleted node
            }
            //correct the depth of the moved node and all the children nodes
            doCorrectDepth = decreaseDepth(nodeList_[parentNodeIndex]->left());	
        }
        else //the other element is a leaf
        {
            if (parentNodeDepth==0)
            {
                root_ = NULL;
                //reset the pointer to the parent node
                nodeList_[parentNodeIndex]->elementLeft()->node()=NULL;
            }			
            else if (parentNodePosition)//position !=0 => on the right
            {
                nodeList_[parentNodeIndex]->parent()->elementRight()=nodeList_[parentNodeIndex]->elementLeft();
                //the pointer to the node is set back to NULL
                nodeList_[parentNodeIndex]->parent()->right()=NULL;					
                nodeList_[parentNodeIndex]->elementLeft()->node()=nodeList_[parentNodeIndex]->parent();
                nodeList_[parentNodeIndex]->elementLeft()->position()=binaryNode<CompType, ThermoType>::RIGHT;
            }
            else //position == 0 => on the left
            {
                nodeList_[parentNodeIndex]->parent()->elementLeft()=nodeList_[parentNodeIndex]->elementLeft();
                //the pointer to the node is set back to NULL					
                nodeList_[parentNodeIndex]->parent()->left()=NULL;					
                nodeList_[parentNodeIndex]->elementLeft()->node()=nodeList_[parentNodeIndex]->parent();
                //the moved leaf was already on the left of the deleted node
            }
            doCorrectDepth = (nodeList_[parentNodeIndex]->depth()+1 == depth_);
        }
    }
    
    //delete the parent node
    deleteDemandDrivenData(nodeList_[parentNodeIndex]) ;	
    if (doCorrectDepth)//the moved branch had a depth equal to depth_
    {
        correctDepth(parentNodeIndex);
    }
}//end else (at leat one node)
//return the deleted node index (-1 when no node)
return parentNodeIndex;
}//end of deleteLeaf

//- deleteLeafCompact
// CAUTION : MRU list is not updated after delete
template<class CompType, class ThermoType>
void binaryTree<CompType, ThermoType>::deleteLeafCompact(label chemPointISATIndex)
{
    label nodeIndex = deleteLeaf(chemPointISATIndex);
    //the tree should at least have two leaves
    //if the last is removed, nothing has to be done except descreasing posChemP_
    if (chemPointISATIndex != chemPointISATList_.size()-1)
    {
        chemPointISATList_[chemPointISATIndex] = chemPointISATList_.remove();
        chemPointISATList_[chemPointISATIndex]->listIndex() = chemPointISATIndex;
    }
    else //if only one leaf, it is deleted and posChemP is decreased by one
    {
        chemPointISATList_.remove();
    }

    if((nodeIndex != -1) && (nodeIndex != nodeList_.size()-1))
    {
        nodeList_[nodeIndex] = nodeList_.remove();
        nodeList_[nodeIndex]->listIndex() = nodeIndex;
    }
    else if((nodeIndex != -1) && (nodeIndex == nodeList_.size()-1))
    {
        nodeList_.remove();
    }
}



//substract 1 to the depth of the binaryNode given in argument
//and all the children node
//return true if the branch had a depth equal to depth_
template<class CompType, class ThermoType>
bool binaryTree<CompType, ThermoType>::decreaseDepth(binaryNode<CompType, ThermoType>*& bn)
{
	bool deeper=(bn->depth()+1 == depth_);
	bn->depth()--;
	if (bn->right()!=NULL) deeper = decreaseDepth(bn->right());
	if (bn->left()!=NULL) deeper = deeper || decreaseDepth(bn->left()); // OR operator if the right branch was deeper
	return deeper;
}



//loop over all the node of the nodeList (except nodeToAvoid)
//and set depth_ to the deeper leaf
template<class CompType, class ThermoType>
void binaryTree<CompType, ThermoType>::correctDepth(label nodeToAvoid)
{
	//the function loop in the nodeList to see if the depth of the tree
	//is changed when removing this leaf (depth of the leaf = depth of the 
	//parent node +1)
	depth_ = 0;
	for (label i=0; i<nodeList_.size(); i++)
	{
		if (i != nodeToAvoid && depth_ < nodeList_[i]->depth()+1)
		{
			depth_ = nodeList_[i]->depth()+1;
		}	
	}
}


//Remove everything entry of the tree and delete the associated objects
template<class CompType, class ThermoType>
void binaryTree<CompType, ThermoType>::clear()
{

	//reset root node
    root_=NULL;
    
    //reset depth
    depth_ = 0;
    
    //delete the node objects and remove from nodeList_
    forAll(nodeList_, nodeI)
    {
        deleteDemandDrivenData(nodeList_[nodeI]);
    }
    nodeList_.clear();

    //delete the chemPoint objects and remove from chemPointList_
    forAll(chemPointISATList_, pointI)
    {
        deleteDemandDrivenData(chemPointISATList_[pointI]) ;   
    }
    chemPointISATList_.clear();
}//end cleanAll

//Check if the tree has reached the maximum number of elements
template<class CompType, class ThermoType>
bool binaryTree<CompType, ThermoType>::isFull()
{
    label n = chemPointISATList_.size();
    
    if(n < maxElements_)
    {
        return false;
    }
    else
    {   
        return true;
    }
      
}//end isListFull

template<class CompType, class ThermoType>
void binaryTree<CompType, ThermoType>::report()
{
	Info << "Binary tree report :" << endl;
	Info << "Total number of leaf : " << size() << endl;
	Info << "Depth of the binary tree : " << depth() << endl;
	Info << "Ratio nbLeaf/depth : ";
	if (depth() == 0 && size() == 0)
	{
		Info << "empty tree" << endl;
	}
	else if (size() == 1 && depth() == 0)
	{
		Info << "one leaf tree" << endl;
	}
	else
	{
		Info << size()/depth() << endl;
	}
	Info << "End of binary tree report" << endl;
}
} // End namespace Foam


