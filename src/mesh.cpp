#include "element.hpp"
#include "mesh.hpp"

/******************************************************************************
 * __Mesh__
 * 
 * @details 	The Mesh constructor, taking 1 argument for the number of elements.
 * 
 * @param[in] a_noElements 		The number of elements in this mesh. 	
 ******************************************************************************/
Mesh::Mesh(const int &a_noElements)
{
	this->noElements = a_noElements;
	this->noNodes    = a_noElements+1;
	this->dimProblem = a_noElements+1;
	this->elements = new Elements(a_noElements);
}

/******************************************************************************
 * __~Mesh__
 * 
 * @details 	Destroys memory for an instance of Mesh.
 ******************************************************************************/
Mesh::~Mesh()
{
	delete this->elements;
}