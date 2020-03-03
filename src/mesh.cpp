/******************************************************************************
 * @details This is a file containing declarations of [mesh].
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
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
	this->ownsElements = true;
}

/******************************************************************************
 * __Mesh__
 * 
 * @details 	The Mesh constructor, taking 1 argument for the elements.
 * 
 * @param[in] a_elements 		The elements to create with the mesh.
 ******************************************************************************/
Mesh::Mesh(Elements* const &a_elements)
{
	this->noElements = a_elements->get_noElements();
	this->noNodes    = a_elements->get_noElements()+1;
	this->dimProblem = a_elements->get_noElements()+1;
	this->elements = a_elements;
	this->ownsElements = false;
}

/******************************************************************************
 * __~Mesh__
 * 
 * @details 	Destroys memory for an instance of Mesh.
 ******************************************************************************/
Mesh::~Mesh()
{
	if (this->ownsElements)
	{
		delete this->elements;
	}
}

/******************************************************************************
 * __get_dimProblem__
 * 
 * @details 	Returns the value of the private variable 'dimProblem'.
 ******************************************************************************/
int Mesh::get_dimProblem() const
{
	return this->dimProblem;
}

/******************************************************************************
 * __get_noElements__
 * 
 * @details 	Returns the value of the private variable 'noElements'.
 ******************************************************************************/
int Mesh::get_noElements() const
{
	return this->noElements;
}

/******************************************************************************
 * __get_noNodes__
 * 
 * @details 	Returns the value of the private variable 'noNodes'.
 ******************************************************************************/
int Mesh::get_noNodes() const
{
	return this->noNodes;
}