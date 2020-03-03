/******************************************************************************
 * @details This is a file containing declarations of [mesh].
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
#ifndef CLASS_MESH
#define CLASS_MESH

#include "element.hpp"

class Mesh
{
	private:
		int noElements;
		int noNodes;
		int dimProblem;
		bool ownsElements;

	public:
		Elements* elements;
		Mesh(const int &a_noElements);
		Mesh(Elements* const &a_elements);
		~Mesh();

		int get_dimProblem() const;
		int get_noElements() const;
		int get_noNodes() const;
		// Also store 'faces', with many 'face's.
};

#endif