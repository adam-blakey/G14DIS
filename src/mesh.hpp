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

	public:
		Elements* elements;
		Mesh(const int &a_noElements);
		~Mesh();

		int get_dimProblem();
		int get_noNodes();
		// Also store 'faces', with many 'face's.
};

#endif