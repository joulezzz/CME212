/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <cassert>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL


/** Traits that MTL used to detmerine propperties of our IdentityMatrix. */
namespace mtl {
namespace ashape {

/** Define IdentityMatri to be a non-scalar type. */
template<>
struct ashape_aux<IdentityMatrix>{
	typedef nonscal type;
};
}

/** IdentityMatrix implements the Collection concept 
 * with value_type and size_type */
template<>
struct Collection<IdentityMatrix> {
	typedef double value_type;
	typedef unsigned size_type;
};
} // end namespace



class IdentityMatrix {
  public:
  	IdentityMatrix(std::size_t nrows, std::size_t ncols) : nrows_(nrows), ncols_(ncols) {}

	/** Compute the product of a vector with this identity matrix
	 */
	template <typename Vector>
	Vector operator*(const Vector& x) const {
		return x;
	}

	/** The numner of elements in the matrix. */
	inline std::size_t size(const IdentityMatrix& A){
		return num_rows_*num_cols_;
	}
	/** The numbner of rows in the matrix. */
	inline std::size_t num_rows(const IdentityMatrix& A){
		return num_rows_;
	}
	/** The number of columns in the matrix. */
	inline std::size_t num_cols(const IdentityMatrix& A){
		return num_cols_;
	}

	/** Helper function to perfom multiplication. Allows for delayed 
	 *  evaluation of results.
	 *  Assign :: apply(a, b) resolves to an assignment operation such as 
	 *     a += b, a -= b, or a = b.
	 *  @pre @a size(v) == size(w) */
	template <typename VectorIn, typename VectorOut, typename Assign>
	void mult (const VectorIn& v, VectorOut& w, Assign) const {
	    assert(size(v) == size(w));
	    Assign::apply(w,v);
	}

	/** Matvec forward to MTL's lazy mat_cvec_multiplier oeprator */
	template <typename Vector> 
	mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
	oeprator*(const Vector& v) const {
		return {*this, v};
	}


  private:
  	// Empty
  	std::size_t nrows_;
  	std::size_t ncols_;

}



int main()
{
  // HW3: YOUR CODE HERE
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver

  return 0;
}
