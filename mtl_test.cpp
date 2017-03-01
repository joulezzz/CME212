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





class IdentityMatrix {
  public:
  	IdentityMatrix(std::size_t n) : n_(n) {}

	/** Compute the product of a vector with this identity matrix
	 */

	/**	
	template <typename Vector>
	Vector operator*(const Vector& x) const {
		return x;
	}
	*/

	/** Get the dimension of the matrix */
	std::size_t get_dim() const{
		return n_;
	}

	/** Helper function to perfom multiplication. Allows for delayed 
	 *  evaluation of results.
	 *  Assign :: apply(a, b) resolves to an assignment operation such as 
	 *     a += b, a -= b, or a = b.
	 *  @pre @a size(v) == size(w) */
	template <typename VectorIn, typename VectorOut, typename Assign>
	void mult (const VectorIn& v, VectorOut& w, Assign) const {
	    assert(size(v) == size(w));
	    for (unsigned int i = 0; i < size(w); i++){
		Assign::apply(w[i],v[i]);
	    }
	}

	/** Matvec forward to MTL's lazy mat_cvec_multiplier oeprator */
	template <typename Vector> 
	mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
	operator*(const Vector& v) const {
		return {*this, v};
	}

  private:
  	// Empty
  	std::size_t n_;

};

inline std::size_t size(const IdentityMatrix& A){
	return A.get_dim()*A.get_dim();
}

inline std::size_t num_rows(const IdentityMatrix& A){
	return A.get_dim();
}

inline std::size_t num_cols(const IdentityMatrix& A){
	return A.get_dim();
}


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


using namespace mtl;
using namespace itl;

int main()
{
    // HW3: YOUR CODE HERE
    // Construct an IdentityMatrix and "solve" Ix = b
    // using MTL's conjugate gradient solver
    std::size_t size = 40, N = size*size;
    IdentityMatrix I(N);

    // Create an ILU(0) preconditioner
    pc::identity<IdentityMatrix>        P(I);

    // Set b such that x == 1 is solution; start with x == 0
    dense_vector<double>          x(N, 1.0), b(N);
    b= I * x; x= 0;

    // Termination criterion: r < 1e-6 * b or N iterations
    //noisy_iteration<double>       iter(b, 500, 1.e-6);
    cyclic_iteration<double> iter(b, 100, 1.e-10, 0, 100);

    // Solve Ax == b with left preconditioner P
    cg(I, x, b, P, iter);
std::cout << x << std::endl;





  return 0;
}
