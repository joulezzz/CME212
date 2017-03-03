/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SFML_Viewer to visualize the solution.
 */

#include <fstream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <cassert>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"
#include <math.h>



#include "Graph.hpp"

// HW3: YOUR CODE HERE
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
using GraphType = Graph<char,char>;  //<  DUMMY Placeholder
using NodeType  = typename GraphType::node_type;


 // Functor that colors the plot in an interesting way
  struct ColorFn {
    // operator
    CME212::Color operator () (NodeType n){ 
      return CME212::Color::make_heat((v_[n.index()] - mtl::min(v_))/(mtl::max(v_) - mtl::min(v_) + 0.00000001));
    }
    // Constructor
    ColorFn(mtl::dense_vector<double>& v) : v_(v) {};
   private :
    mtl::dense_vector<double>& v_;
  };

  // Functor that changes how nodes are plotted, so soln can be viewed on z-axis
  struct NodePosition {
    // Constructor 
    NodePosition(mtl::dense_vector<double>& x) : x_(x) {}
    // operator
    Point operator () (NodeType n){
return Point(n.position().x, n.position().y, x_[n.index()]);
    }
   private:
     mtl::dense_vector<double>& x_;
  };



/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @param[in,out] g  The Graph to remove nodes from
 * @param[in]    bb  The BoundingBox, all nodes inside this box will be removed
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  // HW3: YOUR CODE HERE
  (void) g; (void) bb;   //< Quiet compiler
  for (auto it = g.node_begin(); it != g.node_end();){
    auto n = *it;
    if (bb.contains(n.position())){
      g.remove_node(n);
    }
    else {
      ++it;
    }
  }
  //return;
}


/** g(x), boundary conditions */
double g_boundary(const NodeType n){
  // first check if on boundary
  CME212::BoundingBox<Point> thisbox = Box3D(Point(-0.6,-0.2,-1), Point( 0.6, 0.2,1));
  if (norm_inf(n.position()) == 1){
    return 0.0;
  }
  else if ( (norm_inf(n.position() - Point(0.6,0.6,0)) < 0.2) || (norm_inf(n.position() - Point(-0.6,0.6,0)) < 0.2) || (norm_inf(n.position() - Point(0.6,-0.6,0)) < 0.2) || (norm_inf(n.position() - Point(-0.6,-0.6,0)) < 0.2) ){
    return -0.2;
  }
  else if ( thisbox.contains(n.position()) ) {
    return 1.0;
  }
  // if not on boundary then use forcing function
  else {
    return -1.0;
  }
}

/** f(x), forcing function */
double f(const NodeType& n){
  return 5.0*cos( norm_1(n.position()) );
}

/** b, RHS of system Ax = b */
double b(const NodeType& i, const GraphType& graph){
	double gofx_i = g_boundary(i);
	if (gofx_i != double(-1)){
		return gofx_i;
	}
	else {
		double sum = 0;
		for (auto eit = i.edge_begin(); eit != i.edge_end(); ++eit){
			auto j = (*eit).node2();
			double gofx_j = g_boundary(j);
			if (gofx_j != double(-1)){
				sum += gofx_j;
			}
		}
		//std::cout << graph.edge(0).length() << std::endl;
		return pow( double(graph.edge(0).length()) ,2)*f(i) - sum;
	}
}



class GraphSymmetricMatrix{
  public:
    GraphSymmetricMatrix(GraphType& g) : g_(g) {}

    std::size_t get_dim() const{
    return g_.num_nodes();
    }

    // L(i,j), Discrete Matrix Approximating Laplace Operator
    double L(NodeType i, NodeType j) const{
      if (i == j){
        //std::cout << i.degree() << std::endl;
        return double(-1*int(i.degree()));
      }
      else if ( g_.has_edge(i,j) || g_.has_edge(j,i) ) {
        return double(1);
      }
      else {
        return double(0);
      }
    }

    // A(i,j), Linear System of Equations
    double A(NodeType i, NodeType j) const {
      if ( (i == j) && (g_boundary(i) != -1) ){
        return double(1);
      }
      else if ( (i != j) && ( (g_boundary(i) != -1) || (g_boundary(j) != -1) ) ){
        return double(0);
      }
      else {
        return L(i,j);
      }
    }

    /** Helper function to perfom multiplication. Allows for delayed 
     *  evaluation of results.
     *  Assign :: apply(a, b) resolves to an assignment operation such as 
     *     a += b, a -= b, or a = b.
     *  @pre @a size(v) == size(w) */
    template <typename VectorIn, typename VectorOut, typename Assign>
    void mult (const VectorIn& v, VectorOut& w, Assign) const {
        assert(size(v) == size(w));
        
        for (auto nit = g_.node_begin(); nit != g_.node_end(); ++nit){
          auto i = *nit;
    double temp = 0.0;
          for(auto njt = g_.node_begin(); njt != g_.node_end(); ++njt){
            auto j = *njt;
            //auto j = e.node2();
            temp += A(i,j)*v[j.index()];
          }
          Assign::apply(w[i.index()], temp);
        }
    }

    /** Matvec forward to MTL's lazy mat_cvec_multiplier oeprator */
    template <typename Vector> 
    mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
    operator*(const Vector& v) const {
      return {*this, v};
    }

  private:

    GraphType& g_;    
};


inline std::size_t size(const GraphSymmetricMatrix& M){
	return M.get_dim()*M.get_dim();
}

inline std::size_t num_rows(const GraphSymmetricMatrix& M){
	return M.get_dim();
}

inline std::size_t num_cols(const GraphSymmetricMatrix& M){
	return M.get_dim();
}


/** Traits that MTL used to detmerine propperties of our IdentityMatrix. */
namespace mtl {
namespace ashape {

/** Define IdentityMatri to be a non-scalar type. */
template<>
struct ashape_aux<GraphSymmetricMatrix>{
	typedef nonscal type;
};
}

/** Make comments*/
namespace itl {
  template <class Real, class OStream = std::ostream>
  class visual_iteration : public ::itl::cyclic_iteration <Real>
  {
      typedef ::itl::cyclic_iteration<Real> super;
      typedef visual_iteration self;

    public:

      template <class Vector>
      visual_iteration(GraphType& graph, CME212::SFML_Viewer& viewer, mtl::dense_vector<double>& x_soln, const Vector& r0, int max_iter_, Real tol_, Real atol_ = Real(0), int cycle_ = 100)  
      : super(r0, max_iter_, tol_, atol_, cycle_), graph_(graph), viewer_(viewer), x_soln_(x_soln) {
	node_map_ = viewer_.empty_node_map(graph_);
        viewer_.add_nodes( graph_.node_begin() , graph_.node_end() , ColorFn(x_soln_) , NodePosition(x_soln_), node_map_ );
        viewer_.add_edges( graph_.edge_begin() , graph_.edge_end(), node_map_ );
        viewer_.center_view();
	}
      

      bool finished() { 
	viewer_.add_nodes( graph_.node_begin() , graph_.node_end() , ColorFn(x_soln_) , NodePosition(x_soln_), node_map_ );
	return super::finished(); 
      }
     
      template <typename T>
      bool finished(const T& r){
        viewer_.add_nodes( graph_.node_begin() , graph_.node_end() , ColorFn(x_soln_) , NodePosition(x_soln_), node_map_ );
        //viewer_.add_edges( graph_.edge_begin() , graph_.edge_end(), node_map );
        //viewer_.center_view();
        bool ret= super::finished(r);
        return ret;
      }

    private:
      GraphType& graph_;
      CME212::SFML_Viewer& viewer_;
      std::map<NodeType, unsigned int> node_map_;
      mtl::dense_vector<double>& x_soln_;
  };
}

/** IdentityMatrix implements the Collection concept 
 * with value_type and size_type */
template<>
struct Collection<GraphSymmetricMatrix> {
	typedef double value_type;
	typedef unsigned size_type;
};
} // end namespace


/** MAIN */
int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  {
    // Create a nodes_file from the first input argument
    std::ifstream nodes_file(argv[1]);
    // Interpret each line of the nodes_file as a 3D Point and add to the Graph
    std::vector<NodeType> node_vec;
    Point p;
    while (CME212::getline_parsed(nodes_file, p))
      node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

    // Create a tets_file from the second input argument
    std::ifstream tets_file(argv[2]);
    // Interpret each line of the tets_file as four ints which refer to nodes
    std::array<int,4> t;
    while (CME212::getline_parsed(tets_file, t)) {
      graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
      graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
      graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
      graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
    }
  }

  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());

  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // HW3: YOUR CODE HERE
 
  // Construct the GraphSymmetricMatrix A using the graph
  GraphSymmetricMatrix A(graph);

  // Create an ILU(0) preconditioner
  itl::pc::identity<GraphSymmetricMatrix>        P(A);

  // Define b using the graph, f, and g.
  mtl::dense_vector<double> b_RHS(graph.num_nodes(), 0.0);
  for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit){
    auto i = *nit; // i is my node
    b_RHS[i.index()] = b(i, graph);
  } 

  // Launch the SFML Viewer
  CME212::SFML_Viewer viewer;

  // Set x, Initial Guess
  mtl::dense_vector<double> x_soln(graph.num_nodes(), 0.0);
  
  // cyclic_iteration
  //itl::cyclic_iteration<double> iter(b_RHS, 300, 1.e-10, 0.0, 50);
  
  // visual iteration
  mtl::itl::visual_iteration<double> iter(graph, viewer, x_soln, b_RHS, 300, 1.e-10, 0.0, 50); // graph, viewer, x, b, max_iter, tol
  
  // Solve Ax = b using MTL.
  itl::cg(A, x_soln, b_RHS, P, iter);
	//auto node_map = viewer.empty_node_map(graph);
        //viewer.add_nodes( graph.node_begin() , graph.node_end() , ColorFn(25.0) , NodePosition(x_soln), node_map );
        //viewer.add_edges( graph.edge_begin() , graph.edge_end(), node_map );

  //viewer.center_view();
  //viewer.event_loop();

  return 0;
}
