/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <chrono>
#include <thread>
//#include <cmath>

#include <thrust/for_each.h>
#include <thrust/system/omp/detail/par.h>
#include <thrust/system/omp/execution_policy.h>


#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "SpaceSearcher.hpp"

#include "Graph.hpp"
#include <math.h>


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  //double c; // Node damping coefficient
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double K; // spring constant
  double L; // edge legnth at time 0
  EdgeData(double spring_constant, double init_length){
     K = spring_constant;
     L = init_length;
  }
  EdgeData() {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */

// Functor to update the node positions (i.e. compute t+dt)
struct UpdatePosition {
  void operator()(Node n) {
    n.position() += n.value().vel * dt;
  }
  double dt;
};

template<typename F>
struct UpdateVelocity {
  void operator()(Node n) {
    n.value().vel += force(n, t) * (dt / n.value().mass);
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      n.value().vel = Point(0,0,0);
    }
  }
  double t;
  double dt; 
  F force;
};


template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  // updates all node positions
  

 
  
  //thrust::for_each(thrust::system::omp::par, g.node_begin(), g.node_end(), UpdatePosition{dt});
  thrust::for_each(thrust::system::detail::sequential::seq, g.node_begin(), g.node_end(), UpdatePosition{dt});

 
  
  //for (auto it = g.node_begin(); it != g.node_end(); ++it) {
  //  auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
  //  n.position() += n.value().vel * dt;
  //}

  constraint(g,t); 

  // Compute the t+dt velocity
  // updates all node velocities

 
  //thrust::for_each(thrust::system::omp::par, g.node_begin(), g.node_end(), UpdateVelocity<F>{t, dt, force});
  thrust::for_each(thrust::system::detail::sequential::seq, g.node_begin(), g.node_end(), UpdateVelocity<F>{t, dt, force});

  
  //for (auto it = g.node_begin(); it != g.node_end(); ++it) {
  //  auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
  //  n.value().vel += force(n, t) * (dt / n.value().mass);

  //  if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
	//    n.value().vel = Point(0,0,0);
  //  }
  //}
  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  // Member variables
  // Constructor | Initializer
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    Point force_spring = Point(0,0,0);
    Point force_gravity = n.value().mass*Point(0,0,-grav);
    // Compute x,y,z spring force components 
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
	if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
	    return Point(0,0,0);
        }
        auto e = *it;
	Point diff = n.position() - e.node2().position();
	force_spring += -e.value().K*diff/norm( diff )*( norm(diff) - e.value().L );
    }
    // Compute x,y,z total force components 
    //(void) n; 
    (void) t; // (void) grav;    // silence compiler warnings
    //std::cout << force_spring << std::endl;
    return force_spring + force_gravity;
  }
};

// implements the force of gravity
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force_gravity = n.value().mass*Point(0,0,-grav); 
    (void) t; // silence compiler warnings
    return force_gravity;
  }
};
// implements the spring forces
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force_spring = Point(0,0,0);
    // Compute x,y,z spring force components 
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
        auto e = *it;
	Point diff = n.position() - e.node2().position();
	force_spring += -e.value().K*diff/norm( diff )*( norm(diff) - e.value().L );
    }
    // Compute x,y,z total force components  
    (void) t; // (void) grav;    // silence compiler warnings
    return force_spring;
  }
};
// implements the damping force
struct DampingForce {
  double c; // damping coefficient
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force_damping = -n.value().vel*c; 
    (void) t; // silence compiler warnings
    return force_damping;
  }
};

// implements the damping force
template<typename Force1, typename Force2>
struct CombinedForce {
  Force1 f1;
  Force2 f2;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n,t) + f2(n,t);
  }
};

template<typename Force1, typename Force2>
CombinedForce<Force1, Force2> make_combined_force(Force1 f1, Force2 f2){
 return {f1,f2};
}

template<typename Force1, typename Force2, typename Force3>
CombinedForce<CombinedForce<Force1, Force2>, Force3> make_combined_force(Force1 f1, Force2 f2, Force3 f3){
 return make_combined_force(make_combined_force(f1,f2), f3);
}


struct PlaneConstraint {
  //template<typename GRAPH>
  void operator()(GraphType& graph, double t){
	(void) t;
	//std::cout << "plane constraint" << std::endl;
  	for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
		auto n = *it;
		if (dot(n.position(), Point(0,0,1)) < -0.75){
			n.value().vel.z = 0.0;
			n.position().z = -0.75;
		}
	} 
  }
};

struct SphereConstraint1 {
  //template<typename GRAPH>
  double radius = 0.15;
  Point sphere_center = Point(0.5,0.5,-0.5);
  void operator()(GraphType& graph, double t){
	(void) t;
	for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
		auto n = *it;
		Point dist = n.position() - sphere_center;
		double dis = norm(dist);
		if (dis == 0){
			n.position() = n.position() + Point(0,0,0.15);
			n.value().vel.z = 0;
		}
		else if (dis < radius) {
			n.position() = sphere_center + (dist*(radius/dis));
			n.value().vel = n.value().vel - dot(n.value().vel, dist/dis)*(dist/dis);
		}
	}
  }
};

struct SphereConstraint2 {
  //template<typename GRAPH>
  void operator()(GraphType& graph, double t){
    (void) t;
    //std::cout << "sphere constraint 2" << std::endl;
    double radius = 0.15;
    Point sphere_center = Point(0.5,0.5,-0.5);
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
      auto n = *it;
      Point diff = n.position() - sphere_center;
      if ( norm( diff ) < radius ) {
        graph.remove_node(n);
      }
    }
  }
};

template <typename Constraint1, typename Constraint2>//, typename Constraint3>
struct CombinedConstraints {
  Constraint1 c1;
  Constraint2 c2;
  //Constraint3 c3;

  void operator()(GraphType& graph, double t){
    c1(graph,t);
    c2(graph,t);
    //c3(graph,t);
  }
};

template<typename Constraint1, typename Constraint2>//, typename Constraint3>
CombinedConstraints<Constraint1, Constraint2> make_combined_constraints(Constraint1 c1, Constraint2 c2){
 return {c1,c2};
}

template<typename Constraint1, typename Constraint2, typename Constraint3>
CombinedConstraints<CombinedConstraints<Constraint1, Constraint2>, Constraint3> make_combined_constraints(Constraint1 c1, Constraint2 c2, Constraint3 c3){
 return make_combined_constraints(make_combined_constraints(c1,c2), c3);
}

struct SelfCollisionConstraint {
    void operator()(GraphType& g, double t) const {
        Box3D bigbb(Point(2,2,2), Point(-2,-2,-2));
        auto n2p = [] (const Node& n) {return n.position();};
	SpaceSearcher<Node> ss(bigbb, g.node_begin(), g.node_end(), n2p);
	for (auto it=g.node_begin(); it != g.node_end(); ++it){
	    auto n = *it;
            const Point& center = n.position();
	    double radius2 = std::numeric_limits<double>::max();
	    for (auto eit = n.edge_begin(); eit != n.edge_end(); ++eit){
		auto e = *eit;
		radius2 = std::min(radius2, normSq(e.node2().position()-center));
            }
	    radius2 *= 0.9;
	    Box3D bb(center + sqrt(radius2), center - sqrt(radius2));
	    for (auto nhi = ss.begin(bigbb); nhi != ss.end(bigbb); ++nhi){
                Node n2 = *nhi;
		Point r = center - n2.position();
		double l2 = normSq(r);
		if (n != n2 && l2 < radius2) {
		    n.value().vel -= (dot(r,n.value().vel)/l2)*r;
		}
	    } 
        }
    }
};




int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }
  
  // Construct an empty graph
  GraphType graph;
  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));
  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);
//#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }
  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  double spring_constant = 100.0;
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().vel = Point(0,0,0);
    n.value().mass = 1.0/graph.num_nodes();
  }
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    auto e = *it;
    e.value() = EdgeData(spring_constant, e.length());
  }

  double c = 1.0/graph.num_nodes(); // damping coefficient
  DampingForce force_damping{c};
  MassSpringForce force_spring;
  GravityForce force_gravity;
  auto total_force = make_combined_force(force_damping, force_spring, force_gravity);


  PlaneConstraint constraint1;
  SphereConstraint1 constraint2;
  SphereConstraint2 constraint3;
  SelfCollisionConstraint constraint4 = SelfCollisionConstraint();

  auto all_constraints = make_combined_constraints(constraint1, constraint3, constraint2);
  
  

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 0.001; 
      double t_start = 0;
      double t_end = 5.0;
      CME212::Clock clock;
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, total_force, all_constraints); // replace Problem1Force() with CombinedForce()
	
	// clear the viewer's nodes and edges 
	viewer.clear();
	node_map.clear();	
	
        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
	viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        
	viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 100)
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }
      std::cout << clock.seconds() << std::endl;

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}

