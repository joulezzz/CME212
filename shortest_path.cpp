/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"
#include <queue>

// Define our types
using GraphType = Graph<int>;
using NodeType  = typename GraphType::node_type;
using NodeIter  = typename GraphType::node_iterator;

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */

NodeIter nearest_node(const GraphType& g, const Point& point)
{
  // HW1 #3: YOUR CODE HERE
  NodeIter min_ni = g.node_begin();
  double minDist = norm_2((*min_ni).position() - point);
  for (auto ni = g.node_begin(); ni != g.node_end(); ++ni){
    auto node = *ni;
    double tempDist = norm_2(node.position() - point);
    if (tempDist < minDist){
      minDist  = tempDist;
      min_ni = ni;
    }
  }
  return min_ni;
}

// may not use
void normalize_nodes(const GraphType& g, int normalizer){
  for (auto ni = g.node_begin(); ni != g.node_end(); ++ni){
    node = *ni;
    node.value() = node.value() / normalizer
  }
}

/** Update a graph with the shortest path lengths from a root node.
 * @param[in,out] g     Input graph
 * @param[in,out] root  Root node to start the search.
 * @return The maximum path length found.
 *
 * @post root.value() == 0
 * @post Graph has modified node values indicating the minimum path length
 *           to the root.
 * @post Graph nodes that are unreachable from the root have value() == -1.
 *
 * This sets all nodes' value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(GraphType& g, NodeType& root)
{
  // HW1 #3: YOUR CODE HERE
  //(void) g, (void) root;      // Quiet compiler warnings

  // initialize all node values to -1
  for (auto ni = g.node_begin(); ni != g.node_end(); ++ni){
    auto node1 = *ni;
    node1.value() = -1;
  }
  int current_max = -1;
  std::queue<NodeType> q;
  root.value() = 0;
  q.push(root);

  while(q.size()){
    auto n = q.front();
    q.pop();
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei ){
      auto e = *ei;
      auto node2 = e.node2();
      if (node2.value() == -1){
        node2.value() = n.value() + 1;
        current_max = node2.value();
      }
    }
  }
  return current_max;
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SFML_Viewer
  CME212::SFML_Viewer viewer;

  // HW1 #3: YOUR CODE HERE
  // Use nearest_node and shortest_path_lengths to set the node values
  // Construct a Color functor and view with the SFML_Viewer
  Point point(-1,0,1);
  NodeIter closest_iter = nearest_node(graph, point);
  root = *closest_iter;
  int longest_path = shortest_path_lengths(graph, root);
  auto node_map = empty_node_map(graph)
  // Working on this
  template <typename T> 
  struct ColorFn {
    const T& longest_path_;
    float operator() (const T& ) const {
      float ratio = float(node_value) /  float(longest_path);
      return Color::make_heat(ratio)
    }
  }


  void add_nodes ( graph.node_begin() , graph.node_end() , color_fn , node_map );




  // Center the view and enter the event loop for interactivity
  viewer.center_view();
  viewer.event_loop();

  return 0;
}
