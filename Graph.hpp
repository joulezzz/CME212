
#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
// Allows Nodes to support user specified value of type node_value_type: Graph<V>
template <typename V>
class Graph {
 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of the node values */
  using node_value_type = V;

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes(),edges(),adjacency(){
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes[uid_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /**
     * @return The value of this node which is node_value_type
	 */
    node_value_type & value (){
      return graph_->nodes[uid_].second;
    }
    
   	/**
     * @return The value of this node which is const node_value_type
	 */
    const node_value_type & value () const {
      return graph_->nodes[uid_].second;
    }

   	/**
     * @return The number of nodes adjacent to this node
	 */
    size_type degree() const {
      return graph_->adjacency[uid_].size();

    }

    /**
     * @return An IncidentIterator that points to an edge adjacent to
     * 		   the this node that is the first of the sequence of edges adjacent
     *         to this node
	 */
    IncidentIterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);

    }

    /**
     * @return An IncidentIterator that does not point to an ajacent edge
     * indicating that all adjacent edges have been iterated through
	 */
    IncidentIterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (graph_ == n.graph_) // recall these are pointers so this check their address
        if (uid_ ==  n.uid_)  // this is not a pointer, so it compares numbers
          return true;
      return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      if (graph_ == n.graph_){ // if this belongs to same graph
        if (uid_ <  n.uid_)
          return true;
      }
      else{
        if (graph_ < n.graph_) // if this belongs to previous graph
          return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph_pointer, size_type uid)
        : graph_( const_cast<Graph*>(graph_pointer) ), uid_( uid ) {
    }
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size(); // recall nodes is a vector, size is its method
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size(); // just calls the same function above
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] new_value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& new_value = node_value_type()) {
    //nodes.push_back(position);
    nodes.push_back(std::make_pair(position, new_value));
    adjacency.push_back(std::vector<std::pair<size_type, size_type>> ()); 
    return Node(this,nodes.size()-1); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this)
      if (n.uid_<nodes.size())
        return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_,node1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,node2_id);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ == e.graph_){
        if ((node1_id == e.node1_id and node2_id == e.node2_id) or 
             (node1_id == e.node2_id and node2_id == e.node1_id)){
          return true;
        }
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_){
        if (uid_ < e.uid_)
          return true;
      }
      else{
        if (graph_ < e.graph_)
          return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type uid_;
    size_type node1_id;
    size_type node2_id;
    Edge(const Graph* graph_pointer, size_type uid, size_type node1, size_type node2)
        : graph_(const_cast<Graph*>(graph_pointer)), uid_(uid), node1_id(node1), node2_id(node2){
    }
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i, edges[i].first, edges[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned int i=0; i<adjacency[a.uid_].size(); i++)
       if (b.uid_ == adjacency[a.uid_][i].first)
         return true;
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    for (unsigned int i=0; i<adjacency[a.uid_].size(); i++)
       if (b.uid_ == adjacency[a.uid_][i].first)
         return Edge(this, adjacency[a.uid_][i].second, a.uid_, b.uid_);
    if (a<b)
      edges.push_back(std::make_pair(a.uid_,b.uid_));
    else
      edges.push_back(std::make_pair(b.uid_,a.uid_));
    adjacency[a.uid_].push_back(std::make_pair(b.uid_,edges.size()-1));
    adjacency[b.uid_].push_back(std::make_pair(a.uid_,edges.size()-1));
    return Edge(this,edges.size()-1,a.uid_,b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
    adjacency.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {

   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /**
     * @param[in] g The graph this NodeIterator points to
     * @param[in] i The index associated with this node
     * @return NodeIterator 
     * @pre 0 <= @a i < num_nodes()
     */
    NodeIterator(const Graph* g, size_type i): graph_(g), index_points(i) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** "Dereferences the NodeIterator"
     * @return Node which points to the graph @graph_ and has @index_points as its index
     * @post 0 <= @aindex_points < num_nodes()
     */
    Node operator*() const {
    	return Node(graph_ , index_points);
    }

    /** Increments the NodeIterator to point to the next Node or nullptr
     * @return NodeIterator that points to exactly one of the following:
     *				1.) if index_points < nodes.size(), next Node 
     *				2.) if index_points = nodes.size(), nullptr 
     */
    NodeIterator& operator++(){
      index_points++;
      return  *this;

    }

    /** Checks if two NodeIterators are the same
     * @param[in] nodeIterator A NodeIterator 
     * @return bool value
     * @post bool is true if the following are all true
     *         1. nodeIterator and this point to the same graph
     *         2. if index_points member vairable of nodeIterator and this have same value
     */
    bool operator==(const NodeIterator& nodeIterator) const {
      return ( nodeIterator.graph_ == graph_) && (nodeIterator.index_points == index_points);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph *graph_;
    size_type index_points;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

 /** returns a NodeIterator that points to the first Node in the sets of nodes of this graph
  * @return NodeIterator
  * @post NodeIterator will point to the first Node in this graph if graph is nonempty
  *       otherwise it will return a nullptr
  */
  NodeIterator node_begin() const {
  	return NodeIterator(this, 0);
  }

 /** returns a NodeIterator that indicates the end of the nodes
  * @return NodeIterator
  * @post NodeIterator is assigned nullptr
  */
  NodeIterator node_end() const {
    return NodeIterator(this, nodes.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Constructor for Incident iterator 
     * @param[in] graph An object of Graph
     * @param[in] node_index A value of size_type
     * @param[in] edge_index A value of size_type
     * @return IncidentIterator 
     * @pre 0 <= @a node_index < num_nodes()
     * @pre 0 <= @a edge_index < degree()
     * @post IncidentIterator belongs to @a graph, is for @a node_index, and points to @a edge_index
     */
    IncidentIterator(Graph* graph, size_type node_index, size_type edge_index) : graph_(graph), node_index_(node_index), edge_index_(edge_index) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** 
     * @return The Edge that the this points to
     */
    Edge operator*() const {
      return Edge(graph_, graph_->adjacency[node_index_][edge_index_].second, node_index_, graph_->adjacency[node_index_][edge_index_].first);
    }

    /** Increments the IncidentIterator to point to the next adjacent node or nullptr
     * @return NodeIterator that points to exactly one of the following:
     *				1.) if edge_index_ < degree(), next adjacent node 
     *				2.) if edge_Index_ = degree(), nullptr 
     */
    IncidentIterator& operator++() {
      edge_index_++;
      return *this;
    }

    /** Checks if two IncidentIterators are the same
     * @param[in] incidentIterator A IncidentIterator 
     * @return bool value
     * @post bool is true if the following are all true
     *         1. incidentIterator and this belong to the same graph
     *         2. if node_index_ member variable of this and nodeIterator have same value
     *         3. if edge_index_ member variable of this and nodeIterator have same value
     */
    bool operator==(const IncidentIterator& incidentIterator) const {
      return ( incidentIterator.graph_ == graph_) && (incidentIterator.node_index_ == node_index_) && (incidentIterator.edge_index_ == edge_index_);
    }

   private:
    friend class Graph;
    const Graph *graph_;
    size_type node_index_;
    size_type edge_index_; 


    // HW1 #3: YOUR CODE HERE

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Constructor for EdgeIterator
     * @param[in] graph An object of Graph
     * @param[in] edge_index A value of size_type
     * @return EdgeIterator 
     * @pre 0 <= @a edge_index < num_edges()
     */
    EdgeIterator(const Graph* graph, size_type edge_index) : graph_(graph), edge_index_(edge_index) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

     /** "Dereferences the EdgeIterator"
      * @return Edge that this points 
      */
    Edge operator*() const {
      return Edge(graph_, edge_index_, graph_->edges[edge_index_].first, graph_->edges[edge_index_].second);
    }

     /** Increments the EdgeIterator to point to the next Edge (in a global sense) or the nullptr
      * @return NodeIterator that points to exactly one of the following:
      *				1.) if Edge global number < num_edges.size(), next Node 
      *				2.) if Edge global number = num_edges(), nullptr 
      */
    EdgeIterator& operator++() {
      edge_index_++;
      return *this;
    }

     /** Checks if two EdgeIterators are the same
      * @param[in] edgeIterator An object of EdgeIterator 
      * @return bool value
      * @post bool is true if the following are all hold
      *         1. nodeIterator and this belong to the same graph
      *         2. if *edgeIterator == *this (they point to same edge)
      */
    bool operator==(const EdgeIterator& edge_iterator) const {
      if (edge_iterator.graph_ == graph_){
        return (edge_iterator.edge_index_ == edge_index_);
      }
      return false;
    }
  

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_;
    size_type edge_index_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** returns a EdgeIterator that points to the first Edge in the "set" of all edges of this graph
   * @return EdgeIterator
   * @post EdgeIterator will point to the first Edge in this graph if graph is nonempty and has at least two nodes who are adjacent to each other
   *       otherwise it will return a nullptr
   */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);

  }

  /** returns an EdgeIterator that indicates the end of the "set" of edges
   * @return EdgeIterator
   * @post EdgeIterator is assigned nullptr
   */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, edges.size());

  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  // * @nodes is a list with the position of node @i in index i
  //std::vector<Point> nodes;
  std::vector<std::pair<Point,V>> nodes;
  // * @edges is a list containing pairs of nodes for each edge.
  //   The pair stored in position i corresponds to edge i.
  std::vector<std::pair<size_type, size_type>> edges;
  // * @adjacency is a vector containing vector with pairs. It contains
  //   at position i the pairs of each node adjacent to node i and the
  //   the edge id that connect node i to the adjacent node.
  std::vector<std::vector<std::pair<size_type, size_type>>> adjacency;
};

#endif // CME212_GRAPH_HPP
