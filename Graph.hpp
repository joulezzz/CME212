#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V, typename E>
class Graph {
private:

public:

    //
    // PUBLIC TYPE DEFINITIONS
    //

    /** Type of this graph. */
    using graph_type = Graph;

    typedef V node_value_type;
    typedef E edge_value_type;

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

    using uid_type = size_type;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() {

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
    class Node : private totally_ordered<Node> {
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
        const Point &position() const {
            assert(valid());
            return graph_->nodes_[uid_].position_;
        }

        /** Return this node's position. */
        Point &position() {
            assert(valid());
            return graph_->nodes_[uid_].position_;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            assert(valid());
            return graph_->nodes_[uid_].idx_;
        }

        /** Return this node's value. */
        const node_value_type &value() const {
            assert(valid());
            return graph_->nodes_[uid_].value_;
        }

        /** Return this node's value. */
        node_value_type &value() {
            assert(valid());
            return graph_->nodes_[uid_].value_;
        }

        /** Return this node's degree. */
        size_type degree() const {
            assert(valid());
            return graph_->adjacency_[uid_].size();
        }

        /** Return the beginning of this node's incident iterator. */
        incident_iterator edge_begin() const {
            assert(valid());
            return IncidentIterator(graph_, uid_, 0);
        }

        /** Return the ending of this node's incident iterator. */
        incident_iterator edge_end() const {
            assert(valid());
            return IncidentIterator(graph_, uid_, degree());
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node &n) const {
            assert(valid());
            return (graph_ == n.graph_) and (uid_ == n.uid_);
        }

        /** Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        bool operator<(const Node &n) const {
            assert(valid());
            if (graph_ == n.graph_) {
                return uid_ < n.uid_;
            }
            return (graph_ < n.graph_);
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        Graph *graph_;
        size_type uid_;

        /** Valid Node constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a uid_ < graph_->num_nodes()
         */
        Node(const Graph *graph, size_type uid)
                : graph_(const_cast<Graph *>(graph)), uid_(uid) {
            assert(graph_ != NULL);
            assert(0 <= uid_ and uid_ < graph_->nodes_.size());
        }

        /** Check the representation invariant of the node
         * @return true if node is a valid node.
         *         false if node is invalid.
         */
        bool valid() const {
            return uid_ >= 0 and uid_ < graph_->nodes_.size()
                   and graph_->nodes_[uid_].idx_ < graph_->i2u_.size()
                   and graph_->i2u_[graph_->nodes_[uid_].idx_] == uid_;
        }
    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return i2u_.size();
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }

    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point &position, const node_value_type &node_value = node_value_type()) {
        uid_type  uid = nodes_.size();
        size_type idx = i2u_.size();
        node_info new_node (idx, position, node_value);
        nodes_.push_back(new_node);
        i2u_.push_back(uid);
        std::vector <edge_info> adjacent_edges;
        adjacency_.push_back(adjacent_edges);
        return Node(this, uid);
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node &n) const {
        return (this == n.graph_) and (i2u_[n.index()] == n.uid_);
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        assert(0 <= i and i < num_nodes());
        return Node(this, i2u_[i]);
    }

    /** Remove a node from the graph.
     * @param[in] n The valid node to remove
     * @return The index of the removed node
     *
     * @pre @a n is valid node
     * @pre graph.node(i).index() == i for all 0 <= i <= num_nodes()
     * @pre graph.node(@a n.index()) == @a n
     * @pre has_node(@a n) == true
     * @post @a n is invalidated
     * @post new num_nodes() == old num_nodes() - 1
     * @post has_node(@a n) == false
     * @post For all edges e that adjacent to @a n
     *       ie. e.node1() == @a n or e.node2() == @a n,
     *       e becomes invalid
     * @post new num_edges() == old num_edges() - n.degree()
     * @post Outstanding NodeIterators and EdgeIterators are invalidated
     * @post Outstanding IncidentIterators associated with @a n is invalidated
     *
     * Complexity: O(num_nodes)
     */
    size_type remove_node(const Node &n) {

        // Remove adjacent edges
        for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
            Edge e = *iit;
            Node head = e.node2();
            remove_tail_node(head, n);
        }
        adjacency_[n.uid_].clear();

        // Reindex the idx of the nodes followed the removed node
        // in nodes vector
        size_type idx = n.index();
        for (size_type k = idx; k < i2u_.size(); ++k) {
            --nodes_[i2u_[k]].idx_;
        }

        // Erase the node from i2u
        i2u_.erase(i2u_.begin() + idx);

        return idx;
    }

    /** Remove a node from the graph.
     * @param[in] ni The valid node iterator to remove
     * @return NodeIterator pointing to node followed the erased node
     *
     * @pre @a n is valid node
     * @pre graph.node(i).index() == i for all 0 <= i <= num_nodes()
     * @pre graph.node(@a n.index()) == @a n
     * @pre has_node(@a n) == true
     * @post @a n becomes invalid node
     * @post new num_nodes() == old num_nodes() - 1
     * @post has_node(@a n) == false
     * @post For all edges e that adjacent to @a n
     *       ie. e.node1() == @a n or e.node2() == @a n,
     *       e becomes invalid.
     * @post new num_edges() == old num_edges() - n.degree()
     * @post Outstanding NodeIterators and EdgeIterators are invalidated.
     * @post Outstanding IncidentIterators associated with @a n is invalidated.
     *
     * Complexity: O(num_nodes)
     */
    node_iterator remove_node(NodeIterator ni) {
        Node n = *ni;
        size_type idx = remove_node(n);
        return NodeIterator(this, idx);
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
    class Edge : private totally_ordered<Edge> {
    public:
        /** Construct an invalid Edge. */
        Edge() {

        }

        /** Return a node of this Edge */
        Node node1() const {
            return Node(graph_, node1_uid_);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return Node(graph_, node2_uid_);
        }

        /** Return the length of this Edge */
        double length() const {
            return norm(node1().position() - node2().position());
        }

        /** Return the value of this Edge */
        edge_value_type &value() {
            size_type index;
            if (node1_uid_ > node2_uid_) {
                index = graph_->match_tail_node(Node(graph_, node2_uid_), Node(graph_, node1_uid_));
                return graph_->adjacency_[node2_uid_][index].value_;
            }
            index = graph_->match_tail_node(Node(graph_, node1_uid_), Node(graph_, node2_uid_));
            return graph_->adjacency_[node1_uid_][index].value_;
        }

        /** Return the value of this Edge */
        const edge_value_type &value() const {
            return value();
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge &e) const {
            bool same_graph = (graph_ == e.graph_);
            bool direction1 = ((node1_uid_ == e.node1_uid_) and (node2_uid_ == e.node2_uid_));
            bool direction2 = ((node1_uid_ == e.node2_uid_) and (node2_uid_ == e.node1_uid_));
            return same_graph and (direction1 or direction2);
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge &e) const {
            if (graph_ == e.graph_) {
                if (node1_uid_ == e.node1_uid_) {
                    return (node2_uid_ < e.node2_uid_);
                }
                return (node1_uid_ < e.node1_uid_);
            }
            return (graph_ < e.graph_);
        }

    private:
        friend class Graph;

        Graph *graph_;
        size_type node1_uid_;
        size_type node2_uid_;

        /** Valid Edge constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a node1_uid_, @a node2_uid_ < graph_->num_nodes()
         * @pre @a node1_uid_ != @a node2_uid_
         */
        Edge(const Graph *graph, size_type node1_uid, size_type node2_uid)
                : graph_(const_cast<Graph *>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {
            assert(graph_ != NULL);
            assert(node1_uid_ != node2_uid_);
            assert(0 <= node1_uid_ and node1_uid_ < graph_->nodes_.size());
            assert(0 <= node2_uid_ and node2_uid_ < graph_->nodes_.size());
        }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: O(1)
     */
    size_type num_edges() const {
        return std::distance(edge_begin(), edge_end());
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: O(num_nodes() + num_edges())
     */
    Edge edge(size_type i) const {
        assert(0 <= i and i < num_edges());
        return *std::next(edge_begin(), i);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: O(degree())
     */
    bool has_edge(const Node &a, const Node &b) const {
        size_type index = match_tail_node(a, b);
        return index != a.degree();
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
     * Complexity: O(degree())
     */
    Edge add_edge(const Node &a, const Node &b) {
        size_type node1_uid = a.uid_;
        size_type node2_uid = b.uid_;
        if (!has_edge(a, b)) {
            edge_info new_edge1(node2_uid, edge_value_type());
            edge_info new_edge2(node1_uid, edge_value_type());
            adjacency_[node1_uid].push_back(new_edge1);
            adjacency_[node2_uid].push_back(new_edge2);
        }
        return Edge(this, node1_uid, node2_uid);
    }

    /** Remove an edge from the graph.
     * @param[in] n1 The valid node at one end of the edge to remove
     * @param[in] n2 The valid node at one end of the edge to remove
     * @return The column index on the adjacency list where the edge
     *         is removed if has_edge(@a n1, @a n2) is true.
     *         Otherwise, return the degree of node which has smaller
     *         uid.
     *
     * @pre @a n1 and @a n2 are valid node
     * @post new num_edges() == old num_edges() - 1
     * @post has_edge(@a n1, @a n2) == false
     * @post Outstanding EdgeIterators are invalidated.
     * @post Outstanding IncidentIterators associated with @a n1 and
     *       @a n2 are invalidated.
     *
     * Complexity: O(degree())
     */
    size_type remove_edge(const Node &n1, const Node &n2) {
        if (!has_edge(n1, n2)) {
            return (n1.uid_ < n2.uid_) ? n1.degree() : n2.degree();
        }
        size_type idx1 = remove_tail_node(n1, n2);
        size_type idx2 = remove_tail_node(n2, n1);
        return (n1.uid_ < n2.uid_) ? idx1 : idx2 ;
    }

    /** Remove an edge from the graph.
     * @param[in] e The valid edge to remove
     * @return The column index on the adjacency list where the edge
     *         is removed if has_edge(@a e) is true.
     *         Otherwise, return the degree of the adjacency list that
     *         associated with @a e.
     *
     * @pre @a e is valid edge
     * @post new num_edges() == old num_edges() - 1
     * @post has_edge(@a n1, @a n2) == false
     * @post Outstanding EdgeIterators are invalidated.
     * @post Outstanding IncidentIterators associated with @a e.node1() and
     *       @a e.node2() are invalidated.
     *
     * Complexity: O(degree())
     */
    size_type remove_edge(const Edge &e) {
        return remove_edge(e.node1(), e.node2());
    }

    /** Remove an edge from the graph.
     * @param[in] ei The valid EdgeIterator to remove
     * @return The EdgeIterator pointing to the edge after the removed edge
     *
     * @pre @a e is valid edge
     * @post new num_edges() == old num_edges() - 1
     * @post has_edge(@a n1, @a n2) == false
     * @post Outstanding EdgeIterators are invalidated.
     * @post Outstanding IncidentIterators associated with @a (*ei).node1() and
     *       @a (*ei).node2() are invalidated.
     *
     * Complexity: O(degree())
     */
    edge_iterator remove_edge(edge_iterator ei) {
        Edge e = *ei;
        size_type adj_col = remove_edge(e.node1().uid_, e.node2().uid_);
        return EdgeIterator(this, e.node1().uid_, adj_col);
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     *
     * Complexity: O(num_nodes())
     */
    void clear() {
        nodes_.clear();
        adjacency_.clear();
        i2u_.clear();
    }

    //
    // Node Iterator
    //
/**
    // @class Graph::NodeIterator
    //  @brief Iterator class for nodes. A forward iterator. 
    class NodeIterator : private totally_ordered<NodeIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     // Element type
        using pointer           = Node *;                    // Pointers to elements
        using reference         = Node &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        // Construct an invalid NodeIterator. 
        NodeIterator() {
        }

        /** Return the node at the node iterator.
         * @pre node_uid_ != graph_->num_nodes()
         
        Node operator*() const {
            assert(idx_ != graph_->num_nodes());
            return graph_->node(idx_);
        }

        /** Return the node iterator to the next element. 
        NodeIterator &operator++() {
            ++idx_;
            return *this;
        }

        /** Test whether this node iterator and @a ni are equal.
         *
         * Equal node iterators have the same graph and the same node id.
         
        bool operator==(const NodeIterator &ni) const {
            return (graph_ == ni.graph_) and (idx_ == ni.idx_);
        }

    private:
        friend class Graph;

        Graph *graph_;
        size_type idx_;

        /** Valid NodeIterator constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a node_uid_ < graph_->num_nodes()
         
        NodeIterator(const Graph *graph, size_type idx)
                : graph_(const_cast<Graph *>(graph)), idx_(idx) {
            assert(graph_ != NULL);
            assert(0 <= idx and idx_ <= graph_->num_nodes());
        }

    };
*/

    /** Return a node iterator pointing to the first element of node sequence.
     * @return NodeIterator pointing to the first element
     * @post If node sequence is empty, dereference is invalid for the return value.
     *
     * Complexity: O(1)
     */
    node_iterator node_begin() const {
        return NodeIterator(this, 0);
    }

    /** Return a node iterator pointing to the past-the-end element in node sequence.
     * @return NodeIterator pointing to the past-the-end element
     * @post If node sequence is empty, dereference is invalid for the return value.
     *
     * Complexity: O(1)
     */
    node_iterator node_end() const {
        return NodeIterator(this, num_nodes());
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
        using pointer           = Edge *;                    // Pointers to elements
        using reference         = Edge &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator() {
        }

        /** Return the edge at the incident iterator.
         * @pre adj_row_ != graph_->num_nodes()
         */
        Edge operator*() {
            assert(adj_row_ != graph_->nodes_.size());
            size_type node1_uid = adj_row_;
            size_type node2_uid = graph_->adjacency_[adj_row_][adj_col_].uid2_;
            return Edge(graph_, node1_uid, node2_uid);
        }

        /** Return the incident iterator to the next element. */
        IncidentIterator &operator++() {
            ++adj_col_;
            return *this;
        }

        /** Test whether this incident iterator and @a iit are equal.
         *
         * Equal incident iterators have the same graph and the same row
         * and column on the adjacency list.
         */
        bool operator==(const IncidentIterator &iit) const {
            return (graph_ == iit.graph_) and (adj_row_ == iit.adj_row_) and (adj_col_ == iit.adj_col_);
        }

    private:
        friend class Graph;

        Graph *graph_;
        size_type adj_row_;
        size_type adj_col_;

        /** Valid IncidentIterator constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a adj_row_ <= graph_->num_nodes()
         * @pre 0 <= @a adj_col_ <= graph_->adjacency_[adj_row_].size()
         */
        IncidentIterator(const Graph *graph, size_type adj_row, size_type adj_col)
                : graph_(const_cast<Graph *>(graph)), adj_row_(adj_row), adj_col_(adj_col) {
            assert(graph_ != NULL);
            assert(0 <= adj_row_ and adj_row_ <= graph_->nodes_.size());
            assert(0 <= adj_col_ and adj_col_ <= graph_->adjacency_[adj_row_].size());
        }
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
        using pointer           = Edge *;                    // Pointers to elements
        using reference         = Edge &;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid EdgeIterator. */
        EdgeIterator() {
        }

        /** Return the edge at the edge iterator.
         * @pre adj_row_ != graph_->num_nodes()
         */
        Edge operator*() const {
            assert(adj_row_ != graph_->nodes_.size());
            size_type node1_uid = adj_row_;
            size_type node2_uid = graph_->adjacency_[adj_row_][adj_col_].uid2_;
            return Edge(graph_, node1_uid, node2_uid);
        }

        /** Return the edge iterator to the next element. */
        EdgeIterator &operator++() {
            ++adj_col_;
            fix();
            return *this;
        }

        /** Test whether this edge iterator and @a eit are equal.
         *
         * Equal edge iterators have the same graph and the same row
         * and column on the adjacency list.
         */
        bool operator==(const EdgeIterator &eit) const {
            return (graph_ == eit.graph_) and (adj_row_ == eit.adj_row_) and (adj_col_ == eit.adj_col_);
        }

    private:
        friend class Graph;

        Graph *graph_;
        size_type adj_row_;
        size_type adj_col_;

        /** Valid EdgeIterator constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a adj_row_ <= graph_->num_nodes()
         * @pre 0 <= @a adj_col_ <= graph_->adjacency_[adj_row_].size()
         * @post The iterator will start at the first valid Edge
         */
        EdgeIterator(const Graph *graph, size_type adj_row)
                : graph_(const_cast<Graph *>(graph)), adj_row_(adj_row), adj_col_(0) {
            assert(graph_ != NULL);
            assert(0 <= adj_row_ and adj_row_ <= graph_->nodes_.size());
            assert(0 <= adj_col_ and adj_col_ <= graph_->adjacency_[adj_row_].size());
            fix();
        }

        EdgeIterator(const Graph *graph, size_type adj_row, size_type adj_col)
                : graph_(const_cast<Graph *>(graph)), adj_row_(adj_row), adj_col_(adj_col) {
            assert(graph_ != NULL);
            assert(0 <= adj_row_ and adj_row_ <= graph_->nodes_.size());
            assert(0 <= adj_col_ and adj_col_ <= graph_->adjacency_[adj_row_].size());
            fix();
        }

        /** Move the edge iterator to the next valid element.
         *  (Amy helped me with this. Thanks Amy!)
         */
        void fix() {
            /* Move the row index of adjacency list */
            while (adj_row_ < graph_->adjacency_.size()) {
                /* Move the column index of adjacency list */
                while (adj_col_ < graph_->adjacency_[adj_row_].size()) {
                    /* Edge is valid if node1_uid < node2_uid */
                    size_type node1_uid = adj_row_;
                    size_type node2_uid = graph_->adjacency_[adj_row_][adj_col_].uid2_;
                    if (node1_uid < node2_uid) {
                        return;
                    }
                    ++adj_col_;
                }
                ++adj_row_;
                adj_col_ = 0;
            }
            return;
        }

    };

    /** Return an edge iterator pointing to the first element of edge sequence.
     * @return EdgeIterator pointing to the first element
     * @post If edge sequence is empty, dereference is invalid for the return value.
     *
     * Complexity: O(1)
     */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }

    /** Return an edge iterator pointing to the past-the-end element in edge sequence.
     * @return EdgeIterator pointing to the past-the-end element
     * @post If edge sequence is empty, dereference is invalid for the return value.
     *
     * Complexity: O(1)
     */
    edge_iterator edge_end() const {
        return EdgeIterator(this, adjacency_.size());
    }

private:

    /** Erase the giver index from the vector
     *  @param[in] v The vector to delete the element from
     *  @param[in] i The index of the element to delete
     *
     *  The element at @a i is replaced with the last
     *  element in the vector.
     *
     *  Complexity: O(1)
     */
    template<typename T>
    void erase(std::vector <T> &v, unsigned int i) {
        v[i] = v.back();
        v.pop_back();
    }

    /** Find the adjacent node of head node that have matching
     *  uid with tail node
     *  @param[in] head The valid node that the function will
     *                  access its adjacency list
     *  @param[in] tail The valid node to find
     *  @return The column index on the adjacency list where the
     *          tail node is located, if the tail node is found.
     *          Otherwise, return the degree.
     *
     *  @pre @a head and @a tail are valid nodes.
     *
     *  Complexity: O(degree())
     */
    size_type match_tail_node(const Node &head, const Node &tail) const {
        size_type match_uid = tail.uid_;
        auto is_match = [&match_uid](Edge e) { return e.node2().uid_ == match_uid; };
        IncidentIterator iit = std::find_if(head.edge_begin(), head.edge_end(), is_match);
        return std::distance(head.edge_begin(), iit);
    }

    /** Remove the adjacent node of head node that have matching
    *  uid with tail node
    *  @param[in] head The valid node that the function will
    *                  access its adjacency list
    *  @param[in] tail The valid node to remove
    *  @return The column index on the adjacency list where the
    *          node is deleted, if the tail node is found.
    *          Otherwise, return the degree.
    *
    *  Complexity: O(degree())
    */
    size_type remove_tail_node(const Node &head, const Node &tail) {
        size_type index = match_tail_node(head, tail);
        if (index != head.degree()) {
            erase(adjacency_[head.uid_], index);
        }
        return index;
    }

    /** Structure representing node information */
    struct node_info {
        size_type idx_;
        Point position_;
        node_value_type value_;

        node_info(size_type idx, const Point& position, const node_value_type& value)
                :idx_(idx), position_(position), value_(value){

        }
    };
    /** Structure representing edge information */
    struct edge_info {
        size_type uid2_;
        edge_value_type value_;

        edge_info(size_type uid2, const edge_value_type& value)
        :uid2_(uid2), value_(value){

        }
    };

    std::vector <node_info> nodes_;
    std::vector <uid_type> i2u_;
    std::vector <std::vector<edge_info>> adjacency_;

};

#endif // CME212_GRAPH_HPP

