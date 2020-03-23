#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>

//Forward declarations for classes below so they can be used below without worrying too much about the ordering.
template <typename vertex> class vertex_iterator;
template <typename vertex> class neighbour_iterator;
template <typename vertex> class directed_graph;


template <typename vertex>
class directed_graph {

private:

  //You will need to add some data members here
  //to actually represent the graph internally,
  //and keep track of whatever you need to.
  
  //A map of sets
  //The verticies are the key for the map
  //The set for each key contains the edges of that vertex
  std::map<vertex, std::set<vertex>> map_of_vert;
  


public:


  directed_graph(); //A constructor for directed_graph. The graph should start empty.
  ~directed_graph(); //A destructor. Depending on how you do things, this may
  //not be necessary.
  
  bool contains(const vertex&) const; //Returns true if the given vertex is in the graph, false otherwise.

  bool adjacent(const vertex&, const vertex&) const; //Returns true if the first vertex is adjacent to the second, false otherwise.

  void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
  void add_edge(const vertex&, const vertex&); //Adds an edge from the first vertex to the second.

  void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges.
  void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.

  std::size_t in_degree(const vertex&) const; //Returns number of edges coming in to a vertex.
  std::size_t out_degree(const vertex&) const; //Returns the number of edges leaving a vertex.
  std::size_t degree(const vertex&) const; //Returns the degree of the vertex (both in and out edges).
  
  std::size_t num_vertices() const; //Returns the total number of vertices in the graph.
  std::size_t num_edges() const; //Returns the total number of edges in the graph.

  std::vector<vertex> get_vertices(); //Returns a vector containing all the vertices.
  std::vector<vertex> get_neighbours(const vertex&) ; //Returns a vector containing the neighbours of the given vertex.

  vertex_iterator<vertex> begin(); //Returns a graph_iterator pointing to the start of the vertex set.
  vertex_iterator<vertex> end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

  neighbour_iterator<vertex> nbegin(const vertex&); //Returns a neighbour_iterator pointing to the start of the neighbour set for the given vertex.
  neighbour_iterator<vertex> nend(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end of the neighbour set for the given vertex.

  std::vector<vertex> depth_first(const vertex&) ; //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
  std::vector<vertex> breadth_first(const vertex&) ; //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

  directed_graph<vertex> out_tree(const vertex&) ; //Returns a spanning tree of the graph starting at the given vertex using the out-edges.
  directed_graph<vertex> in_tree(const vertex&) ; //Returns a spanning tree of the graph starting at the given vertex using the in-edges.

  bool reachable(const vertex&, const vertex&) const; //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.
};

//The vertex_iterator class provides an iterator
//over the vertices of the graph.
//This is one of the harder parts, so if you're
//not too comfortable with C++ leave this for last.
//If you are, there are many ways of doing this,
//as long as it passes the tests, it's okay.
//You may want to watch the videos on iterators before starting.
template <typename vertex>
class vertex_iterator : public std::iterator<std::forward_iterator_tag, directed_graph<vertex>>{

private:
    //You may need data members here.
    directed_graph<vertex> d_graph;
    std::size_t pos;

public:
    vertex_iterator(const vertex_iterator<vertex>& other) : d_graph(other.d_graph), pos(other.pos) {}
    vertex_iterator(const directed_graph<vertex>& graph, std::size_t position = 0) : d_graph(graph), pos(position) {}
    ~vertex_iterator() {}
    vertex_iterator<vertex> operator=(const vertex_iterator<vertex>& other) { 
		d_graph = other.d_graph;
		pos = other.pos;
	}
  	bool operator==(const vertex_iterator<vertex>& other) const { return (pos == other.pos); }
    bool operator!=(const vertex_iterator<vertex>& other) const { return (pos != other.pos); }
    vertex_iterator<vertex> operator++() { 
		++pos;
		return *this; 
	}
    vertex_iterator<vertex> operator++(int) { 
		auto temp = this;
		++pos;
		return temp; 
	}
	//getting the graphs vertices and iterating over them
	//using count to then pick the right vertex in the chain
    const vertex operator*() { 
		int count = 0;
		for(auto& vert : d_graph.get_vertices()){
			if(pos == count){
				return vert;
			}
			count++;
		}
		return vertex();
	}
	//getting the graphs vertices and iterating over them
	//using count to then pick the right vertex in the chain
    const vertex* operator->() { 
		int count = 0;
		for(auto& vert : d_graph.get_vertices()){
			if(pos == count){
				return &vert;
			}
			count++;
		}
		return &vertex();
	}
};

//The neighbour_iterator class provides an iterator
//over the neighbours of a given vertex. This is
//probably harder (conceptually) than the graph_iterator.
//Unless you know how iterators work.
template <typename vertex>
class neighbour_iterator : public std::iterator<std::forward_iterator_tag, directed_graph<vertex>>{

private:

  	//You may need data members here.
  	directed_graph<vertex> d_graph;
	vertex start;
    std::size_t pos;

public:
	neighbour_iterator(const neighbour_iterator<vertex>& other) : d_graph(other.d_graph), start(other.start), pos(other.pos) {}
	neighbour_iterator(const directed_graph<vertex>& graph, const vertex& u, std::size_t position) : d_graph(graph), start(u), pos(position) {}
	~neighbour_iterator() {}
	neighbour_iterator<vertex> operator=(const neighbour_iterator<vertex>& other) { 
		d_graph = other.graph;
		start = other.start;
		pos = other.pos;
	}
 	bool operator==(const neighbour_iterator<vertex>& other) const { return (pos == other.pos); }
 	bool operator!=(const neighbour_iterator<vertex>& other) const { return (pos != other.pos); }
	neighbour_iterator<vertex> operator++() { 
		++pos;
		return *this; 
	}
	neighbour_iterator<vertex> operator++(int) { 
		auto temp = this;
		++pos;
		return temp; 
	}		
	vertex operator*() { 
		//getting the starts neighbouring vertices and iterating over them
		//using count to then pick the right vertex in the chain
		int count = 0;
		for(auto& vert : d_graph.get_neighbours(start)){
			if(pos == count){
				return vert;
			}
			count++;
		}
		return vertex();
	}
	vertex* operator->() { 
		//getting the starts neighbouring vertices and iterating over them
		//using count to then pick the right vertex in the chain
		int count = 0;
		for(auto& vert : d_graph.get_neighbours(start)){
			if(pos == count){
				return &vert;
			}
			count++;
		}
		return &vertex();
	}
};


//Define all your methods down here (or move them up into the header, but be careful you don't double up). If you want to move this into another file, you can, but you should #include the file here.
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.

template <typename vertex> directed_graph<vertex>::directed_graph() {}
template <typename vertex> directed_graph<vertex>::~directed_graph() {}
template <typename vertex> bool directed_graph<vertex>::contains(const vertex& u) const { 
	//iterating over the map to get a vertex and then comparing it to 'u' to see if it exits
	for(auto& map : map_of_vert) {
		if(map.first == u){
			return true;
		}
	}
	return false; 
}
template <typename vertex> bool directed_graph<vertex>::adjacent(const vertex& u, const vertex& v) const { 
	//iterating over the map to get a vertex and then 
	//checking to see if the vertex is 'u'
	//and if it contains an edge to 'v'
	for(auto& map : map_of_vert) {
		if(map.second.count(v) > 0 && map.first == u){
			return true;
		}
	}
	return false; 
}
template <typename vertex> void directed_graph<vertex>::add_vertex(const vertex& u) {
	map_of_vert[u] = std::set<vertex>();
}
template <typename vertex> void directed_graph<vertex>::add_edge(const vertex& u, const vertex& v) {
	map_of_vert[u].insert(v);
}
template <typename vertex> void directed_graph<vertex>::remove_vertex(const vertex& u) {
	map_of_vert.erase(u);
	//removing all possible edges of u from the map
	for(auto& map : map_of_vert) {
		remove_edge(map.first,u);
	}
}
template <typename vertex> void directed_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
	map_of_vert[u].erase(v);
}
template <typename vertex> std::size_t directed_graph<vertex>::in_degree(const vertex& u) const {
	int count = 0;
	//checks to see how many instances of u
	//there are in each vertex
	for(auto& map : map_of_vert) {
		count += map.second.count(u);
	}
	return count; 
}
template <typename vertex> std::size_t directed_graph<vertex>::out_degree(const vertex& u) const { 
	return map_of_vert.at(u).size(); 
}
template <typename vertex> std::size_t directed_graph<vertex>::degree(const vertex& u) const { 
	return in_degree(u) + out_degree(u); 
}
template <typename vertex> std::size_t directed_graph<vertex>::num_vertices() const { 
	return map_of_vert.size(); 
}
template <typename vertex> std::size_t directed_graph<vertex>::num_edges() const { 
	int count = 0;
	//iterating over each vertex and counting each of it's edges
	for(auto& map : map_of_vert) {
		count += map.second.size();
	}
	return count; 
}
template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_vertices() { 
	std::vector<vertex> vertices;
	//iterating over each vertex in the map 
	//and adding it to a vector
	for(auto& map : map_of_vert) {
		vertices.push_back(map.first);
	}
	return vertices; 
}
template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_neighbours(const vertex& u) { 
	std::vector<vertex> neighbours;
	//iterating over the set at key 'u'
	//then adding each vertex to a vector
	for(auto& set : map_of_vert.at(u)) {
		neighbours.push_back(set);
	}
	return neighbours; 
}
template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::begin() { 
	return vertex_iterator<vertex>(*this, 0); 
}
template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::end() { 
	return vertex_iterator<vertex>(*this, num_vertices()); 
}
template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nbegin(const vertex& u) {
	return neighbour_iterator<vertex>(*this, u, 0);
}
template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nend(const vertex& u) {
	return neighbour_iterator<vertex>(*this, u, get_neighbours(u).size());
}
template <typename vertex> std::vector<vertex> directed_graph<vertex>::depth_first(const vertex& u) { 
	std::map<vertex, bool> visited; //used to store visited verticies
	std::vector<vertex> depth; //ordered vector for dfs
	std::stack<vertex> stack;
	stack.push(u);
	//repeats indefinitely until stack is empty
	while(!stack.empty()){
		vertex ts = stack.top(); //top of stack
		stack.pop();
		//double checking that the top of the stack has not been visited
		if(!visited[ts]){
			visited[ts] = true;
			depth.push_back(ts);
			//using a reverse_iterator to reverse the order of get_neighbours
			typename std::vector<vertex>::reverse_iterator rit;
			std::vector<vertex> neighbours = get_neighbours(ts);
			for(rit = neighbours.rbegin(); rit != neighbours.rend(); ++rit){
				if(!visited[*rit]){
					stack.push(*rit);
				}
			}
		}
	}
	return depth;
}
template <typename vertex> std::vector<vertex> directed_graph<vertex>::breadth_first(const vertex& u) { 
	std::map<vertex, bool> visited; //used to store visited verticies
	std::vector<vertex> breadth; //ordered vector for bfs
	std::queue<vertex> queue;
	queue.push(u);
	//repeats indefinitely until queue is empty
	while(!queue.empty()){
		vertex fq = queue.front(); //front of queue
		queue.pop();
		//double checking that the top of the stack has not been visited
		if(!visited[fq]){
			visited[fq] = true;
			breadth.push_back(fq);
			for(auto& edge : get_neighbours(fq)) {
				if(!visited[edge]){
					queue.push(edge);
				}
			}
		}
	}
	return breadth;
}
template <typename vertex> directed_graph<vertex> directed_graph<vertex>::out_tree(const vertex& u) { 
	directed_graph<vertex> tree;
	std::map<vertex, bool> visited; //used to store visited verticies
	std::queue<std::pair<vertex,vertex>> queue;
	queue.push({u, u});
	//repeats indefinitely until queue is empty
	while (!queue.empty()){
		std::pair<vertex,vertex> fq = queue.front(); //front of queue
		queue.pop();
		tree.add_vertex(fq.second);
		if (!visited[fq.second]){
			visited[fq.second] = true;
			//making sure that we aren't adding and edge to itself
			if(fq.first != fq.second){
				tree.add_edge(fq.first, fq.second);
			}
			//getting edges of current vertex
			std::vector<vertex> neighbours = get_neighbours(fq.second);
			for(auto& edge : neighbours){
				if(!visited[edge]){
					queue.push({fq.second, edge});
				}
			}
		}
	}
	return tree;
}
template <typename vertex> directed_graph<vertex> directed_graph<vertex>::in_tree(const vertex& u) { 
	directed_graph<vertex> tree;
	std::map<vertex, bool> visited; //used to store visited verticies
	std::queue<std::pair<vertex,vertex>> queue;
	queue.push({u, u});
	//repeats indefinitely until queue is empty
	while(!queue.empty()){
		std::pair<vertex,vertex> fq = queue.front(); //front of queue
		queue.pop();
		tree.add_vertex(fq.first);
		//making sure that we aren't adding and edge to itself
		if(fq.first != fq.second){
			tree.add_edge(fq.first,fq.second);
		}
		if(!visited[fq.first]){
			visited[fq.first] = true;
			//checking for incoming edges to current vertex
			for(auto& map : map_of_vert){
				if(adjacent(map.first, fq.first) && !visited[map.first]){
					queue.push({map.first,fq.first});
				}
			}
		}
	}
	return tree;
}
template <typename vertex> bool directed_graph<vertex>::reachable(const vertex& u, const vertex& v) const { 
	std::map<vertex, bool> visited; //used to store visited verticies
	std::vector<vertex> depth;
	std::queue<vertex> queue;
	//start from vertex 'u'
	queue.push(u);
	//repeats indefinitely until queue is empty
	while(!queue.empty()){
		vertex fq = queue.front(); //front of queue
		queue.pop();
		if(!visited[fq]){
			visited[fq] = true;
			depth.push_back(fq);
			//if the front of the queue is equal to 'v'
			//then we have successfully traversed the graph
			//from 'u' to 'v'
			if(fq == v){
				return true;
			}
			std::vector<vertex> neighbours;
			for(auto& set : map_of_vert.at(fq)) {
				neighbours.push_back(set);
			}
			//traversing neighbours and pushing non visited edges
			for(auto& edge : neighbours) {
				if(!visited[edge]){
					queue.push(edge);
				}
			}
		}
	}
	return false;
}

#endif