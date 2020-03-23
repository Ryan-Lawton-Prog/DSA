/*
 * Notice that the list of included headers has
 * expanded a little. As before, you are not allowed
 * to add to this.
 */
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
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
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>

#include "directed_graph.hpp"

/*
 * Computes whether the input is a Directed Acyclic Graph (DAG).
 * A digraph is a DAG if there is no vertex that has a cycle.
 * A cycle is a non-empty set of [out-]edges that starts at one 
 * vertex, and returns to it.
 */
template <typename vertex>
bool is_dag(const directed_graph<vertex> & d) {
  typename std::unordered_set<vertex>::const_iterator jt;
  for(auto it : d){
    std::unordered_map<vertex, bool> visited;
    std::queue<vertex> queue;
    vertex start = it;
    queue.push(start);
    visited[start] = true;
    //repeats indefinitely until queue is empty
    while(!queue.empty()){
      vertex next = queue.front();
      queue.pop();
      //iterates over the front of queues neighbours
      for(jt = d.nbegin(next); jt != d.nend(next); jt++){
        if(*jt == start){
          return false;
        }else{
          //if neighbour has not been visited, add to queue
          if(!visited[*jt]){
            queue.push(*jt);
            visited[*jt] = true;
          }
        }
      } 
    }
  }
  return true;
}

/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 */
template <typename vertex>
std::list<vertex> topological_sort(const directed_graph<vertex> & d) {
  //Kahn's Algorithm
  typename std::unordered_set<vertex>::const_iterator jt;
  directed_graph<vertex> graph = d;
  std::list<vertex> sorted_graph;
  std::queue<vertex> no_edge;
  
  //Push all vertices with no incoming edges to queue
  for(auto it : d){
    if(d.in_degree(it) == 0){
      no_edge.push(it);
    }
  }
  
  //repeats indefinitely until queue is empty
  while(!no_edge.empty()){
    vertex next = no_edge.front();
    no_edge.pop();
    sorted_graph.push_back(next);
    //Iterates over the front of queues neighbours
    for(jt = d.nbegin(next); jt != d.nend(next); jt++){
      graph.remove_edge(next, *jt);
      //If neighbour has no incoming edges now, add the queue
      if(graph.in_degree(*jt) == 0){
        no_edge.push(*jt);
      }
    }
  }
  return sorted_graph;
}

/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {
  //Sort graph using topological sort
  std::list<vertex> topo_graph = topological_sort(d);
  std::queue<vertex> queue;
  //If graph has only one or no vertex, graph is a Hamiltonian DAG
  //otherwise add top of sorted graph to queue
  if(topo_graph.size() > 1){
    queue.push(topo_graph.front());
    topo_graph.pop_front();
  }else{
    return true;
  }
  
  //repeats indefinitely until queue is empty
  while(!queue.empty()){
    vertex next = queue.front();
    queue.pop();
    //if next vertex in sorted graph is a neighbour to current vertex
    //add to queue, or if no remaining vertex in sorted graph, graph 
    //is a Hamiltonian DAG
    if(d.adjacent(next, topo_graph.front())){
      queue.push(topo_graph.front());
      topo_graph.pop_front();
    }else if(topo_graph.size() == 0){
      return true;
    }
  }
  
  return false;
}

/*
 * Computes the weakly connected components of the graph.
 * A [weak] component is the smallest subset of the vertices
 * such that the in and out neighbourhood of each vertex in
 * the set is also contained in the set.
 */
template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
  directed_graph<vertex> graph = d;
  std::vector<std::vector<vertex>> components;
  typename std::unordered_set<vertex>::const_iterator jt;
  std::queue<vertex> queue;
  std::unordered_map<vertex, bool> visited;
  
  //Converting graph from directed to normal
  for(auto it : d){
    for(jt = d.nbegin(it); jt != d.nend(it); jt++){
      graph.add_edge(*jt, it);
    }
    //initializing visited
    visited[it] = false;
  }
  
  
  for(auto it : graph){
    std::vector<vertex> vertices;
    if(!visited[it]){
      queue.push(it);
      visited[it] = true;
      vertices.push_back(it);
      //DFS to find loosley connected components
      while(!queue.empty()){
        vertex next = queue.front();
        queue.pop();
        for(jt = graph.nbegin(next); jt != graph.nend(next); jt++){
          if(!visited[*jt]){
            vertices.push_back(*jt);
            queue.push(*jt);
            visited[*jt] = true;
          }
        }
      }
      //adding component
      components.push_back(vertices);
    }
  }
  
  return components;
}

/*
 * Tarjan's Algorithm to find Strongly Connect Components
 */
template <typename vertex>
void strongconnect(const directed_graph<vertex> & d, vertex v, std::unordered_map<vertex, int>*index, std::unordered_map<vertex, int>*low, 
    std::stack<vertex>*stack, std::unordered_map<vertex, bool>*on_stack, std::vector<std::vector<vertex>>*components){
  typename std::unordered_set<vertex>::const_iterator it;
  static int i = 0;
  
  //inilializing index
  (*index)[v] = (*low)[v] = ++i;
  stack->push(v);
  (*on_stack)[v] = true;
  
  //iterate over neighbours of vertex v
  for(it = d.nbegin(v); it != d.nend(v); it++){
    //If we can't find neighbour in index add it
    if(index->find(*it) == index->end()){
      //recursively find next vertex
      strongconnect(d, *it, index, low, stack, on_stack, components);
      (*low)[v] = std::min((*low)[v], (*low)[*it]);
    } else if ((*on_stack)[*it]){ //looped back
      (*low)[v] = std::min((*low)[v], (*index)[*it]);
    }
  }
  
  //if index is equal to low we found a component
  if((*low)[v] == (*index)[v]){
    std::vector<vertex> component;
    vertex next;
    //adding each component on stack to vector
    do {
      next = stack->top(); stack->pop();
      (*on_stack)[next] = false;
      component.push_back(next);
    } while(next != v);
    //adding vector components to a vector of vectors
    components->push_back(component);
  }
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */
template <typename vertex>
std::vector<std::vector<vertex>> strongly_connected_components(const directed_graph<vertex> & d) {
  std::vector<std::vector<vertex>>*components = new std::vector<std::vector<vertex>>(); //final list of components
  std::unordered_map<vertex, int>*index = new std::unordered_map<vertex, int>(); //tracks vertex's index in traversal
  std::unordered_map<vertex, int>*low = new std::unordered_map<vertex, int>(); //tracks lowest index
  std::unordered_map<vertex, bool>*on_stack = new std::unordered_map<vertex, bool>(); //tracks if a vertex is on the stack
  std::stack<vertex>*stack = new std::stack<vertex>(); //stack of vertices
  
  //loop over each vertex in graph
  for(auto it : d){
    if(index->find(it) == index->end()){
      strongconnect(d, it, index, low, stack, on_stack, components);
    }
  }
  return *components;
}

/*
 * Computes the shortest distance from u to every other vertex
 * in the graph d. The shortest distance is the smallest number
 * of edges in any path from u to the other vertex.
 * If there is no path from u to a vertex, set the distance to
 * be the number of vertices in d plus 1.
 */
template <typename vertex>
std::unordered_map<vertex, std::size_t> shortest_distances(const directed_graph<vertex> & d, const vertex & u) {
  //Dijkstra's Algorithm
  std::unordered_map<vertex, std::size_t> distances;
  std::queue<vertex> queue;
  std::unordered_map<vertex, std::pair<bool,int>> unvisited; //map of pairs containing if it has been visited and it's depth
  typename std::unordered_set<vertex>::const_iterator it;

  //initializing map
  for(vertex it : d){
    unvisited[it] = {false, 0};
  }
  queue.push(u);
  unvisited[u].first = true;
  //loops indefinetly until queue is empty
  while(!queue.empty()){
    vertex next = queue.front();
    queue.pop();
    distances[next] = unvisited[next].second;
    //iterates over neighbours
    for(it = d.nbegin(next); it != d.nend(next); it++){
      //check to see if visited
      if(!unvisited[*it].first){
        queue.push(*it);
        //if neighbour, it's depth is current vertex depth + 1
        unvisited[*it] = {true, unvisited[next].second + 1};
      }
    }
  }
  
  //setting all none connected vertices to num of vertices in graph + 1
  for(auto it : d){
    if(!unvisited[it].first){
      distances[it] = d.num_vertices()+1;
    }
  }
  
  return distances;
}

