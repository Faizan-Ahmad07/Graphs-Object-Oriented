# Graph-OOPS

A C++ project demonstrating object-oriented programming (OOP) concepts for various types of graphs, including unweighted/weighted and directed/undirected graphs. The project provides implementations and tests for:

- **Graph**: Undirected, unweighted graph
- **DirectedGraph**: Directed, unweighted graph
- **WeightedGraph**: Undirected, weighted graph
- **WeightedDirectedGraph**: Directed, weighted graph

## Features
- Add and remove edges
- Breadth-First Search (BFS)
- Depth-First Search (DFS)
- Cycle detection
- Connected components count
- Shortest path (Dijkstra for weighted, BFS for unweighted)
- Minimum Spanning Tree (Prim's algorithm)
- Topological sort (for directed graphs)

## File Structure
- `graph_oops.cpp` — All class definitions, implementations, and a comprehensive set of tests in `main()
## How to Build and Run
1. **Compile:**
   ```sh
   g++ -std=c++17 -o graph_oops.exe graph_oops.cpp
   ```
2. **Run:**
   ```sh
   ./graph_oops.exe
   ```
   or on Windows:
   ```sh
   graph_oops.exe
   ```

## Example Output
The program runs a suite of tests for all graph types, showing:
- Adjacency lists
- BFS/DFS traversals
- Cycle checks
- Connected components
- Shortest paths
- Minimum spanning tree edges
- Topological sort (for directed graphs)

## Customization
You can modify or extend the tests in `main()` to experiment with different graph structures and algorithms.

## Requirements
- C++17 or later
- Standard C++ libraries (no external dependencies)