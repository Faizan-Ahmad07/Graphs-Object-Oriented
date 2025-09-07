
#include <bits/stdc++.h>
using namespace std;

void printVector(const vector<int> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " ";
    }
    cout << endl;
}
// This class is the base class for other graph-related classes such as DirectedGraph, WeightedGraph, and WeightedDirectedGraph
class Graph
{
private:
    // making this private so only friend class i.e. DirectedGraph can access it
    vector<vector<int>> adjList;
    void numConnectedComponentsHelper(int node, vector<int> &visited);
    bool isCyclicHelper(int node, vector<int> &visited, int parent);

protected:
    //  making this protected as numNodes will be common for all graph types
    int numNodes;

public:
    Graph(int nodes); // constructor
    virtual ~Graph() {}
    virtual void addEdge(int v, int w);
    virtual void removeEdge(int v, int w);
    // this display function will be used both in Graph and DirectedGraph classes
    virtual void display();
    virtual vector<int> bfs(int node);
    virtual vector<int> dfs(int node);

    // FOR YOU----> // used call by reference, for just using the address, instead of unnecessary copying
    virtual void dfsHelper(int node, vector<int> &visited, vector<int> &ans);

    bool isCyclic();
    int numConnectedComponents();
    virtual int ShortestPath(int start, int destination);
    // this way only DerictedGraph class can access adjList
    friend class DirectedGraph;
    friend ostream & operator << (ostream &out, const Graph &graph);
};

ostream & operator << (ostream &out, const Graph &graph)
{
    for (int node = 0; node < graph.numNodes; node++){
        out << node << " -> ";
        for(auto neighbour : graph.adjList[node]){
            out << neighbour << " ";
        }
        out << endl;
    }
    return out;
}
// Graph::Graph(int nodes): This is the constructor for the Graph class, which takes an integer argument vertices to specify the number of nodes in the graph.
// : numNodes(nodes): This is a member initialization list, which initializes the numNodes attribute with the value of the vertices argument.
Graph::Graph(int nodes) : numNodes(nodes)
{
    // Resizes the adjacency list to accommodate the specified number of vertices
    adjList.resize(nodes);
}
void Graph::addEdge(int v, int w)
{
    adjList[v].push_back(w);
    adjList[w].push_back(v);
}
void Graph::removeEdge(int v, int w)
{
    adjList[v].erase(find(adjList[v].begin(), adjList[v].end(), w));
    adjList[w].erase(find(adjList[w].begin(), adjList[w].end(), v));
}
void Graph::display()
{
    for (int i = 0; i < numNodes; i++)
    {
        cout << i << " -> ";
        for (int j = 0; j < adjList[i].size(); j++)
        {
            cout << adjList[i][j] << " ";
        }
        cout << endl;
    }
}
vector<int> Graph::bfs(int vertex)
{
    vector<int> ans;
    queue<int> q;
    q.push(vertex);
    vector<int> visited(numNodes, 0);
    visited[vertex] = 1;
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        ans.push_back(node);
        for (auto neighbour : adjList[node])
        {
            if (!visited[neighbour])
            {
                visited[neighbour] = 1;
                q.push(neighbour);
            }
        }
    }
    return ans;
}
void Graph::dfsHelper(int node, vector<int> &visited, vector<int> &ans)
{
    visited[node] = 1;
    ans.push_back(node);
    for (auto neighbour : adjList[node])
    {
        if (!visited[neighbour])
        {
            dfsHelper(neighbour, visited, ans);
        }
    }
}
vector<int> Graph::dfs(int node)
{
    vector<int> ans;
    vector<int> visited(numNodes, 0);
    visited[node] = 1;
    ans.push_back(node);
    for (auto neighbour : adjList[node])
    {
        dfsHelper(neighbour, visited, ans);
    }
    return ans;
}
bool Graph::isCyclicHelper(int node, vector<int>& visited, int parent)
{
    visited[node] = 1;
    for(auto neighbour : adjList[node]){
        if(!visited[neighbour]){
            if(isCyclicHelper(neighbour, visited, node)){
                return true;
            }
        }else if(neighbour != parent){
            return true;
        }
    }
    return false;
}
bool Graph::isCyclic()
{
    vector<int> visited(numNodes, 0);
    for(int node = 0; node < numNodes; node++){
        if(!visited[node])
        {
            if(isCyclicHelper(node, visited, -1)==true)
            {
                return true;
            }
        }
    }
    return false;
}
void Graph::numConnectedComponentsHelper(int node, vector<int> &visited)
{
    visited[node] = 1;
    for (auto neighbour : adjList[node])
    {
        if (!visited[neighbour])
        {
            numConnectedComponentsHelper(neighbour, visited);
        }
    }
}
int Graph::numConnectedComponents()
{
    vector<int> visited(numNodes, 0);
    int result = 0;
    for (int node = 0; node < numNodes; node++)
    {
        if (!visited[node])
        {
            numConnectedComponentsHelper(node, visited);
            result++;
        }
    }
    return result;
}
int Graph::ShortestPath(int start, int destination)
{
    // FOR YOU---> // this array is used to store the minimum distance from start to any node encountered during the loop
    vector<int> dp(numNodes, 1e9);
    dp[start] = 0;
    priority_queue <pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push(make_pair(0, start));
    while(!q.empty())
    {
        // lambda capture 
        pair<int, int> top = q.top();
        int dist = top.first;
        int node = top.second;
        q.pop();
        if (node == destination) break;
        for (auto neighbour : adjList[node])
        {
            if (dp[neighbour] > dist + 1)
            {
                dp[neighbour] = dist + 1;
                q.push(make_pair(dp[neighbour], neighbour));
            }
        }
    }
    // -1 for not reachable from start to destination
    if (dp[destination] == 1e9) return -1;
    return dp[destination];
}
class DirectedGraph : public Graph
{
private:
    bool isCyclicHelper(int node, vector<int>& visited);
public:
    DirectedGraph(int nodes);
    void addEdge(int v, int w);
    void removeEdge(int v, int w);
    // display function of Graph class will be used
    // topological sort using Kahn algorithm
    virtual vector<int> topoSort();
    bool isCyclic();
};
DirectedGraph::DirectedGraph(int nodes) : Graph(nodes) {}
void DirectedGraph::addEdge(int v, int w)
{
    adjList[v].push_back(w); // v:source  w:destination
}
void DirectedGraph::removeEdge(int v, int w)
{
    adjList[v].erase(find(adjList[v].begin(), adjList[v].end(), w));
}

vector<int> DirectedGraph::topoSort()
{
    vector<int> ans;
    vector<int> indegree(numNodes, 0);
    queue<int> q;
    for (int node = 0; node < numNodes; node++)
    {
        for (auto neighbour : adjList[node])
        {
            indegree[neighbour]++;
        }
    }
    for (int node = 0; node < numNodes; node++)
    {
        if (indegree[node] == 0)
        {
            q.push(node);
        }
    }
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        ans.push_back(node);
        for (auto neighbour : adjList[node])
        {
            indegree[neighbour]--;
            if (indegree[neighbour] == 0)
            {
                q.push(neighbour);
            }
        }
    }
    return ans;
}
bool DirectedGraph::isCyclicHelper(int node, vector<int>& visited)
{
    visited[node] = 1;
    for (auto neighbour : adjList[node]){
        if(!visited[neighbour]){
            if(isCyclicHelper(neighbour, visited)){
                return true;
            }
        }else if (visited[neighbour] == 1){
            return true;
        }
    }
    visited[node] = 2;
    return false;
}
bool DirectedGraph::isCyclic()
{
    vector<int> visited(numNodes, 0);
    for(int node = 0; node < numNodes; node++){
        if (!visited[node])
        {
            if(isCyclicHelper(node, visited)==true)
            {
                return true;
            }
        }
    }
    return false;
}

class WeightedGraph : public Graph
{
protected:
    vector<vector<pair<int, int>>> weightedAdjList;
    void numConnectedComponenetsHelper(int node, vector<int> &visited);
    bool isCyclicHelper(int node, vector<int> & visited, int parent);

public:
    WeightedGraph(int nodes);
    virtual void addEdge(int v, int w, int weight);
    virtual void removeEdge(int v, int w);
    vector<int> bfs(int node);
    vector<int> dfs(int node);
    void dfsHelper(int node, vector<int> &visited, vector<int> &ans);

    int numConnectedComponents();
    virtual int ShortestPath(int start, int target);
    vector<pair<int,int>> MinimumSpanningTree();

    void display();
    friend ostream & operator << (ostream &out, const WeightedGraph &graph);
};

vector<int> WeightedGraph::bfs(int vertex)
{
    vector<int> ans;
    queue<int> q;
    q.push(vertex);
    vector<int> visited(numNodes, 0);
    visited[vertex] = 1;
    while (!q.empty())
    {
        int node = q.front();
        q.pop();
        ans.push_back(node);
        for (auto neighbour : weightedAdjList[node])
        {
            if (!visited[neighbour.first])
            {
                visited[neighbour.first] = 1;
                q.push(neighbour.first);
            }
        }
    }
    return ans;
}

void WeightedGraph::dfsHelper(int node, vector<int> &visited, vector<int> &ans)
{
    visited[node] = 1;
    ans.push_back(node);
    for (auto neighbour : weightedAdjList[node])
    {
        if (!visited[neighbour.first])
        {
            dfsHelper(neighbour.first, visited, ans);
        }
    }
}

vector<int> WeightedGraph::dfs(int node)
{
    vector<int> ans;
    vector<int> visited(numNodes, 0);
    visited[node] = 1;
    ans.push_back(node);
    for (auto neighbour : weightedAdjList[node])
    {
        dfsHelper(neighbour.first, visited, ans);
    }
    return ans;
}

ostream & operator << (ostream &out, const WeightedGraph &graph){
    for (int node = 0; node < graph.numNodes; node++){
        out << node << " -> ";
        for(auto neighbour : graph.weightedAdjList[node]){
            out << neighbour.first << " " << neighbour.second << "  ";
        }
        out << endl;
    }
    return out;
}

WeightedGraph::WeightedGraph(int nodes) : Graph(nodes)
{
    weightedAdjList.resize(numNodes);
}
void WeightedGraph::addEdge(int v, int w, int weight)
{
    weightedAdjList[v].push_back(make_pair(w, weight));
    weightedAdjList[w].push_back(make_pair(v, weight));
}
void WeightedGraph::removeEdge(int v, int w)
{
    for (auto itr = weightedAdjList[v].begin(); itr != weightedAdjList[v].end(); itr++)
    {
        if (itr->first == w)
        {
            weightedAdjList[v].erase(itr);
            break;
        }
    }
    for (auto itr = weightedAdjList[w].begin(); itr != weightedAdjList[w].end(); itr++)
    {
        if (itr->first == v)
        {
            weightedAdjList[w].erase(itr);
            break;
        }
    }
}
void WeightedGraph::numConnectedComponenetsHelper(int node, vector<int> & visited)
{
    visited[node] = 1;
    for (auto neighbour : weightedAdjList[node])
    {
        if (!visited[neighbour.first])
        {
            WeightedGraph::numConnectedComponenetsHelper(neighbour.first, visited);
        }
    }
}
int WeightedGraph::numConnectedComponents()
{
    vector<int> visited(numNodes, 0);
    int result = 0;
    for (int node = 0; node < numNodes; node++)
    {
        if (!visited[node])
        {
            WeightedGraph::numConnectedComponenetsHelper(node, visited);
            result++;
        }
    }
    return result;
}
bool WeightedGraph::isCyclicHelper(int node, vector<int> & visited, int parent)
{
    visited[node] = 1;
    for(auto neighbour : weightedAdjList[node]){
        if(!visited[neighbour.first]){
            if(WeightedGraph::isCyclicHelper(neighbour.first, visited, node)){
                return true;
            }
        }else if(neighbour.first != parent){
            return true;
        }
    }
    return false;
}
int WeightedGraph::ShortestPath(int start, int destination)
{
    vector<int> dp(numNodes, 1e9);
    dp[start] = 0;
    priority_queue <pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push(make_pair(0, start));
    while(!q.empty())
    {
        auto top = q.top();
        int dist = top.first;
        int node = top.second;
        q.pop();
        if (node == destination) break;
        for (auto neighbour : weightedAdjList[node])
        {
            if (dp[neighbour.first] > dist + neighbour.second)
            {
                dp[neighbour.first] = dist + neighbour.second;
                q.push(make_pair(dp[neighbour.first], neighbour.first));
            }
        }
    }
    if (dp[destination] == 1e9) return -1;
    return dp[destination];
}
vector<pair<int,int>> WeightedGraph::MinimumSpanningTree()
{
    vector<int> included(numNodes, 0);
    vector<int> dp(numNodes, 1e9);
    dp[0] = 0;
    vector<int> parent(numNodes);
    parent[0] = 0;
    priority_queue <pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> q;
    q.push(make_pair(0, 0));
    while(!q.empty())
    {
        auto top = q.top();
        int dist = top.first;
        int node = top.second;
        q.pop();
        if (included[node] == 1) continue;
        included[node] = 1;
        for (auto neighbour : weightedAdjList[node])
        {
            if (!included[neighbour.first] && dp[neighbour.first] > neighbour.second) {
                dp[neighbour.first] = neighbour.second;
                parent[neighbour.first] = node;
                q.push(make_pair(dp[neighbour.first], neighbour.first));
            }
        }
    }
    vector<pair<int,int>> ans;
    for (int node = 1; node < numNodes; node++)
    {
        ans.push_back(make_pair(node, parent[node]));
    }
    return ans;
}
void WeightedGraph::display()
{
    for (int i = 0; i < numNodes; i++)
    {
        cout << i << " -> ";
        for (int j = 0; j < weightedAdjList[i].size(); j++)
        {
            cout << weightedAdjList[i][j].first << " " << weightedAdjList[i][j].second << "  ";
        }
        cout << endl;
    }
}

class WeightedDirectedGraph : public WeightedGraph
{
public:
    WeightedDirectedGraph(int vertices);
    void addEdge(int v, int w, int weight);
    void removeEdge(int v, int w);
};
WeightedDirectedGraph::WeightedDirectedGraph(int vertices) : WeightedGraph(vertices) {}
void WeightedDirectedGraph::addEdge(int v, int w, int weight)
{
    weightedAdjList[v].push_back(make_pair(w, weight));
}
void WeightedDirectedGraph::removeEdge(int v, int w)
{
    for (auto itr = weightedAdjList[v].begin(); itr != weightedAdjList[v].end(); itr++)
    {
        if (itr->first == w)
        {
            weightedAdjList[v].erase(itr);
            break;
        }
    }
}

int main()
{
    cout << "==== Test 1: Undirected Unweighted Graph ====" << endl;
    Graph ug(6);
    ug.addEdge(0, 1);
    ug.addEdge(0, 2);
    ug.addEdge(1, 3);
    ug.addEdge(2, 3);
    ug.addEdge(4, 5);
    ug.display();
    cout << "BFS from 0: "; printVector(ug.bfs(0));
    cout << "DFS from 0: "; printVector(ug.dfs(0));
    cout << "Is Cyclic: " << ug.isCyclic() << endl;
    cout << "Number of Connected Components: " << ug.numConnectedComponents() << endl;
    cout << "Shortest Path 0->3: " << ug.ShortestPath(0, 3) << endl;

    cout << "\n==== Test 2: Directed Graph ====" << endl;
    DirectedGraph dg(6);
    dg.addEdge(0, 1);
    dg.addEdge(0, 2);
    dg.addEdge(1, 3);
    dg.addEdge(2, 3);
    dg.addEdge(3, 4);
    dg.addEdge(4, 5);
    dg.display();
    cout << "Topological Sort: "; printVector(dg.topoSort());
    cout << "Is Cyclic: " << dg.isCyclic() << endl;

    cout << "\n==== Test 3: Weighted Undirected Graph ====" << endl;
    WeightedGraph wg(5);
    wg.addEdge(0, 1, 2);
    wg.addEdge(0, 2, 4);
    wg.addEdge(1, 2, 1);
    wg.addEdge(1, 3, 7);
    wg.addEdge(2, 4, 3);
    wg.display();
    cout << "BFS from 0: "; printVector(wg.bfs(0));
    cout << "DFS from 0: "; printVector(wg.dfs(0));
    cout << "Number of Connected Components: " << wg.numConnectedComponents() << endl;
    cout << "Shortest Path 0->4: " << wg.ShortestPath(0, 4) << endl;
    cout << "Minimum Spanning Tree (edges): ";
    for (auto p : wg.MinimumSpanningTree()) cout << "(" << p.first << "," << p.second << ") ";
    cout << endl;

    cout << "\n==== Test 4: Weighted Directed Graph ====" << endl;
    WeightedDirectedGraph wdg(5);
    wdg.addEdge(0, 1, 10);
    wdg.addEdge(0, 2, 3);
    wdg.addEdge(1, 2, 1);
    wdg.addEdge(2, 1, 4);
    wdg.addEdge(2, 3, 2);
    wdg.addEdge(3, 4, 2);
    wdg.addEdge(4, 3, 9);
    wdg.display();
    cout << "BFS from 0: "; printVector(wdg.bfs(0));
    cout << "DFS from 0: "; printVector(wdg.dfs(0));
    cout << "Shortest Path 0->4: " << wdg.ShortestPath(0, 4) << endl;

    cout << "\n==== Test 5: Edge Removal and Display ====" << endl;
    ug.removeEdge(0, 1);
    ug.display();
    wg.removeEdge(0, 1);
    wg.display();
    wdg.removeEdge(0, 1);
    wdg.display();

    return 0;
}
