#include <algorithm>
#include <iostream>
#include <vector>
#include <queue>
#include <cmath>

using namespace std;

typedef vector<int> TOUR;
typedef vector<vector<int>> MATRIX;

class Graph {
private:
	typedef struct Edge {
		int u;
		int v;
		int w;
		Edge() {};
		Edge(int _u, int _v, int _w) {
			u = _u; v = _v; w = _w;
		};
	} Edge;

	int e, v;
	vector<Edge> edge_list;
	MATRIX dist;

	/*
	 * For DFS 
	 */
	TOUR visited;
	TOUR result;

	class UnionFind {
	public:
		vector<int> p;
		vector<int> r;
		
		UnionFind(int n) {
			p = vector<int>(n);
			r = vector<int>(n);
			for (int i = 0; i < n; i++) {
				p[i] = -1;
				r[i] = -1;
			}
		}
		
		void makeSet(int v) {
			p[v] = v;
			r[v] = 0;
		}
		
		int _find(int v) {
			if (p[v] != v) {
				p[v] = _find(p[v]);
			}
			return p[v];
		}
		
		void _union(int u, int v) {
			int ur = _find(u);
			int vr = _find(v);
			// Union only if in different sets
			if (ur != vr) {
				if (ur < vr) {
					p[ur] = vr;
				} else if (ur > vr) {
					p[vr] = ur;
				} else {
					p[ur] = vr;
					r[ur] += 1;
				}
			}
		}
	};

	bool compareEdge(Edge e1, Edge e2) {
		if (e1.w != e2.w) {
			return e1.w < e2.w;
		} else {
			if (e1.u != e2.u) {
				return e1.u < e2.u;
			} else {
				return e1.v < e2.v;
			}
		}
	}

	struct comp {
		bool operator() (Edge e1, Edge e2) {
			if (e1.w != e2.w) {
				return e1.w < e2.w;
			} else {
				if (e1.u != e2.u) {
					return e1.u < e2.u;
				} else {
					return e1.v < e2.v;
				}
			}
		}
	} edge_comp;

	vector<Edge> kruskals() {
		int n = edge_list.size();
		vector<Edge> result(0);
		// Create the forest of trees
		UnionFind f = UnionFind(n);
		for (int i = 0; i < n; i++) {
			f.makeSet(i);
		}
		// Sort the graph in order of ascending edge weights
		sort(edge_list.begin(), 
			 edge_list.end(), 
			 edge_comp);

		int added_edges = 0;
		for (int i = 0; i < edge_list.size(); i++) {
			Edge e = edge_list[i];
			int u = e.u;
			int v = e.v;
			if (f._find(u) != f._find(v)) {
				// Append to MST here
				result.push_back(e);
				added_edges++;
				f._union(u, v);
			}
			if (added_edges == n-1) {
				break;
			}
		}
		return edge_list;
	}

	MATRIX to_adj_list(vector<Edge> edge_list) {
		MATRIX adj_list(v, TOUR(0));
		for (int i = 0; i < edge_list.size(); i++) {
			int u = edge_list[i].u;
			int v = edge_list[i].v;
			adj_list[u].push_back(v);
			adj_list[v].push_back(u);
		}
		return adj_list;
	}

	void clear_DFS() {
		visited = TOUR(v, 0);
		result = TOUR(0);
	}

	void DFS(MATRIX adj_list, int start) {
		visited[start] = 1;
		result.push_back(start);
		for (int i = 0; i < adj_list[start].size(); i++) {
			int next = adj_list[start][i];
			if (visited[next] == 0) {
				DFS(adj_list, next);
			}
		}
	}

	TOUR get_tour() {
		TOUR new_tour(0);
		vector<bool> mask(v, false);
		for (int i = 0; i < result.size(); i++) {
			if (mask[result[i]] == false) {
				new_tour.push_back(result[i]);
				mask[result[i]] = true;
			}
		}
		return new_tour;
	}

	int tour_length(TOUR tour) {
		int length = 0;
		for (int i = 0; i < tour.size()-1; i++) {
			length += dist[tour[i]][tour[i+1]];
		}
		length += dist[tour[tour.size()-1]][tour[0]];
		return length;
	}

public:
	// Initialize a Graph with a distance matrix
	Graph(MATRIX _dist) {
		dist = _dist;
		v = dist.size();
		// Initialize the edge list
		edge_list = vector<Edge>(0);
		for (int i = 0; i < v; i++) {
			for (int j = i; j < v; j++) {
				if (i != j) {
					edge_list.push_back(Edge(i, j, dist[i][j]));
				}
			}
		}
	}

	TOUR two_approx_tour() {
		// Run MST on the graph and get the edges
		vector<Edge> mst = kruskals();
		// Run DFS on the MST
		clear_DFS();
		DFS(to_adj_list(mst), 0);
		// Remove duplicate vertices
		return get_tour();
	}

	TOUR multi_frag_tour() {
		MATRIX adj_list(v, TOUR(0));
		// Sort the edge list
		sort(edge_list.begin(), edge_list.end(), edge_comp);
		// Use a bitmask to ensure outdegree <= 2
		vector<int> mask(v, 0);

		UnionFind f = UnionFind(v);
		for (int i = 0; i < v; i++) {
			f.makeSet(i);
		}

		// Add the cheapest edges while maintaining
		// 1. No internal cycles
		// 2. No vertices with outdegree > 2
		int sum = 0;
		for (int i = 0; i < edge_list.size(); i++) {
			Edge e = edge_list[i];
			int u = e.u, v = e.v;
			int uf = f._find(u);
			int vf = f._find(v);
			//printf("uf = %d, vf = %d\n", uf, vf);
			if (uf != vf) {
				if (mask[u] <= 1 && mask[v] <= 1) {
					adj_list[u].push_back(v);
					adj_list[v].push_back(u);
					f._union(u, v);
					mask[u]++; mask[v]++;
					sum += dist[u][v];
					//printf("adding (%d,%d)\n", u, v);					
				}
			}
		}

		// Finding the first tail end value and break
		int tail;
		for (int i = 0; i < mask.size(); i++) {
			if (mask[i] == 1) {
				tail = i;
				break;
			}
		}

		// Iterative DFS starting from a tail end of the tour
		TOUR mf_tour(0);
		vector<bool> taken(v, false);
		int n_taken = 1, start = tail;
		mf_tour.push_back(start);
		taken[start] = true;
		while (n_taken < v) {
			int next = 0;
			if (adj_list[start].size() == 1) {
				next = adj_list[start][0];
				mf_tour.push_back(next);
				taken[adj_list[start][0]] = true;
			} else {
				for (int j = 0; j < adj_list[start].size(); j++) {
					if (taken[adj_list[start][j]] == false) {
						next = adj_list[start][j];
						mf_tour.push_back(next);
						taken[adj_list[start][j]] = true;
						break;
					}
				}			
			}
			n_taken++;
			start = next;
		}

		return mf_tour;
	}
};

