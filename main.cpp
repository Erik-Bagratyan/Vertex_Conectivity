#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

const unsigned int MAXN = 1005; 
const int INFTY = 1000000;

int flow[MAXN][MAXN];	
int cap[MAXN][MAXN];

int pv[MAXN];

int q[MAXN];
int q_front,q_back;

vector<vector<int>> graph;
int new_nodes[MAXN];

void print_matrix(int n, int m)
{
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < m; ++j)
		{
			cout << cap[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

bool bfs(int s, int t,const int& n)
{
	int u = s;
	q_front = q_back = 0;
	for(int i=0;i<=n;i++) pv[i] = -1;
	pv[u] = -2;
	q[q_front++] = u;
	while(q_back < q_front && pv[t] == -1)
	{ 
		u = q[q_back++]; 
		for(int v=1;v<=n;v++)
		{
			if(pv[v] == -1 && flow[u][v] < cap[u][v])
			{ 
				q[q_front++] = v;
				pv[v] = u;
			}
		}
	}
	return pv[t] != -1;
}

int maxflow(int s,int t, int n)
{
	int result = 0,u,v;
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<=n;j++)
		{
			flow[i][j] = 0;
		}
	}


	int min_res_cap;
	while(bfs(s,t,n))
	{
		u = pv[t]; v = t;
		if(u<0) break;
		min_res_cap = INFTY;
		while(u>0)
		{
			min_res_cap = min(min_res_cap, cap[u][v]-flow[u][v]);
			v = u;
			u = pv[u];
		}
		u = pv[t]; v = t;
		while(u>0)
		{
			flow[u][v] += min_res_cap;
			flow[v][u] -= min_res_cap;
			v = u;
			u = pv[u];
		}
		result += min_res_cap;

	}
	return result;
}

int transform_flow(int s, int t, int n)
{
	int ns,nt,cnt=0;
	for(int i=0;i<=2*(n+1);i++)
	{
		for(int j=0;j<=2*(n+1);j++){
			cap[i][j] = 0;
		}
		new_nodes[i] = 0;
	}
	for(int i=0;i<=n;i++)
	{
		if(i==s) {
			ns = cnt;
			new_nodes[2*i] = cnt;
			new_nodes[2*i+1] = cnt++;
		} else if(i==t){
			nt = cnt;
			new_nodes[2*i] = cnt;
			new_nodes[2*i+1] = cnt++;
		} else {
			new_nodes[2*i] = cnt++;
			new_nodes[2*i+1] = cnt++;
			cap[new_nodes[2*i]][new_nodes[2*i+1]] = 1;
		}
	}
	
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<graph[i].size();j++)
		{
			cap[new_nodes[graph[i][j]*2+1]][new_nodes[i*2]] = 1;
		}
		for(int j=0;j<graph[i].size();j++)
		{
			cap[new_nodes[i*2+1]][new_nodes[graph[i][j]*2]] = 1;
		}
	}
	print_matrix(2*n+1, 2*n+1);
	return maxflow(ns,nt,cnt);
}

int max_connectivity(int n)
{
	int minc = INFTY;
	for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)
		{
			if(i!=j)
			{
				minc = min(minc, transform_flow(i,j,n));
			}
		}
	}
	return minc;
}

void input_graph(vector<vector<int>>& graph, int& n, int& m)
{
	graph.resize(MAXN);	
    std::ifstream input("input.txt");
    if (!input.is_open())
    {
        std::cerr << "File is not open !!! " << std::endl;
        exit(0);
    }
	input >> n; // Vertex
	input >> m; // Egde
	input.get();
	int i = 1;
    	while(true)
	{
		std::vector<int> elem;
		std::string temp;
		while(input.peek() != '\n' && input >> temp)
		{
			elem.push_back(std::stoi(temp));
		}
		if(!input.eof())
		{
			graph[elem[0]].push_back(elem[1]);
			graph[elem[1]].push_back(elem[1]);
			std::cout << i++ << ") " << elem[0] << " - " << elem[1] << std::endl;
			input.get();
		}
		else break;
	}
}

int main(){

	int n,m,u,v,res;
	n = m = 0;
	input_graph(graph, n, m);
	res = max_connectivity(n);
	std::cout << "Min vertex : " << res << std::endl;
	return 0;
}
