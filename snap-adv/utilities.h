#ifndef utilities_h
#define utilities_h


#include "datastructures.h"
#include "temporalmotiftypes.h"


#include "Snap.h"

#include <random>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/functional/hash.hpp>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iostream>
#include <cstdlib> 



// TODO: REPLACE THIS WITH SNAP STRUCTURES

using DataStructures::TEdge;


namespace Utils {
	int LoadEdgeList(const std::string& filenameG, std::vector<TEdge>& edges_, std::vector<size_t>& temporal_degrees_);
	int LoadEdgeList(const std::string& filenameG, TVec< THash<TInt, TVec<TInt64> >> &temporal_data, std::vector<size_t>& temporal_degrees_, long long int& totedges, timestamp_t& timespan);
	int LoadTemporalMotif(const std::string& filenameM, std::vector<std::pair<int, int>>& edgesM_, int& Vm_, int& l_);

	int RemapNodes(std::vector<TEdge>& edgeList);
	int PrintGraph(std::vector<TEdge>& edgeList);
	int GetTrianglesFromNode(PUNGraph graph, const int startingNode, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles);
	int GetTrianglesFromEdge(const PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles);
	int GetTrianglesFromEdgeThreads(const PUNGraph& graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles);
	int GetWedgesFromEdge(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& wedges);
	int GetSquaresFromEdge(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& squares);
	int GetSquaresFromEdgeMin(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& squares);
}


// Object to represent a static directed edge src -> dst and all the timestamps that are mapped on such static edges
class StaticEdge{
	private:
		std::pair<node_t, node_t> m_srcdst;

	public:
		StaticEdge(const TEdge edge) {
			m_srcdst.first = edge.src;
			m_srcdst.second = edge.dst;
		}
		node_t getSrc() { return m_srcdst.first; }
		node_t getDst() { return m_srcdst.second; }

		/*
		void addTimestamp(const timestamp_t time) { m_timestamps.insert(time); };
		std::vector<TEdge> getTemporalEdges() 
		{
			std::vector<TEdge> temporal_edges;
			long long ID = 0;
			for(auto it = m_timestamps.begin(); it != m_timestamps.end(); ++it)
			{
				TEdge tmp = {.src = m_src, .dst = m_dst, .tim = *it, .id=ID++};
				temporal_edges.push_back(tmp);
			}
			std::sort(temporal_edges.begin(), temporal_edges.end());
			return temporal_edges;
		} 
		*/
};
struct hash_pair
{ 
	template <class T1, class T2> 
	size_t operator()(const std::pair<T1, T2>& p) const
	{ 
		auto hash1 = std::hash<T1>{}(p.first); 
		auto hash2 = std::hash<T2>{}(p.second); 
		size_t seed = 0;
		boost::hash_combine(seed, hash1);
		boost::hash_combine(seed, hash2);
		return seed;
		//return hash1 ^ hash2; 
	} 
};

struct hash_vector
{
	template <class T1> 
	size_t operator()(const std::vector<T1>& v) const
	{
		std::hash<T1> hasher;
		size_t seed = 0;
		for(T1 el : v)
		{
			seed ^= hasher(el) + 0x9e3779b9 + (seed<<6) + (seed>>2);
		}
		return seed;
	}
};



class StaticGraph{
	private:
		//std::unordered_map<std::pair<node_t, node_t>, std::vector<timestamp_t>, boost::hash<std::pair<node_t, node_t>> > m_edges;
		std::unordered_map<std::pair<node_t, node_t>, std::vector<timestamp_t>, hash_pair> m_edges;

	public:
		StaticGraph() {};
		StaticGraph(const std::vector<TEdge>& temporal_graph) 
		{
			for(auto& edge : temporal_graph)
			{
				std::pair<node_t, node_t> srcdst = std::make_pair(edge.src, edge.dst);
				if(m_edges.find(srcdst) != m_edges.end())
				{
					m_edges.at(srcdst).push_back(edge.tim);
				}
				else
				{
					std::vector<timestamp_t> mytime{edge.tim};
					m_edges.insert(std::make_pair(srcdst, mytime));
				}
			}
		};
		long long getNumberOfEdges() { return m_edges.size(); };
		std::vector<std::pair<std::pair<node_t, node_t>, std::vector<timestamp_t>> > getFullVector() 
		{
			std::vector<std::pair<std::pair<node_t, node_t>, std::vector<timestamp_t>> > toreturn;
			for(auto& entry : m_edges)
			{
				toreturn.push_back(entry);
			}
			for(auto& element : toreturn)
			{
				std::sort(element.second.begin(), element.second.end());
			}
			return toreturn;
		}

};

class Permutation{
	private:
		std::random_device m_rd;
		std::mt19937 m_g;
	public:
		Permutation() : m_g{m_rd()} {};
		
		void GetPermutation(std::vector<int>& input){
			std::shuffle(input.begin(), input.end(), m_g);
		}
};

#endif  // utiilities_h
