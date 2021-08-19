#include "datastructures.h"

#include <unordered_map>
#include <vector>
#include <algorithm>

void DataStructures::Datum::InitializeStructures(std::vector< TEdge >& edgeList, int Vm) { 
	std::unordered_map<node_t, node_t> remap;
	node_t id = 0;
	for (auto e : edgeList) {
		if (!remap.count(e.src)) {
			remap[e.src] = id++;
		}
		if (!remap.count(e.dst)) {
			remap[e.dst] = id++;
		}
	}
	for (auto& e : edgeList) {
		e.src = remap[e.src];
		e.dst = remap[e.dst];
	}

	/*
	max_node_id = id;
	*/
	sort(edgeList.begin(), edgeList.end());

	// FIXME: only up to 2B edges otherwise error!!
	for (int i = 0; i < (int) edgeList.size(); i++) {
		edgeList[i].id = i;
	}

	adj_list.resize(id);
	revadj_list.resize(id);
	adjMap.resize(id);
	for (const auto& edge : edgeList) {
		adj_list[edge.src].push_back(edge);
		revadj_list[edge.dst].push_back(edge);
		adjMap[edge.src][edge.dst].push_back(edge);
	}
	edgeCount.resize(id, 0);
	mapGM.resize(id, -1);
	mapMG.resize(Vm, -1);
	//TODO: Remove this
	for(auto el : remap)
	{
		remaptmp.emplace(el.second, el.first);
	}
} 

