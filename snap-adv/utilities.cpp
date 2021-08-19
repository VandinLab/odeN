#include "utilities.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <queue>

int Utils::LoadEdgeList(const std::string& filenameG, std::vector<TEdge>& edges_, std::vector<size_t>& temporal_degrees_){
	std::ifstream fileG(filenameG);
	if(fileG.is_open())
	{
		std::string line;
		long int ID = 0;
		int self_edges = 0;
		bool first = true;
		long long prev = 0;
		int max_node_id {0};
			while(getline(fileG, line))
			{
				std::istringstream iss(line);
				std::vector<std::string> results(std::istream_iterator<
							std::string>{iss},std::istream_iterator<std::string>());
				if(results.size() != 3)
				{
					std::cerr << "ERROR, dataset is not in format src dst timestamp" << std::endl;
					exit(1);
				}
				else
				{
					int src =stoi(results[0]);
					int dst =stoi(results[1]);
					max_node_id = std::max(std::max(max_node_id, src), dst);
					long long int timestamp = stoll(results[2]);
					if(src != dst)
					{
						TEdge edge = {.src=src, .dst=dst, .tim=timestamp, .id=ID++};
						edges_.push_back(edge);
					}
					else
						self_edges++;
				}
			}
			std::cout << "removed edges is: " << self_edges << std::endl;
			std::cout << "maximum id of a node: " << max_node_id << std::endl;
	}
	fileG.close();
	return 0;
}

int Utils::LoadTemporalMotif(const std::string& filenameM, std::vector<std::pair<int, int>>& edgesM_, int& Vm_, int& l_){

	std::ifstream fileM(filenameM);
	if(fileM.is_open())
	{
		std::string line;
		int max_nodes = 0;
		while(getline(fileM, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> results(std::istream_iterator<
						std::string>{iss},std::istream_iterator<std::string>());
			if(results.size() != 3)
			{
				std::cerr << "ERROR, motif is not in format src dst order" << std::endl;
				exit(1);
			}
			else
			{
				// Assume the edges are already ordered
				int src =stoi(results[0]);
				int dst =stoi(results[1]);
				if(src != dst)
				{
      				edgesM_.push_back(std::make_pair(src, dst));
	  				l_++; // counting the number of edges
					max_nodes = std::max(max_nodes, src + 1);
					max_nodes = std::max(max_nodes, dst + 1);
				}
			}
		}
		Vm_ = max_nodes;
	}

	fileM.close();

	return 0;
}

int Utils::LoadEdgeList(const std::string& filenameG, TVec< THash<TInt, TVec<TInt64> >> &temporal_data, std::vector<size_t>& temporal_degrees_, long long int& totedges, timestamp_t& timespan){
	std::ifstream fileG(filenameG);
	timestamp_t t1;
	timestamp_t tm;
	const char* DELIM = " "; // TODO: Possibly add custom delims
	char currentString[150]; // Bound on the number of characters on a line: 149
	if(fileG.is_open())
	{
		//std::string line;
		long int ID = 0;
		int self_edges = 0;
		bool first = true;
		long long prev = 0;
		int max_node_id {0};
			//while(getline(fileG, line))
			while (fileG.getline(currentString, 150, '\n'))
			{
				char* chunk = strtok(currentString, DELIM);
				int src = atoi(chunk);
				chunk = strtok(NULL, DELIM);
				int dst = atoi(chunk);
				max_node_id = std::max(std::max(max_node_id, src), dst);
				chunk = strtok(NULL, DELIM);
				long long int timestamp = atoll(chunk); 
				chunk = strtok(NULL, DELIM);
				if(src != dst)
				{
					temporal_data[src](dst).Add(timestamp);
					//Trial
					temporal_degrees_[src]++;
					temporal_degrees_[dst]++;
					totedges++;
					tm = timestamp;
					if(first)
					{
						first = false;
						t1 = timestamp;
					}
				}
				else
					self_edges++;
				delete chunk;
				/*
				std::istringstream iss(line);
				std::vector<std::string> results(std::istream_iterator<
							std::string>{iss},std::istream_iterator<std::string>());
				if(results.size() != 3)
				{
					std::cerr << "ERROR, dataset is not in format src dst timestamp" << std::endl;
					continue;
				}
				else
				{
					int src =stoi(results[0]);
					int dst =stoi(results[1]);
					max_node_id = std::max(std::max(max_node_id, src), dst);
					long long int timestamp = stoll(results[2]);
					if(src != dst)
					{
						temporal_data[src](dst).Add(timestamp);
						//Trial
						temporal_degrees_[src]++;
						temporal_degrees_[dst]++;
						totedges++;
						tm = timestamp;
						if(first)
						{
							first = false;
							t1 = timestamp;
						}
					}
					else
						self_edges++;
				}
				*/
			}
			std::cout << "removed edges is: " << self_edges << std::endl;
			std::cout << "maximum id of a node: " << max_node_id << std::endl;
	}
	//delete DELIM;
	fileG.close();
	timespan = tm -t1;
	return 0;
}

int Utils::RemapNodes(std::vector<TEdge>& edgeList){
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
	return 0;
}

int Utils::PrintGraph(std::vector<TEdge>& edgeList){
	for(auto& edge : edgeList)
		std::cout << edge.src << " " << edge.dst << " " << edge.tim << '\n';
	return 0;
}


int Utils::GetTrianglesFromNode(PUNGraph graph, const int startingNode, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles){
	TUNGraph::TNodeI sNode = graph->GetNI(startingNode);
	//int bag{0};
	for(int i=0; i < sNode.GetDeg(); i++)
	{
		int firstNeighbor{ sNode.GetNbrNId(i) };
		TUNGraph::TNodeI secondNode = graph->GetNI(firstNeighbor);

		//TODO: improve this by cycling over the min degree e.g. as it is now or the original node if it has a lower degree
		for(int j=0; j < secondNode.GetDeg(); j++)
		{
			int secondNeighbor{ secondNode.GetNbrNId(j) };
			if( (secondNeighbor != startingNode) && (firstNeighbor < secondNeighbor) && (graph->IsEdge(startingNode, secondNeighbor)) )
			{
				//if(bag < 1000)
				//{
						std::vector<std::pair<node_t, node_t>> triangle;
						triangle.push_back(std::make_pair<node_t, node_t>(startingNode, firstNeighbor) );
						triangle.push_back(std::make_pair<node_t, node_t>(firstNeighbor, secondNeighbor) );
						triangle.push_back(std::make_pair<node_t, node_t>(startingNode, secondNeighbor) );
						triangles.push_back(triangle);
						//bag++;
				//}
				//else
					//return 0;

				/*
				triangles.push_back(std::vector<std::pair<node_t, node_t>>{ std::make_pair<node_t, node_t>(startingNode, firstNeighbor) } );
				triangles.back().push_back(std::make_pair<node_t, node_t>(firstNeighbor, secondNeighbor));
				triangles.back().push_back(std::make_pair<node_t, node_t>(startingNode, secondNeighbor));
				*/
			}
		}
	}

	return 0;
}

int Utils::GetTrianglesFromEdge(const PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles)
{
	TUNGraph::TNodeI firstNode = graph->GetNI(sampledEdge.first);
	TUNGraph::TNodeI secondNode = graph->GetNI(sampledEdge.second);
	//bool first = (firstNode.GetDeg() < secondNode.GetDeg()) ? true : false;
	TUNGraph::TNodeI cycleNode = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.first) : graph->GetNI(sampledEdge.second);
	TUNGraph::TNodeI toMatch = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.second) : graph->GetNI(sampledEdge.first);

	for(int i=0; i < cycleNode.GetDeg(); i++)
	{
		int firstNeighbor{ cycleNode.GetNbrNId(i) };

		if( (firstNeighbor != toMatch.GetId()) && (graph->IsEdge(firstNeighbor, toMatch.GetId())) )
		{
			std::vector<std::pair<node_t, node_t>> triangle;
			//triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), toMatch.GetId()) );
			triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), firstNeighbor) );
			triangle.push_back(std::make_pair<node_t, node_t>(firstNeighbor, toMatch.GetId()) );
			triangles.push_back(triangle);

		}
	}

	return 0;

}

int Utils::GetTrianglesFromEdgeThreads(const PUNGraph& graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& triangles)
{
	TUNGraph::TNodeI firstNode = graph->GetNI(sampledEdge.first);
	TUNGraph::TNodeI secondNode = graph->GetNI(sampledEdge.second);
	//bool first = (firstNode.GetDeg() < secondNode.GetDeg()) ? true : false;
	TUNGraph::TNodeI cycleNode = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.first) : graph->GetNI(sampledEdge.second);
	TUNGraph::TNodeI toMatch = (firstNode.GetDeg() < secondNode.GetDeg()) ? graph->GetNI(sampledEdge.second) : graph->GetNI(sampledEdge.first);

	for(int i=0; i < cycleNode.GetDeg(); i++)
	{
		int firstNeighbor{ cycleNode.GetNbrNId(i) };

		if( (firstNeighbor != toMatch.GetId()) && (graph->IsEdge(firstNeighbor, toMatch.GetId())) )
		{
			std::vector<std::pair<node_t, node_t>> triangle;
			//triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), toMatch.GetId()) );
			triangle.push_back(std::make_pair<node_t, node_t>(cycleNode.GetId(), firstNeighbor) );
			triangle.push_back(std::make_pair<node_t, node_t>(firstNeighbor, toMatch.GetId()) );
			triangles.push_back(triangle);

		}
	}

	return 0;

}


int Utils::GetWedgesFromEdge(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& wedges)
{
	TUNGraph::TNodeI firstNode = graph->GetNI(sampledEdge.first);
	TUNGraph::TNodeI secondNode = graph->GetNI(sampledEdge.second);

	for(int i=0; i < firstNode.GetDeg(); i++)
	{
		int firstNeighbor{ firstNode.GetNbrNId(i) };

		if( (firstNeighbor != sampledEdge.second))
		{
			std::vector<std::pair<node_t, node_t>> wedge;
			//wedge.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(sampledEdge.first), static_cast<node_t>(sampledEdge.second)) );
			wedge.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(sampledEdge.first), firstNeighbor) );
			wedges.push_back(wedge);

		}
	}

	for(int i=0; i < secondNode.GetDeg(); i++)
	{
		int firstNeighbor{ secondNode.GetNbrNId(i) };

		if( (firstNeighbor != sampledEdge.first))
		{
			std::vector<std::pair<node_t, node_t>> wedge;
			//wedge.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(sampledEdge.first), static_cast<node_t>(sampledEdge.second)) );
			wedge.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(sampledEdge.second), firstNeighbor) );
			wedges.push_back(wedge);

		}
	}
	return 0;

}


int Utils::GetSquaresFromEdge(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& squares)
{
	node_t src = static_cast<node_t>(sampledEdge.first);
	node_t dst = static_cast<node_t>(sampledEdge.second);
	TUNGraph::TNodeI firstNode = graph->GetNI(src);
	TUNGraph::TNodeI secondNode = graph->GetNI(dst);

	for(int i=0; i < firstNode.GetDeg(); i++)
	{
		int firstNeighbor{ firstNode.GetNbrNId(i) };
		TUNGraph::TNodeI firstNeighborNode = graph->GetNI(firstNeighbor);

		if(firstNeighbor != dst)
		{
			for(int j=0; j < firstNeighborNode.GetDeg(); j++)
			{
				int secondNeighbor{ firstNeighborNode.GetNbrNId(j) };

				if( (firstNeighbor < secondNeighbor) && (secondNeighbor != dst) && (secondNeighbor != src) && (graph->IsEdge(secondNeighbor, dst)) )
				{
					std::vector<std::pair<node_t, node_t>> square;
					//square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(firstNeighbor)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(secondNeighbor)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(secondNeighbor), static_cast<node_t>(dst)) );
					squares.push_back(square);
				}
			}
		}
	}

	for(int i=0; i < secondNode.GetDeg(); i++)
	{
		int firstNeighbor{ secondNode.GetNbrNId(i) };
		TUNGraph::TNodeI firstNeighborNode = graph->GetNI(firstNeighbor);

		if(firstNeighbor != src)
		{
			for(int j=0; j < firstNeighborNode.GetDeg(); j++)
			{
				int secondNeighbor{ firstNeighborNode.GetNbrNId(j) };

				if( (firstNeighbor < secondNeighbor) && (secondNeighbor != src) && (secondNeighbor != dst) && (graph->IsEdge(secondNeighbor, src)) )
				{
					std::vector<std::pair<node_t, node_t>> square;
					//square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(dst), static_cast<node_t>(firstNeighbor)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(secondNeighbor)) );
					square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(secondNeighbor), static_cast<node_t>(src)) );
					squares.push_back(square);
				}
			}
		}
	}
	return 0;


}

int Utils::GetSquaresFromEdgeMin(PUNGraph graph, const std::pair<node_t,node_t>& sampledEdge, std::vector<std::vector<std::pair<node_t, node_t> > >& squares)
{
	node_t src = static_cast<node_t>(sampledEdge.first);
	node_t dst = static_cast<node_t>(sampledEdge.second);
	TUNGraph::TNodeI firstNode = graph->GetNI(src);
	TUNGraph::TNodeI secondNode = graph->GetNI(dst);
	bool minFirst = (firstNode.GetDeg() < secondNode.GetDeg());


	if(minFirst)
	{
		for(int i=0; i < firstNode.GetDeg(); i++)
		{
			int firstNeighbor{ firstNode.GetNbrNId(i) };
			TUNGraph::TNodeI firstNeighborNode = graph->GetNI(firstNeighbor);

			if(firstNeighbor != dst)
			{
				for(int j=0; j < firstNeighborNode.GetDeg(); j++)
				{
					int secondNeighbor{ firstNeighborNode.GetNbrNId(j) };

					if( (secondNeighbor != dst) && (secondNeighbor != src) && (graph->IsEdge(secondNeighbor, dst)) )
					{
						std::vector<std::pair<node_t, node_t>> square;
						//square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(firstNeighbor)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(secondNeighbor)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(secondNeighbor), static_cast<node_t>(dst)) );
						squares.push_back(square);
					}
				}
			}
		}
	}
	else
	{
		for(int i=0; i < secondNode.GetDeg(); i++)
		{
			int firstNeighbor{ secondNode.GetNbrNId(i) };
			TUNGraph::TNodeI firstNeighborNode = graph->GetNI(firstNeighbor);

			if(firstNeighbor != src)
			{
				for(int j=0; j < firstNeighborNode.GetDeg(); j++)
				{
					int secondNeighbor{ firstNeighborNode.GetNbrNId(j) };

					if((secondNeighbor != src) && (secondNeighbor != dst) && (graph->IsEdge(secondNeighbor, src)) )
					{
						std::vector<std::pair<node_t, node_t>> square;
						//square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(dst), static_cast<node_t>(firstNeighbor)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(firstNeighbor), static_cast<node_t>(secondNeighbor)) );
						square.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(secondNeighbor), static_cast<node_t>(src)) );
						squares.push_back(square);
					}
				}
			}
		}
	}
	return 0;
}

