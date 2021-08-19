#ifndef custom_datastructures_h
#define custom_datastructures_h

#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/functional/hash.hpp>
#include <unordered_set>

// Type defintions, they are used to improve readability of thew code by employing domain specific names
using timestamp_t = long long int;
using node_t = long long int;

using staticEdge_t = std::pair<node_t, node_t>; // directed static edge <src, dst>

// Datastructures used in the project
namespace DataStructures {

	struct TEdge {
		node_t src;
   		node_t dst;
		timestamp_t tim;
   		long long id;
		const bool operator<(const TEdge& o) const {
			if (o.tim != this->tim) return this->tim < o.tim;
			if (o.id != this->id) return this->id < o.id;
			return false; 
		}
		friend bool operator<(const TEdge &e, const timestamp_t &time) 
		{
			return (e.tim < time);
		}
		friend bool operator<(const timestamp_t &time ,const TEdge &e) 
		{
			return (time < e.tim);
		}
		const bool operator==(const TEdge& e) const{
				if( (e.tim==tim) && (e.src==src)  && (e.dst==dst) ) return true;
				return false;
		}
	
		friend std::ostream& operator<<(std::ostream& in, const TEdge& o) {
			in << "(" << o.src << "," << o.dst << "," << o.tim << "," << o.id << ")";
			return in;
		}
	};


	// Used for 2-delta patch, every partition is an interval of some length, patch = true iff the partition is a patch
	struct Partition
	{
		unsigned int start;
		timestamp_t tStart;
		unsigned int end;
		timestamp_t tEnd;
		bool patch;
	
		const bool operator<(const Partition& p) const
		{
			if (tStart != p.tStart)
				return tStart < p.tStart;
			if ((tStart == p.tStart) && (tEnd != p.tEnd))
				return tEnd < p.tEnd;
			return false;
		}
	};

	struct Datum {
		// Vertex based data structures
		std::vector<int> edgeCount;
    	std::vector<int> mapGM;
    	std::vector<int> mapMG;
    	// Edge index based adjacency list
		std::vector< std::vector< TEdge > > adj_list;
		std::vector< std::vector< TEdge > > revadj_list;
		std::vector< std::unordered_map<int, std::vector<TEdge> > > adjMap;
		// for each node u keep the list of indexes of (u,v,t) or (v,u,t)
		std::vector< std::vector<int> > node_edge_index_;
		node_t max_node_id; // This is also the number of nodes in edgeList
		std::unordered_map<node_t, node_t> remaptmp;

		//void InitializeStructures(std::vector< TEdge >& edgeList, int Vm); 
		void InitializeStructures(std::vector< TEdge >& edgeList, int Vm); 
	};

	//inline size_t HashEdge(const std::map<node_t, node_t>& edges){
	inline size_t HashEdge(const std::vector<std::pair<node_t, node_t>>& edges){
		std::size_t seed = 0;	
			for (auto& entry : edges)
			{
				boost::hash_combine(seed, entry.first * 2654435761);
				boost::hash_combine(seed, entry.second * 2654435761);
			}
			return seed;
		}

	template < typename Graph1, typename Graph2 > 
	struct vf2_print_custom_callback
	{

			vf2_print_custom_callback(const Graph1& graph1, const Graph2& graph2)
			: graph1_(graph1), graph2_(graph2), c{0}, count{ &c }, occurr{ new std::vector<std::vector<std::pair<node_t, node_t>>> } , hashes{ new std::unordered_set<size_t> }
			{
				//motif.push_back(std::make_pair<node_t,node_t>(0,1));
				//motif.push_back(std::make_pair<node_t,node_t>(1,0));
				BGL_FORALL_EDGES_T(ed, graph1_, Graph1)
					motif.push_back(std::make_pair<node_t, node_t>(source(ed,graph1_), target(ed,graph1_)));
				/*
				BGL_FORALL_EDGES_T(ed, graph1_, Graph1)
					std::cout << "Motif: src boost: " <<  source(ed,graph1_) << " dst boost:" << target(ed,graph1_) << '\n';
				BGL_FORALL_EDGES_T(ed, graph2_, Graph2)
					std::cout <<  source(ed,graph2_) << " " << target(ed,graph2_) << '\n';
				*/
			}

			template < typename CorrespondenceMap1To2, typename CorrespondenceMap2To1 >
			bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) const
			{

				//std::unordered_map<node_t, node_t> mapping;
				// Print (sub)graph isomorphism map
				/*
				BGL_FORALL_VERTICES_T(v, graph1_, Graph1)
					//get(vertex_index_t(), graph1_, v) << ", "  << get(vertex_index_t(), graph2_, get(f, v)) << ") \n" ;
					mapping.insert(std::make_pair<node_t, node_t>( boost::get(boost::vertex_index_t(), graph1_, v), boost::get(boost::vertex_index_t(), graph2_, boost::get(f, v)) ) ) ; 
				*/
				//std::map<node_t, node_t> occurrence;
				std::vector<std::pair<node_t, node_t>> occurrence;
				BGL_FORALL_EDGES_T(ed, graph1_, Graph1)
					occurrence.push_back(std::make_pair<node_t,node_t>(boost::get(boost::vertex_index_t(), graph2_, boost::get(f, source(ed, graph1_)) ), boost::get(boost::vertex_index_t(), graph2_, boost::get(f, target(ed, graph1_)) ) ) );
					//occurrence.insert(std::make_pair<node_t,node_t>(boost::get(f, source(ed, graph1_)), boost::get(f, target(ed, graph1_)) ));
					//occurrence.push_back(std::make_pair<node_t,node_t>(boost::get(f, source(ed, graph1_)), boost::get(f, target(ed, graph1_)) ));

				sort(occurrence.begin(), occurrence.end());
				
				for( auto& el : occurrence)
				{
					std::cout << "src: " << el.first << " dst: " << el.second << '\n';
				}
				std::cout << "----------------\n";
				
				//int i= 0;
				/*
				for( const staticEdge_t &edge : motif)
				{
					//occurrence[i++] = std::make_pair<node_t, node_t>(src, dst);
					occurrence.insert(std::make_pair<node_t,node_t>(int(mapping.at(edge.first)), int(mapping.at(edge.second)))) ;

				}
				*/
				auto myhash = HashEdge(std::cref(occurrence));
				//std::cout << myhash << '\n';
				if(hashes->find(myhash) != hashes->end())
				{
				}
				else
				{
						/*
					for( const staticEdge_t &edge : motif)
					{
						std::cout << int(mapping.at(edge.first)) << ", " <<int(mapping.at(edge.second)) << '\n';
					}
					*/
					*count = (*(count)) + 1;
					occurr->push_back(occurrence);
					hashes->insert(myhash);
				}
				return true;
			}
			inline const std::vector<std::vector<std::pair<node_t, node_t>>> & getVector(){ return std::cref((*occurr)); }
			inline const double getCount(){ return c; }
			

		private:
			const Graph1& graph1_;
			const Graph2& graph2_;
			std::vector<staticEdge_t> motif;
			mutable std::vector<std::vector<std::pair<node_t,node_t>>>* occurr;
			mutable double* count;
			mutable double c;
			mutable std::unordered_set<size_t>* hashes;
		};

		struct hash_pair
		{ 
			template <class T1, class T2> 
			size_t operator()(const std::pair<T1, T2>& p) const
			{ 
				auto hash1 = std::hash<T1>{}(p.first); 
				auto hash2 = std::hash<T2>{}(p.second); 
				size_t seed = 0;
				boost::hash_combine(seed, hash1 * 2654435761);
				boost::hash_combine(seed, hash2 * 2654435761);
				return seed;
				//return hash1 ^ hash2; 
			} 
		};

}


#endif // custom_datastructures_h
