
/* Lemon stuff */
#include <lemon/list_graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/vf2pp.h>

#include "Snap.h"

#include "utilities.h"
#include "temporalmotifstaticsampler.h"
#include <random>
#include "ctpl.h"
#include <mutex>

#include <vector>
#include <algorithm>
#include <iomanip>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/functional/hash.hpp>
#include <chrono>
#include <unordered_set>
#include <numeric>


typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS> graph_type;

///////////////////////////////////////////////////////////////////////////////

TempMotifSamplerStatic::TempMotifSamplerStatic(const TStr& filenameG, const TStr& filenameM) { 
	static_motif_ = TSnap::LoadEdgeList<PNGraph>(filenameM, 0, 1, ' '); 
	static_undirected_graph_ = TSnap::LoadEdgeList<PUNGraph>(filenameG, 0, 1, ' '); 
	static_undirected_motif_ = TSnap::LoadEdgeList<PUNGraph>(filenameM, 0, 1, ' '); 
	staticEdgesMotif_m = static_undirected_motif_->GetEdges();
	std::cout << "loaded undirected static MOTIF has: " << static_undirected_motif_->GetNodes() << " nodes and: " << static_undirected_motif_->GetEdges() <<'\n';
	std::cout << "Loaded static undirected graph has: " << static_undirected_graph_->GetMxNId() << " max node id and: " << static_undirected_graph_->GetEdges() << " edges, " << static_undirected_graph_->GetNodes() << " nodes" << '\n' ;
	temporal_degrees_.resize(static_undirected_graph_->GetNodes(), 0);

	ems_ = static_motif_->GetEdges();
	
	temporal_data_ = TVec< THash<TInt, TVec<TInt64> >>(static_undirected_graph_->GetNodes());
	std::cout << "temporal data instatiated" << std::endl;

	m_ = 0;
	Utils::LoadEdgeList(filenameG.CStr(), temporal_data_, std::ref(temporal_degrees_), m_, timespan_m);	
	std::cout << "temporal data filled" << std::endl;
	
	}



TempMotifSamplerStatic::TempMotifSamplerStatic(std::string filenameG, std::string filenameM, bool parallel)
{
	Utils::LoadEdgeList(filenameG, std::ref(edges_), std::ref(temporal_degrees_));	

	//Utils::RemapNodes(std::ref(edges_));
	//Utils::PrintGraph(std::ref(edges_));
	//exit(0);

	//std::cout << "Number on nodes from struct : " << allgraph.max_node_id << '\n';

	std::cout << "READ FILE DONE" << std::endl;
	sort(edges_.begin(), edges_.end());
	std::cout << "Number of temporal edges: " << edges_.size() << '\n';
	StaticGraph sg(std::cref(edges_)); 
	std::cout << "Number of static edges: " << sg.getNumberOfEdges() << '\n';
	staticandtemporalgraph_ = sg.getFullVector();
	std::cout << "STRUCTURE BUILT CORRECTLY" << std::endl;
	size_t max_temporal_edges = 0;
	for(const auto& edge : staticandtemporalgraph_)
	{
		size_t num = edge.second.size();
		if(num > max_temporal_edges)
			max_temporal_edges = num;
	}
	std::cout << "--Maxedges-- mapped on a static one is: " << max_temporal_edges << std::endl;

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
}

int TempMotifSamplerStatic::setMotif(std::string filenameM)
{
	std::ifstream fileM(filenameM);
	edgesM_.clear();
	l_ = 0;
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

double TempMotifSamplerStatic::get_wall_time()
{
	struct timeval time;
	if (gettimeofday(&time,NULL)){
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


double TempMotifSamplerStatic::EnumerateOnStaticGraph(int delta, int samp, int seed, double timeLimit, int ell, int triORsq)
{
	

	l_ = ell;
	std::random_device rd;
	std::mt19937 generator(seed);
	
	std::vector<double> mass; 
	double maxdegst { 0};
	double avgdegst {0};
	double alphamax {0};
	double denom{0};
	TUNGraph::TEdgeI edgeIt = static_undirected_graph_->BegEI();
	std::vector<std::pair<node_t, node_t>> edges_exp;
	TVec< THash<TInt, TInt64 >> probsEdges = TVec< THash<TInt, TInt64 >>(static_undirected_graph_->GetNodes());
	double firsttime = get_wall_time();
	for(long long int i=0; i< static_undirected_graph_->GetEdges(); i++)
	{
		node_t src = edgeIt.GetSrcNId();
		node_t dst = edgeIt.GetDstNId();

		(src < dst) ? edges_exp.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst))) : edges_exp.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(dst), static_cast<node_t>(src)));
		edgeIt++;
		double currmass { 0 };
		if(temporal_data_[src].IsKey(dst))
		{
			currmass += temporal_data_[src](dst).Len();
		}
		double mval = currmass;
		if(temporal_data_[dst].IsKey(src))
		{
			currmass += temporal_data_[dst](src).Len();
		}
		

		mass.push_back(static_cast<double>(currmass));
		(src < dst) ? probsEdges[src].AddDat(dst, i) : probsEdges[dst].AddDat(src, i);
	}

	std::discrete_distribution<long long int> importance(mass.begin(), mass.end());
	std::vector<double> sampling_probs{ importance.probabilities() };
	double randomtime = get_wall_time() - firsttime;

	double totestimate{ 0};
	int s{ samp };
	double timeEnumerateTri{ 0 };
	double timeBuildSample{ 0 };
	double timeEnumerateTemporal{ 0 };
	double timenow{ 0 };

	double timeSequence{ 0 };
	double timeIso{ 0 };
	double totInstances{ 0 };
	double totStaticInstances{ 0 };
	std::vector<int> numEdges;
	std::vector<int> staticInsts;
	int pruned{ 0 };
	std::unordered_map<std::uint_fast64_t, double> counts;
	double artifactNorm { 0 };
	long long samples1 = 0;
	double init = get_wall_time();
	while( (get_wall_time() - init) < timeLimit)
	{
			samples1++;
			long long int samplednode { importance(generator) };
			
			std::vector<std::vector<std::pair<node_t, node_t> > > static_graphlets;

			double timemine = get_wall_time();
			
			if(triORsq == 0)
				Utils::GetWedgesFromEdge(static_undirected_graph_, edges_exp[samplednode], std::ref(static_graphlets) );
			else if(triORsq == 1)
				Utils::GetTrianglesFromEdge(static_undirected_graph_, edges_exp[samplednode], std::ref(static_graphlets) );
			else if(triORsq == 2)
				Utils::GetSquaresFromEdgeMin(static_undirected_graph_, edges_exp[samplednode], std::ref(static_graphlets) );
			else if(triORsq == 3) // Single edge motif
			{
				std::vector<std::pair<node_t,node_t>> mytmpvec;
				mytmpvec.push_back(std::make_pair(0,1)); //placeholder
				static_graphlets.push_back(mytmpvec);
			}
			timeEnumerateTri += (get_wall_time()-timemine);
			timemine = get_wall_time() - timemine;

			if(static_graphlets.size() == 0)
				continue;

			timemine = get_wall_time();
			std::vector<std::pair<node_t, node_t>> staticEdgesInSample;
			double sampled{ 0 };
			
			totStaticInstances += static_cast<double>(static_graphlets.size());
			std::vector<std::pair<timestamp_t,size_t>> times; //vector of <timestamp, idx_static_edge_In_sample (which corresponds to the timestamp)> 

			//Computing on sampled edge to avoid this repeated computation
			std::vector<std::pair<timestamp_t,size_t>> timesSampledEdge;
			node_t sampledSrc = static_cast<node_t>(edges_exp[samplednode].first);
			node_t sampledDst = static_cast<node_t>(edges_exp[samplednode].second);
			std::pair<node_t, node_t> P1{ std::make_pair(sampledSrc, sampledDst) };
			std::pair<node_t, node_t> P2{ std::make_pair(sampledDst, sampledSrc) };
			bool exsistenceSD {false};
			bool exsistenceDS {false};
			if(temporal_data_[sampledSrc].IsKey(sampledDst))
			{
				exsistenceSD = true;
				const TVec<TInt64>& timestamps = temporal_data_[sampledSrc].GetDat(sampledDst);
				for(size_t j=0; j< timestamps.Len(); j++)
				{
					timesSampledEdge.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], 0));
				}
			}

			if(temporal_data_[sampledDst].IsKey(sampledSrc))
			{
				exsistenceDS = true;
				const TVec<TInt64>& timestamps = temporal_data_[sampledDst].GetDat(sampledSrc);
				for(size_t j=0; j< timestamps.Len(); j++)
				{
					if(exsistenceSD)
						timesSampledEdge.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], 1));
					else
						timesSampledEdge.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], 0));
				}
			}
			timeBuildSample += (get_wall_time() - timemine);


			for(size_t k{ 0 }; k < static_graphlets.size(); k++)
			{
				//timenow = get_wall_time();
				timemine = get_wall_time();
				staticEdgesInSample.clear();
				times.clear();

				if(exsistenceSD)
					staticEdgesInSample.push_back(P1);
				if(exsistenceDS)
					staticEdgesInSample.push_back(P2);

				times.insert(times.begin(), timesSampledEdge.begin(), timesSampledEdge.end());
				



				//Collecting graphlet edges
				long id{0};
					
				if(triORsq != 3)
				{
					for(const auto& pair : static_graphlets[k])
					{
						node_t src = pair.first;
						node_t dst = pair.second;
						if(temporal_data_[src].IsKey(dst))
						{
							staticEdgesInSample.push_back(std::make_pair<node_t,node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)));
							size_t currentStaticEdges{ staticEdgesInSample.size() };
							const TVec<TInt64>& timestamps = temporal_data_[src].GetDat(dst);
							for(size_t j=0; j< timestamps.Len(); j++)
							{
								times.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], currentStaticEdges-1));
							}
						}

						if(temporal_data_[dst].IsKey(src))
						{
							staticEdgesInSample.push_back(std::make_pair<node_t,node_t>(static_cast<node_t>(dst), static_cast<node_t>(src)));
							size_t currentStaticEdges{ staticEdgesInSample.size() };
							const TVec<TInt64>& timestamps = temporal_data_[dst].GetDat(src);
							for(size_t j=0; j< timestamps.Len(); j++)
							{
								times.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], currentStaticEdges-1));
							}
						}
					}
				}

					std::sort(times.begin(), times.end());
					timeBuildSample += (get_wall_time() - timemine);
					
					timemine = get_wall_time();
					bool foundCandidate{ false };

					// Pruning
					for(long long int t{ 0 }; t < static_cast<long long int>(times.size())-l_+1; t++)
					{
						if( ((times[t+l_-1].first - times[t].first) <= delta ))
						{
								foundCandidate = true;
								break;
						}
					}
					timeEnumerateTemporal += (get_wall_time() - timemine);
					if(!foundCandidate)
					{
						pruned++;
					}
					else
					{
							timemine = get_wall_time();
							double graphlet_prob{0};
							graphlet_prob += sampling_probs[samplednode];
							

							double weight { static_cast<double>(1./graphlet_prob) }; 

							double ressamp  = FastCountSequence(delta, std::ref(times), std::ref(staticEdgesInSample), timeSequence, timeIso, std::ref(counts), graphlet_prob);
							totInstances += ressamp;
							numEdges.push_back(static_cast<int>(times.size()));

							ressamp = ressamp * weight;
	
							sampled++;

							totestimate += ressamp;
							timeEnumerateTemporal += (get_wall_time() - timemine);
					}

			} // TODO: Delete tmp pointers!!
		
	}

	double tottime = get_wall_time() - firsttime;
	std::cout.precision(2);
	for(auto& el : counts)
	{
		std::cout << " count is: " <<  std::fixed << el.second/(1.*samples1*static_undirected_motif_->GetEdges()) << " motif is: " << std::hex << el.first << std::dec;
		std::cout << '\n';
	}
	std::cout << "T0 time spent in building prob distrib: " << std::fixed << randomtime  << '\n';
	std::cout << "T1 time spent in static enumeration: " << std::fixed << timeEnumerateTri  << '\n';
	std::cout << "T2 time spent in building sample: " << std::fixed << timeBuildSample  << '\n';
	std::cout << "T3 time spent in temporal enumeration: " << std::fixed << timeEnumerateTemporal  << '\n';
	std::cout << "T3-a time spent in temporal enumeration when enumerating sequences: " << std::fixed << timeSequence  << '\n';
	std::cout << "T3-b time spent in temporal enumeration with isomorphism: " << std::fixed << timeIso  << '\n';
	std::cout << "T4 time spent in whole precedure: " << std::fixed << tottime  << '\n';

	std::cout << "samples s examined: " << samples1 << '\n';
	return totestimate/s;
}


double TempMotifSamplerStatic::EnumerateOnStaticGraphParallel(int delta, int samp, int seed, double timeLimit, int ell, int triORsq, int nthreads)
{
	

	l_ = ell;
	double firsttime = get_wall_time();
	std::random_device rd;
	std::mt19937 generator(seed);
	
	std::vector<double> mass; 
	std::vector<double> mass2v; 
	std::vector<double> mass3v; 
	TUNGraph::TEdgeI edgeIt = static_undirected_graph_->BegEI();
	std::vector<std::pair<node_t, node_t>> edges_exp;
	TVec< THash<TInt, TInt64 >> probsEdges = TVec< THash<TInt, TInt64 >>(static_undirected_graph_->GetNodes());
	for(long long int i=0; i< static_undirected_graph_->GetEdges(); i++)
	{
		node_t src = edgeIt.GetSrcNId();
		node_t dst = edgeIt.GetDstNId();
		(src < dst) ? edges_exp.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst))) : edges_exp.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(dst), static_cast<node_t>(src)));
		edgeIt++;
		double currmass { 0 };
		if(temporal_data_[src].IsKey(dst))
		{
			currmass += temporal_data_[src](dst).Len();
		}
		if(temporal_data_[dst].IsKey(src))
		{
			currmass += temporal_data_[dst](src).Len();
		}

		mass.push_back(static_cast<double>(currmass));
		(src < dst) ? probsEdges[src].AddDat(dst, i) : probsEdges[dst].AddDat(src, i);
	}
	
	std::discrete_distribution<long long int> importance(mass.begin(), mass.end());
	std::vector<double> sampling_probs{ importance.probabilities() };
	double randomtime = get_wall_time() - firsttime;

	double totestimate{ 0};
	int s{ samp };
	double timeEnumerateTri{ 0 };
	double timeBuildSample{ 0 };
	double timeEnumerateTemporal{ 0 };
	double timenow{ 0 };
	ctpl::thread_pool p(nthreads);
	std::mutex insertSolution;
	std::mutex printLock;
	std::mutex useRandom;
	std::vector<std::unordered_map<std::uint_fast64_t, double>> results;
	std::vector<std::future<bool>> futures;

	std::vector<int> numEdges;
	std::vector<int> staticInsts;
	long long samples1 = 0;
	double init = get_wall_time();

	//DEFINING LAMBDA
	auto sampleAndEstimate{ [&](int ID) mutable -> bool {

	std::unordered_map<std::uint_fast64_t, double> counts;
	useRandom.lock();
	long long int samplednode { importance(generator) };
	useRandom.unlock();

	
	std::vector<std::vector<std::pair<node_t, node_t> > > static_graphlets;

	if(triORsq == 0)
		Utils::GetWedgesFromEdge(static_undirected_graph_, edges_exp[samplednode], std::ref(static_graphlets) );
	else if(triORsq == 1)
		Utils::GetTrianglesFromEdgeThreads(std::cref(static_undirected_graph_), std::cref(edges_exp[samplednode]), std::ref(static_graphlets) );
	else if(triORsq == 2)
		Utils::GetSquaresFromEdgeMin(static_undirected_graph_, edges_exp[samplednode], std::ref(static_graphlets) );
	
	if(static_graphlets.size() == 0)
	{
		
		return true;
	}

	std::vector<std::pair<node_t, node_t>> staticEdgesInSample;
	double sampled{ 0 };
	
	std::vector<std::pair<timestamp_t,size_t>> times; //vector of <timestamp, idx_static_edge_In_sample (which corresponds to the timestamp)> 

	//Computing on sampled edge to avoid this repeated computation
	std::vector<std::pair<timestamp_t,size_t>> timesSampledEdge;
	node_t sampledSrc = static_cast<node_t>(edges_exp[samplednode].first);
	node_t sampledDst = static_cast<node_t>(edges_exp[samplednode].second);
	std::pair<node_t, node_t> P1{ std::make_pair(sampledSrc, sampledDst) };
	std::pair<node_t, node_t> P2{ std::make_pair(sampledDst, sampledSrc) };
	bool exsistenceSD {false};
	bool exsistenceDS {false};
	if(temporal_data_[sampledSrc].IsKey(sampledDst))
	{
		exsistenceSD = true;
		const TVec<TInt64>& timestamps = temporal_data_[sampledSrc].GetDat(sampledDst);
		for(size_t j=0; j< timestamps.Len(); j++)
		{
			timesSampledEdge.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], 0));
		}
	}

	if(temporal_data_[sampledDst].IsKey(sampledSrc))
	{
		exsistenceDS = true;
		const TVec<TInt64>& timestamps = temporal_data_[sampledDst].GetDat(sampledSrc);
		for(size_t j=0; j< timestamps.Len(); j++)
		{
			if(exsistenceSD)
				timesSampledEdge.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], 1));
			else
				timesSampledEdge.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], 0));
		}
	}


	for(size_t k{ 0 }; k < static_graphlets.size(); k++)
	{
		staticEdgesInSample.clear();
		times.clear();

		if(exsistenceSD)
			staticEdgesInSample.push_back(P1);
		if(exsistenceDS)
			staticEdgesInSample.push_back(P2);

		times.insert(times.begin(), timesSampledEdge.begin(), timesSampledEdge.end());
		

		//Collecting graphlet edges
		long id{0};
			
			for(const auto& pair : static_graphlets[k])
			{
				node_t src = pair.first;
				node_t dst = pair.second;
				if(temporal_data_[src].IsKey(dst))
				{
					staticEdgesInSample.push_back(std::make_pair<node_t,node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)));
					size_t currentStaticEdges{ staticEdgesInSample.size() };
					const TVec<TInt64>& timestamps = temporal_data_[src].GetDat(dst);
					for(size_t j=0; j< timestamps.Len(); j++)
					{
						times.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], currentStaticEdges-1));
					}
				}

				if(temporal_data_[dst].IsKey(src))
				{
					staticEdgesInSample.push_back(std::make_pair<node_t,node_t>(static_cast<node_t>(dst), static_cast<node_t>(src)));
					size_t currentStaticEdges{ staticEdgesInSample.size() };
					const TVec<TInt64>& timestamps = temporal_data_[dst].GetDat(src);
					for(size_t j=0; j< timestamps.Len(); j++)
					{
						times.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], currentStaticEdges-1));
					}
				}
			}

			std::sort(times.begin(), times.end());
			
			bool foundCandidate{ false };

			for(long long int t{ 0 }; t < static_cast<long long int>(times.size())-l_+1; t++)
			{
				if( ((times[t+l_-1].first - times[t].first) <= delta ))
				{
						foundCandidate = true;
						break;
				}
			}
			if(!foundCandidate)
			{
				continue;
			}
			else
			{
					double timeSequence{0};
					double timeIso{0};
					//Computing probability and checking if the selected vertex is in the graphlet
					double graphlet_prob{0};
					graphlet_prob += sampling_probs[samplednode];
					
					double ressamp  = FastCountSequence(delta, std::ref(times), std::ref(staticEdgesInSample), timeSequence, timeIso, std::ref(counts), graphlet_prob);
			}

	} 
	if(counts.size() != 0)
	{
		insertSolution.lock();
		results.push_back(counts);
		insertSolution.unlock();
	}
	return true;
	} };


	for(long long int ii{ 0 }; ii < s; ii++)
	{
			futures.push_back(p.push(std::ref(sampleAndEstimate)));
	}

	for(int i{ 0}; i < futures.size(); i++)
	{
		bool correct = futures[i].get();
		if(!correct)
			std::cout << "Thread: " << " experienced an error!" << std::endl;
	}
	std::unordered_map<std::uint_fast64_t, double> counts;
	for(auto& map : results)
	{
		for(auto& pair : map)
		{
			if(counts.find(pair.first) != counts.end())
				counts[pair.first] += pair.second;
			else
				counts[pair.first] = pair.second;
		}
	}
				
	

	double tottime = get_wall_time() - firsttime;
	std::cout.precision(2);
	std::cout << "Time to build the probability distribution: " << randomtime << '\n';
	std::cout << "Innercycle launching time: " << tottime << '\n';
	for(auto& el : counts)
	{
		std::cout << " count is: " <<  std::fixed << el.second/(1.*s*static_undirected_motif_->GetEdges()) << " motif is: " << std::hex << el.first << std::dec;
		std::cout << '\n';
	}
	return 0;
	
}




double TempMotifSamplerStatic::FastCountSequence(int delta, std::vector<std::pair<timestamp_t, size_t>>& tmp, std::vector<std::pair<node_t, node_t>>& staticEdgesInSample, double& timeSequence, double& timeIso, std::unordered_map<std::uint_fast64_t, double>& counts, double instanceProbability)
{
	//double startT{ get_wall_time() };
	//std::sort(tmp.begin(), tmp.end());
	
	/*
	std::cout << "\n--------------------------------\n";
	for(size_t i{ 0 }; i < tmp->size(); i++)
		std::cout << "source: " << (*tmp)[i].src << " dst: " << (*tmp)[i].dst << " tim: " << (*tmp)[i].tim << '\n';
	std::cout << "--------------------------------\n";
	for(auto& pair : staticEdgesInSample)
		std::cout << "source: " << pair.first << " dst: " << pair.second << '\n';
	*/
	
	//std::unordered_map<std::pair<node_t,node_t>, int, DataStructures::hash_pair> hashID; //hashmap: <edge, edgeID>
	std::unordered_map<std::pair<node_t,node_t>, std::uint_fast64_t, DataStructures::hash_pair> hashID; //hashmap: <edge, edgeID>
	std::vector<std::pair<node_t,node_t>> IDReverse; //Needed to reconstruct the motif graph
	std::uint_fast64_t ID { 1 }; // WARNING: Up to 15 static edges ranging from 1 to 15 otherwise enlarge edgeMask or use integer >64 bits
	for(auto tuple : staticEdgesInSample) // cycling over the static directed edges of the sample
	{
		IDReverse.push_back(tuple);
		hashID[tuple] = ID++;
	}
	std::vector<std::unordered_set<std::uint_fast64_t>> keysLength; // each entry i=0,...,l-1 keeps the keys of length i+1 that have count > 0; 
	for(int i{ 0}; i < l_; i++)
	{
		std::unordered_set<std::uint_fast64_t> newset;
		keysLength.push_back(newset);
	}
	std::unordered_map<std::uint_fast64_t, long long> hashCounts; //hasmap["sequence"] = "count of the sequence"
	size_t start { 0 };
	for(const auto& time : tmp )
	{
		while( (time.first - tmp[start].first) > delta )
		{
		
			std::uint_fast64_t decID { hashID[staticEdgesInSample[tmp[start].second]] } ;
			DecrementCounts(decID, std::ref(hashCounts), std::ref(keysLength));
			start++;
		}

		std::uint_fast64_t incID { hashID[staticEdgesInSample[time.second]] } ;
		IncrementCounts(incID, std::ref(hashCounts), std::ref(keysLength));

	}

	std::vector<std::pair<node_t,node_t>> thisSequence;	
	double result{ 0 };
	for(auto& key : keysLength[l_-1])
	{
		thisSequence.clear();
		for(std::uint_fast64_t i{ 0 }; i < static_cast<std::uint_fast64_t>(l_); i++)
		{
			int keyedge = getDigit(key, l_-i);
			std::pair<node_t, node_t> edge{ IDReverse[keyedge-1] };	
			thisSequence.push_back(edge);
		}
		int edgeStatic = GetStaticEdgesInstance(std::ref(thisSequence));
		//NOTE: test for isomorphism here is not needed since we sampled only static graphs isomorphic to H
		// thus by the algorithm above each static projection of a motif can have up to |E_H| edges, but we know
		// that we require exactly |E_H| edges otherwise for sure the two motifs are non isomorphic!
		if(edgeStatic == staticEdgesMotif_m) // The instance is isomorphic to the given static graph (?) || here one can implement the constraint on q
		{	
			uint_fast64_t encode = ComputeEncoding(std::ref(thisSequence));
			if(counts.find(encode) == counts.end())
			{
				counts[encode] = (hashCounts[key] * 1/instanceProbability) ;
			}
			else
				counts[encode] += (hashCounts[key] * 1/instanceProbability);
		}
		//else
		//	continue;
	}
	
	return result;
}


bool TempMotifSamplerStatic::testIsomorphism(std::vector<std::pair<node_t, node_t>>& nodeSequence, std::vector<std::pair<int, int>>& motif)
{
	std::vector<node_t> mapSequence(Vm_, -1);
	int currentEdge{ 0 };
	std::unordered_set<node_t> alreadyMapped;
	for(const auto& edge : motif)
	{
		int srcM = edge.first;
		int dstM = edge.second;
		node_t srcSeq = nodeSequence[currentEdge].first;
		node_t dstSeq = nodeSequence[currentEdge].second;
		if(mapSequence[srcM] == -1) // node not mapped
		{
			if(alreadyMapped.find(srcSeq) == alreadyMapped.end())
			{
				mapSequence[srcM] = srcSeq;
				alreadyMapped.insert(srcSeq);
			}
			else // Node in sequence was already assigned
				return false;
		}
		else // node is mapped
		{
			if(!(mapSequence[srcM] == srcSeq))
				return false;
		}
		if(mapSequence[dstM] == -1) // node not mapped
		{
			if(alreadyMapped.find(dstSeq) == alreadyMapped.end())
			{
				mapSequence[dstM] = dstSeq;
				alreadyMapped.insert(dstSeq);
			}
			else // Node in sequence was already assigned
				return false;
		}
		else // node is mapped
		{
			if(!(mapSequence[dstM] == dstSeq))
				return false;
		}
		currentEdge++;
	}
	return true;
}

int TempMotifSamplerStatic::getLength(const long long number)
{
	return (std::floor(log10(number)) + 1);
}

/*
int TempMotifSamplerStatic::getDigit(const long long number,const int digit)
{
	return static_cast<int>( number / static_cast<int>(std::pow(10, digit - 1))) % 10;
}
*/

int TempMotifSamplerStatic::getDigit(const std::uint_fast64_t number,const int digit)
{
	return static_cast<int>( (number >> ((digit-1)*edgeMaskLength_m) ) & edgeMask_m);
}

void TempMotifSamplerStatic::DecrementCounts(const std::uint_fast64_t ID, std::unordered_map<std::uint_fast64_t, long long>& hashCounts, std::vector<std::unordered_set<std::uint_fast64_t>>& keysLength){

	//if(hashCounts.find(ID) == hashCounts.end()) std::cout << "ERRORE SYNCH DECREMENT ID!" << std::endl ; 
	hashCounts.at(ID) -= 1;
	if(hashCounts.at(ID) == 0) 
	{
		keysLength[0].erase(ID);
		//hashCounts.erase(ID);
	}
	for(int i{ 0 }; i < l_-2; i++)
	{
		for(const auto& key : keysLength[i])
		{
			//long long finalkey = CombineKeys(ID, key);
			std::uint_fast64_t finalkey = CombineKeys(ID, key, i+1);
			//if(hashCounts.find(finalkey) == hashCounts.end()) std::cout << "ERRORE SYNCH DECREMENT FINALKEY!" << std::endl; 
			//if(hashCounts.find(key) == hashCounts.end()) std::cout << "ERRORE SYNCH DECREMENT KEY!" << std::endl ; 
			hashCounts.at(finalkey) -= hashCounts[key];
			if(hashCounts.at(finalkey) == 0) 
			{
				keysLength[i+1].erase(finalkey);
				//hashCounts.erase(finalkey);
			}
		}
	}
	
	/*
	std::cout << "++++++++++++++++++++++++++++++++\n";
	std::cout << "DECREMENT: \n";
	for(auto& entry : hashCounts)
	{
		std::cout << "key: " << entry.first << " count: " << entry.second << '\n';
	}
	int i{0};
	for(auto list : keysLength)
	{
		std::cout << "i: " << i+1 << " list entries: ";
		for(auto& key : list)
		{
			std::cout << key << ", ";
		}
		i++;
		std::cout << '\n';
	}
	std::cout << "++++++++++++++++++++++++++++++++\n";
	*/
	
					
}


void TempMotifSamplerStatic::IncrementCounts(const std::uint_fast64_t ID, std::unordered_map<std::uint_fast64_t, long long>& hashCounts, std::vector<std::unordered_set<std::uint_fast64_t>>& keysLength){

	for(int i{ l_-2 }; i >= 0; i--)
	{
		for(const auto& key : keysLength[i])
		{
			std::uint_fast64_t finalkey = CombineKeys(key, ID, 1);
			long long currentCount{ hashCounts[key] };
			if(keysLength[i+1].find(finalkey) != keysLength[i+1].end()) // key exists
			{
				//if(hashCounts.find(finalkey) == hashCounts.end()) std::cout << "ERRORE SYNCH INCREMENT FINALKEY!" << std::endl ; 
				//if(hashCounts.find(key) == hashCounts.end()) std::cout << "ERRORE SYNCH INCREMENT KEY!" << std::endl ; 
				//hashCounts.at(finalkey) += hashCounts[key];
				hashCounts.at(finalkey) += currentCount;
			}
			else
			{
				//hashCounts[finalkey] = hashCounts[key];
				hashCounts[finalkey] = currentCount;
				keysLength[i+1].insert(finalkey);
			}
		}
	}
	if(keysLength[0].find(ID) != keysLength[0].end()) 
	{
		//if(hashCounts.find(ID) == hashCounts.end()) std::cout << "ERRORE SYNCH INCREMENT ID!" << std::endl ; 
		hashCounts.at(ID) += 1;// key exists
	}
	else
	{
		hashCounts[ID] = 1;	
		keysLength[0].insert(ID);
	}
	
	/*
	std::cout << "********************************\n";
	std::cout << "INCREMENT: \n";
	for(auto& entry : hashCounts)
	{
		std::cout << "key: " << entry.first << " count: " << entry.second << '\n';
	}
	int i{0};
	for(auto list : keysLength)
	{
		std::cout << "i: " << i+1 << " list entries: ";
		for(auto& key : list)
		{
			std::cout << key << ", ";
		}
		i++;
		std::cout << '\n';
	}
	std::cout << "********************************\n";
	*/
	

}

int TempMotifSamplerStatic::GetStaticEdgesInstance(std::vector<std::pair<node_t, node_t>>& instance)
{
	std::unordered_set<std::pair<node_t,node_t>, DataStructures::hash_pair> staticEdges;
	//int totalEdges{ 0 };
	for( const auto & edge : instance)
	{
		// Graph orientation -> transforming the instance in an undirected graph, considering only the underlying static structure
		std::pair<node_t, node_t> staticEdge = std::make_pair<node_t,node_t>(static_cast<node_t>(std::min(edge.first, edge.second)), static_cast<node_t>(std::max(edge.first, edge.second)) );
		staticEdges.emplace(staticEdge);
		/*
		if(staticEdges.find(staticEdge) == staticEdges.end())
		{
			staticEdges.insert(staticEdge);
			totalEdges++;
		}
		*/
	}
	//return totalEdges;
	return staticEdges.size();
}

std::uint_fast64_t TempMotifSamplerStatic::ComputeEncoding(std::vector<std::pair<node_t, node_t>>& instance)
{
	int ID { 1};
	std::unordered_map<node_t, node_t> mapping; 
	//Remapping nodes to be in 1,...,15 according to their first appeareance
	for(auto & edge : instance)
	{
		node_t src = edge.first;
		node_t dst = edge.second;

		if(mapping.find(src) == mapping.end())
		{
			edge.first = ID;
			mapping[src] = ID++;
		}
		else
			edge.first = mapping[src];

		if(mapping.find(dst) == mapping.end())
		{
			edge.second = ID;
			mapping[dst] = ID++;
		}
		else
			edge.second = mapping[dst];
	}
	//Computing encoding
	std::uint_fast64_t encode { static_cast<std::uint_fast64_t>(instance[0].first) };
	encode =( (encode << edgeMaskLength_m) | (static_cast<std::uint_fast64_t>(instance[0].second) & edgeMask_m));
	for(auto it{instance.begin()+1}; it <  instance.end(); it++)
	{
		encode = ( (encode << edgeMaskLength_m) | (static_cast<std::uint_fast64_t>(it->first) & edgeMask_m));
		encode = ( (encode << edgeMaskLength_m) | (static_cast<std::uint_fast64_t>(it->second) & edgeMask_m));
	}
	return encode;				
}

long long TempMotifSamplerStatic::CombineKeys(long long key1, long long key2)
{
	int toshift = getLength(key2);
	return (static_cast<long long>( key1 * pow(10.0, toshift)) + key2);
}

// lengthKey2: number of edges that Key2 represents
std::uint_fast64_t TempMotifSamplerStatic::CombineKeys(const std::uint_fast64_t key1, const std::uint_fast64_t key2, const int lengthKey2)
{
	std::uint_fast64_t combinedKey { key1 << (lengthKey2 * edgeMaskLength_m) };
	return (combinedKey | key2);
}


double TempMotifSamplerStatic::EnumerateOnStaticGraphWithLemon(int delta, int samp, int seed, double timeLimit, int ell, int triORsq)
{
	std::vector<std::pair<node_t, node_t>> motifSeq;
	lemon::ListGraph pattern = lemon::ListGraph();
	std::vector<lemon::ListGraph::Node> nodesPattern;
	for(int i{0}; i < static_undirected_motif_->GetNodes(); i++)
	{
		nodesPattern.push_back(pattern.addNode());
	}
	TUNGraph::TEdgeI edgeItP = static_undirected_motif_->BegEI();
	for(long long int i=0; i< static_undirected_motif_->GetEdges(); i++)
	{
		pattern.addEdge(nodesPattern[edgeItP.GetSrcNId()], nodesPattern[edgeItP.GetDstNId()]);
		motifSeq.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(edgeItP.GetSrcNId()), static_cast<node_t>(edgeItP.GetDstNId())));
		edgeItP++;
	}
		
	// Adding colors of the motifs in lemon
	lemon::ListGraph::NodeMap<int> mapColorsPattern(pattern);
	lemon::ListGraph::NodeMap<int> mapIDPattern(pattern);
	for(int i{0}; i< nodesPattern.size(); i++)
	{
		mapColorsPattern[nodesPattern[i]] = 1;
		mapIDPattern[nodesPattern[i]] = i;
	}
		
	// Adding nodes of the temporal network in lemon
	lemon::ListGraph target = lemon::ListGraph();
	std::vector<lemon::ListGraph::Node> nodesTarget;
	for(int i{0}; i < static_undirected_graph_->GetNodes(); i++)
	{
		nodesTarget.push_back(target.addNode());
	}
	TUNGraph::TEdgeI edgeItT = static_undirected_graph_->BegEI();
	for(long long int i=0; i< static_undirected_graph_->GetEdges(); i++)
	{
		target.addEdge(nodesTarget[edgeItT.GetSrcNId()], nodesTarget[edgeItT.GetDstNId()]);
		edgeItT++;
	}
		
	// Adding colors of the motifs in lemon
	lemon::ListGraph::NodeMap<int> mapColorsTarget(target);
	lemon::ListGraph::NodeMap<int> mapIDTarget(target);
	for(int i{0}; i< nodesTarget.size(); i++)
	{
		mapColorsTarget[nodesTarget[i]] = 1;
		mapIDTarget[nodesTarget[i]] = i;
	}

	// -- Computing Automorphisms of the motif 
	std::vector<int> equivalenceClasses;
	for(int i{0}; i < static_undirected_motif_->GetNodes(); i++)
	{
		equivalenceClasses.push_back(i);
	}
	
   	lemon::ListGraph::NodeMap<lemon::ListGraph::Node> mm(pattern);
	auto* myVf2ppo = vf2pp(pattern,pattern).mapping(mm).nodeLabels(mapColorsPattern,mapColorsPattern).getPtrToVf2ppObject();
	int countAuto = 0;
   	while(myVf2ppo->find()){
		//process the current mapping mm
		countAuto++;
		for(int ii{0}; ii < static_undirected_motif_->GetNodes(); ii++)
		{
			int nodeCurrMap = equivalenceClasses[ii];
			int nodeMap = equivalenceClasses[mapIDPattern[mm[nodesPattern[ii]]]];
			int newclass = std::min(nodeCurrMap, nodeMap);
			int oldclass = std::max(nodeCurrMap, nodeMap);
			if(newclass == nodeCurrMap)
				equivalenceClasses[mapIDPattern[mm[nodesPattern[ii]]]] = newclass;
			else
 				equivalenceClasses[ii] = newclass;
			// Quadratic cost in the number of edges of the pattern is not a problem since very fast on SMALL pattern graphs
			for(int ij{0}; ij < static_undirected_motif_->GetNodes(); ij++)
			{
				if(equivalenceClasses[ij] == oldclass)
					equivalenceClasses[ij] = newclass;
			}
		}
		
  	}
	delete myVf2ppo;

	std::unordered_set<std::pair<int, int>, DataStructures::hash_pair> edgesPatternColors;
	std::vector<std::pair<int, int>> edgesToFix;
	for(long long int i=0; i< static_cast<int>(motifSeq.size()); i++)
	{
		int equivSrc = equivalenceClasses[motifSeq[i].first];
		int equivDst = equivalenceClasses[motifSeq[i].second];
		auto currpair = std::make_pair<int, int>(static_cast<int>(std::min(equivSrc, equivDst)), static_cast<int>(std::max(equivSrc, equivDst)) );
		if(edgesPatternColors.find(currpair) == edgesPatternColors.end())
		{
			edgesPatternColors.emplace(currpair);
			edgesToFix.push_back(std::make_pair<int,int>(static_cast<int>(motifSeq[i].first), static_cast<int>(motifSeq[i].second)));
		}
	}

	l_ = ell;
	std::random_device rd;
	std::mt19937 generator(seed);
	
	std::vector<double> mass; 
	double maxdegst { 0};
	double avgdegst {0};
	double alphamax {0};
	double denom{0};
	TUNGraph::TEdgeI edgeIt = static_undirected_graph_->BegEI();
	std::vector<std::pair<node_t, node_t>> edges_exp;
	TVec< THash<TInt, TInt64 >> probsEdges = TVec< THash<TInt, TInt64 >>(static_undirected_graph_->GetNodes());
	double firsttime = get_wall_time();
	for(long long int i=0; i< static_undirected_graph_->GetEdges(); i++)
	{
		node_t src = edgeIt.GetSrcNId();
		node_t dst = edgeIt.GetDstNId();

		(src < dst) ? edges_exp.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst))) : edges_exp.push_back(std::make_pair<node_t, node_t>(static_cast<node_t>(dst), static_cast<node_t>(src)));
		edgeIt++;
		double currmass { 0 };
		if(temporal_data_[src].IsKey(dst))
		{
			currmass += temporal_data_[src](dst).Len();
			denom++;
		}
		double mval = currmass;
		if(temporal_data_[dst].IsKey(src))
		{
			currmass += temporal_data_[dst](src).Len();
			denom++;
		}

		mass.push_back(static_cast<double>(currmass));
		(src < dst) ? probsEdges[src].AddDat(dst, i) : probsEdges[dst].AddDat(src, i);
	}

	std::discrete_distribution<long long int> importance(mass.begin(), mass.end());
	std::vector<double> sampling_probs{ importance.probabilities() };
	double randomtime = get_wall_time() - firsttime;

	double totestimate{ 0};
	int s{ samp };
	double timeEnumerateTri{ 0 };
	double timeBuildSample{ 0 };
	double timeEnumerateTemporal{ 0 };
	double timenow{ 0 };

	double timeSequence{ 0 };
	double timeIso{ 0 };
	double totInstances{ 0 };
	double totStaticInstances{ 0 };
	std::vector<int> numEdges;
	std::vector<int> staticInsts;
	int pruned{ 0 };
	std::unordered_map<std::uint_fast64_t, double> counts;
	double artifactNorm { 0 };
	//std::unordered_map<int, int> lvals;
	long long samples1 = 0;
	double init = get_wall_time();
	while(samples1 < s)
	//while( (get_wall_time() - init) < timeLimit)
	{
			samples1++;
			long long int samplednode { importance(generator) };
			
			std::vector<std::vector<std::pair<node_t, node_t> > > static_graphlets;
			
			double initVF2pp = get_wall_time();

			double timemine = get_wall_time();

			int countvf2 = 0;

			mapColorsTarget[nodesTarget[edges_exp[samplednode].first]] = 2;
			mapColorsTarget[nodesTarget[edges_exp[samplednode].second]] = 2;
			node_t srcSampled = std::min(edges_exp[samplednode].first, edges_exp[samplednode].second) ;
			node_t dstSampled = std::max(edges_exp[samplednode].first, edges_exp[samplednode].second) ;
			std::vector<std::vector<std::pair<node_t, node_t> > > static_graphletsVF2pp;

			for(auto toFix : edgesToFix){
			mapColorsPattern[nodesPattern[toFix.first]] = 2;
			mapColorsPattern[nodesPattern[toFix.second]] = 2;
    		lemon::ListGraph::NodeMap<lemon::ListGraph::Node> m(pattern);

			std::unordered_set<std::string> foundEmbeddings;
			auto* myVf2pp = vf2pp(pattern,target).mapping(m).nodeLabels(mapColorsPattern,mapColorsTarget).getPtrToVf2ppObject();
    		while(myVf2pp->find()){
				//process the current mapping m
				countvf2++;

				std::vector<std::pair<node_t, node_t>> embedding;
				std::string encode = "";
				//Get the current embedding without the sampled edge
				for(long long int ix=0; ix< static_cast<long long>(motifSeq.size()); ix++)
				{
					node_t src = mapIDTarget[m[nodesPattern[motifSeq[ix].first]]];
					node_t dst = mapIDTarget[m[nodesPattern[motifSeq[ix].second]]];
					if(src > dst)
					{
						node_t tmp = src;
						src = dst;
						dst = tmp;
					}
					//encode = encode + std::to_string(src) + std::to_string(dst);
					auto currpair = std::make_pair<node_t, node_t>(static_cast<node_t>(src), static_cast<node_t>(dst));
					
					if((src != srcSampled) || (dst != dstSampled))
						embedding.push_back(currpair);
				}
				// Sort for the encoding to be unique (i.e., permutation invariant)
				std::sort(embedding.begin(), embedding.end(), [](auto &left, auto&right){
					if(left.first == right.first)
						return left.second < right.second;
					return left.first < right.first;
				});
				//Compute encoding
				for(auto e : embedding)
				{
					encode = encode + std::to_string(e.first) + std::to_string(e.second);
				}
						
				//If not already seen this map under a different permutation then save and remember it. 
				if(foundEmbeddings.find(encode) == foundEmbeddings.end())
				{
					static_graphletsVF2pp.push_back(embedding);
					foundEmbeddings.emplace(encode);
				}
    		}
			delete myVf2pp;
			mapColorsPattern[nodesPattern[toFix.first]] = 1;
			mapColorsPattern[nodesPattern[toFix.second]] = 1;
			}

			mapColorsTarget[nodesTarget[edges_exp[samplednode].first]] = 1;
			mapColorsTarget[nodesTarget[edges_exp[samplednode].second]] = 1;
			double endVF2pp = get_wall_time();

			for(auto emb : static_graphletsVF2pp)
				static_graphlets.push_back(emb);

			timeEnumerateTri += (get_wall_time()-timemine);
			if(static_graphlets.size() == 0)
				continue;

			timemine = get_wall_time();
			std::vector<std::pair<node_t, node_t>> staticEdgesInSample;
			double sampled{ 0 };
			
			totStaticInstances += static_cast<double>(static_graphlets.size());
			std::vector<std::pair<timestamp_t,size_t>> times; //vector of <timestamp, idx_static_edge_In_sample (which corresponds to the timestamp)> 

			//Computing on sampled edge to avoid this repeated computation
			std::vector<std::pair<timestamp_t,size_t>> timesSampledEdge;
			node_t sampledSrc = static_cast<node_t>(edges_exp[samplednode].first);
			node_t sampledDst = static_cast<node_t>(edges_exp[samplednode].second);
			std::pair<node_t, node_t> P1{ std::make_pair(sampledSrc, sampledDst) };
			std::pair<node_t, node_t> P2{ std::make_pair(sampledDst, sampledSrc) };
			bool exsistenceSD {false};
			bool exsistenceDS {false};
			if(temporal_data_[sampledSrc].IsKey(sampledDst))
			{
				exsistenceSD = true;
				const TVec<TInt64>& timestamps = temporal_data_[sampledSrc].GetDat(sampledDst);
				for(size_t j=0; j< timestamps.Len(); j++)
				{
					timesSampledEdge.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], 0));
				}
			}

			if(temporal_data_[sampledDst].IsKey(sampledSrc))
			{
				exsistenceDS = true;
				const TVec<TInt64>& timestamps = temporal_data_[sampledDst].GetDat(sampledSrc);
				for(size_t j=0; j< timestamps.Len(); j++)
				{
					if(exsistenceSD)
						timesSampledEdge.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], 1));
					else
						timesSampledEdge.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], 0));
				}
			}
			timeBuildSample += (get_wall_time() - timemine);


			for(size_t k{ 0 }; k < static_graphlets.size(); k++)
			{
				//timenow = get_wall_time();
				timemine = get_wall_time();
				staticEdgesInSample.clear();
				times.clear();

				if(exsistenceSD)
					staticEdgesInSample.push_back(P1);
				if(exsistenceDS)
					staticEdgesInSample.push_back(P2);

				times.insert(times.begin(), timesSampledEdge.begin(), timesSampledEdge.end());
				



				//Collecting graphlet edges
				long id{0};
					
				if(triORsq != 3)
				{
					for(const auto& pair : static_graphlets[k])
					{
						node_t src = pair.first;
						node_t dst = pair.second;
						if(temporal_data_[src].IsKey(dst))
						{
							staticEdgesInSample.push_back(std::make_pair<node_t,node_t>(static_cast<node_t>(src), static_cast<node_t>(dst)));
							size_t currentStaticEdges{ staticEdgesInSample.size() };
							const TVec<TInt64>& timestamps = temporal_data_[src].GetDat(dst);
							for(size_t j=0; j< timestamps.Len(); j++)
							{
								times.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], currentStaticEdges-1));
							}
						}

						if(temporal_data_[dst].IsKey(src))
						{
							staticEdgesInSample.push_back(std::make_pair<node_t,node_t>(static_cast<node_t>(dst), static_cast<node_t>(src)));
							size_t currentStaticEdges{ staticEdgesInSample.size() };
							const TVec<TInt64>& timestamps = temporal_data_[dst].GetDat(src);
							for(size_t j=0; j< timestamps.Len(); j++)
							{
								times.push_back(std::make_pair<timestamp_t, size_t>(timestamps[j], currentStaticEdges-1));
							}
						}
					}
				}

					std::sort(times.begin(), times.end());
					timeBuildSample += (get_wall_time() - timemine);
					
					timemine = get_wall_time();
					bool foundCandidate{ false };

					// Pruning
					for(long long int t{ 0 }; t < static_cast<long long int>(times.size())-l_+1; t++)
					{
						if( ((times[t+l_-1].first - times[t].first) <= delta ))
						{
								foundCandidate = true;
								break;
						}
					}
					//timeBuildSample += (get_wall_time() - timenow);
					timeEnumerateTemporal += (get_wall_time() - timemine);
					if(!foundCandidate)
					{
						pruned++;
					}
					else
					{
	
							//Computing probability and checking if the selected vertex is in the graphlet
							timemine = get_wall_time();
							double graphlet_prob{0};
							graphlet_prob += sampling_probs[samplednode];
							

							double weight { static_cast<double>(1./graphlet_prob) }; 

							double ressamp  = FastCountSequence(delta, std::ref(times), std::ref(staticEdgesInSample), timeSequence, timeIso, std::ref(counts), graphlet_prob);
							totInstances += ressamp;
							numEdges.push_back(static_cast<int>(times.size()));

							ressamp = ressamp * weight;
	
							sampled++;

							totestimate += ressamp;
							timeEnumerateTemporal += (get_wall_time() - timemine);
					}

			} // TODO: Delete tmp pointers!!
		
	}

	double tottime = get_wall_time() - firsttime;
	std::cout.precision(2);
	for(auto& el : counts)
	{
		std::cout << " count is: " <<  std::fixed << el.second/(1.*samples1*static_undirected_motif_->GetEdges()) << " motif is: " << std::hex << el.first << std::dec;
		std::cout << '\n';
	}
	std::cout << "T0 time spent in building prob distrib: " << std::fixed << randomtime  << '\n';
	std::cout << "T1 time spent in static enumeration: " << std::fixed << timeEnumerateTri  << '\n';
	std::cout << "T2 time spent in building sample: " << std::fixed << timeBuildSample  << '\n';
	std::cout << "T3 time spent in temporal enumeration: " << std::fixed << timeEnumerateTemporal  << '\n';
	std::cout << "T3-a time spent in temporal enumeration when enumerating sequences: " << std::fixed << timeSequence  << '\n';
	std::cout << "T3-b time spent in temporal enumeration with isomorphism: " << std::fixed << timeIso  << '\n';
	std::cout << "T4 time spent in whole precedure: " << std::fixed << tottime  << '\n';

	std::cout << "samples s examined: " << samples1 << '\n';
	return totestimate/s;
	
}
