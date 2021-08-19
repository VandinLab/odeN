#ifndef snap_temporalmotifstaticsampler_h
#define snap_temporalmotifstaticsampler_h


#include "Snap.h"

#include "temporalmotiftypes.h"
#include "temporalmotifs.h"

#include "datastructures.h"

#include "utilities.h"

#include <vector>
#include <map>
#include <stack>
#include <queue>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <thread>
#include <future>
#include <fstream>
#include <sstream>
#include <iterator>
#include <chrono>
#include <ctime>
#include <cstdint>
#include <memory>


using std::chrono::operator""s;

#define VERBOSE 100
#define THREADS 1 
#define SINGLETHREAD 1 


using DataStructures::TEdge;
using DataStructures::Datum;
using DataStructures::Partition;

template<typename R>
  bool is_ready(std::future<R> const& f)
  { return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready; }


class TempMotifSamplerStatic {
public:
	double get_wall_time();
	

    TempMotifSamplerStatic(const TStr& filenameG, const TStr& filenameM);
	TempMotifSamplerStatic(std::string filenameG, std::string filenameM, bool parallel=false);
	
	// USING STATIC INFORMATION
	double EnumerateOnStaticGraph(int delta, int samp, int seed, double timeLimit, int ell, int triORsq);

	long long CombineKeys(const long long key1, const long long key2);

	int getLength(const long long number);

	//double FastCountSequence(int delta, std::vector<TEdge>& tmp, std::vector<std::pair<node_t, node_t>>& staticEdgesInSample, double& timeSequence, double& timeIso); 

	double FastCountSequence(int delta, std::vector<std::pair<timestamp_t, size_t>>& tmp, std::vector<std::pair<node_t, node_t>>& staticEdgesInSample, double& timeSequence, double& timeIso, std::unordered_map<std::uint_fast64_t, double>& counts, double instanceProbability);

	//void DecrementCounts(const int ID, std::unordered_map<long long, long long>& hashCounts, std::vector<std::unordered_set<long long>>& keysLength);
	void DecrementCounts(const std::uint_fast64_t ID, std::unordered_map<std::uint_fast64_t, long long>& hashCounts, std::vector<std::unordered_set<std::uint_fast64_t>>& keysLength);

	//void IncrementCounts(const int ID, std::unordered_map<long long, long long>& hashCounts, std::vector<std::unordered_set<long long>>& keysLength);
	void IncrementCounts(const std::uint_fast64_t ID, std::unordered_map<std::uint_fast64_t, long long>& hashCounts, std::vector<std::unordered_set<std::uint_fast64_t>>& keysLength);

	int getDigit(long long number, int digit);

	bool testIsomorphism(std::vector<std::pair<node_t, node_t>>& nodeSequence, std::vector<std::pair<int, int>>& motif);

	std::uint_fast64_t CombineKeys(const std::uint_fast64_t key1, const std::uint_fast64_t key2, const int lengthKey2);

	int getDigit(const std::uint_fast64_t number,const int digit);
	
	int GetStaticEdgesInstance(std::vector<std::pair<node_t, node_t>>& instance);

	std::uint_fast64_t ComputeEncoding(std::vector<std::pair<node_t, node_t>>& instance);

	int setMotif(std::string filenameM);

	double EnumerateOnStaticGraphParallel(int delta, int samp, int seed, double timeLimit, int ell, int triORsq, int nthreads);

	double EnumerateOnStaticGraphWithLemon(int delta, int samp, int seed, double timeLimit, int ell, int triORsq);
	
private:
	double get_wall_time_inner(){

		struct timeval time;
		if (gettimeofday(&time,NULL)){
			//  Handle error
			return 0;
		}
		return (double)time.tv_sec + (double)time.tv_usec * .000001;
	}
	std::vector< std::vector<TEdge> > all_components;
	std::vector< TEdge > edges_;
	std::vector< std::pair<int, int> > edgesM_;
    std::default_random_engine generatorLiu;
	int staticEdgesMotif_m {0};
	int Vm_;
	int l_ = 0;
	int ems_; // edges in the directed static motif
	long long int m_;
	Permutation permut_;
	std::vector<std::pair<std::pair<node_t, node_t>, std::vector<timestamp_t>>> staticandtemporalgraph_;
	PNGraph static_directed_graph_;
	PNGraph static_motif_;
	PUNGraph static_undirected_graph_;
	PUNGraph static_undirected_motif_;
	TVec< THash<TInt, TVec<TInt64> > > temporal_data_;
	std::vector<size_t> temporal_degrees_;
	const std::uint_fast64_t edgeMask_m{ 0x0F };
	const int edgeMaskLength_m{ 4 };
	timestamp_t timespan_m;
};

#endif  // snap_temporalmotifstaticsampler_h
