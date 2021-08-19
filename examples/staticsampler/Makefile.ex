#
#	configuration variables for the example

## Main application file
MAIN=staticsamplermain
DEPH=$(EXSNAPADV)/temporalmotifstaticsampler.h
DEPCPP=$(EXSNAPADV)/temporalmotifstaticsampler.cpp $(EXSNAPADV)/datastructures.cpp $(EXSNAPADV)/utilities.cpp $(EXSNAPADV)/temporalmotifs.cpp
# $(EXGLASGLIB)/formats/csv.cc $(EXGLASGLIB)/formats/dimacs.cc $(EXGLASGLIB)/formats/graph_file_error.cc $(EXGLASGLIB)/formats/input_graph.cc $(EXGLASGLIB)/formats/lad.cc $(EXGLASGLIB)/formats/read_file_format.cc $(EXGLASGLIB)/formats/vfmcs.cc $(EXGLASGLIB)/cheap_all_different.cc $(EXGLASGLIB)/clique.cc $(EXGLASGLIB)/common_subgraph.cc $(EXGLASGLIB)/configuration.cc $(EXGLASGLIB)/graph_traits.cc $(EXGLASGLIB)/homomorphism.cc $(EXGLASGLIB)/homomorphism_domain.cc $(EXGLASGLIB)/homomorphism_model.cc $(EXGLASGLIB)/homomorphism_searcher.cc $(EXGLASGLIB)/homomorphism_traits.cc $(EXGLASGLIB)/lackey.cc $(EXGLASGLIB)/proof.cc $(EXGLASGLIB)/restarts.cc $(EXGLASGLIB)/sip_decomposer.cc $(EXGLASGLIB)/svo_bitset.cc $(EXGLASGLIB)/symmetries.cc $(EXGLASGLIB)/thread_utils.cc $(EXGLASGLIB)/timeout.cc $(EXGLASGLIB)/verify.cc $(EXGLASGLIB)/watches.cc 
