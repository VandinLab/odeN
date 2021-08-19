#include "temporalmotifstaticsampler.h"

#include <omp.h>

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mutex>
#include <ctime>

double get_wall_time(){

    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int main(int argc, char* argv[]) {
	const rlim_t kStackSize = RLIM_INFINITY;   // min stack size = 16 MB
    struct rlimit rl;
    int status;

	status = getrlimit(RLIMIT_STACK, &rl);
    if (status == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            status = setrlimit(RLIMIT_STACK, &rl);
            if (status != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", status);
            }
        }
    }
	

  Env = TEnv(argc, argv, TNotify::StdNotify);

  const TStr temporal_graph_filename =
    Env.GetIfArgPrefixStr("-i:", "simple-example.txt",
			  "Input directed temporal graph file");
  const TStr motif_graph_filename =
    Env.GetIfArgPrefixStr("-m:", "simple-motif.txt",
        "Input directed motif graph file");
  const TFlt delta =
    Env.GetIfArgPrefixFlt("-delta:", 4096, "Time window delta");
  const int num_threads =
    Env.GetIfArgPrefixInt("-nt:", 4, "Number of threads for parallelization");
  const TInt type =
    Env.GetIfArgPrefixFlt("-type:", 0, "Which alg?");
  const int samp=
    Env.GetIfArgPrefixFlt("-samp:", 100, "Time window delta");
  const TFlt timeLimit =
    Env.GetIfArgPrefixFlt("-tl:", 10, "Time limit in seconds for running the sampling algorithm ");
  const TInt execs =
    Env.GetIfArgPrefixFlt("-ex:", 1, "How may times to exec the procedure? integer between [1, 10] !");
  const TInt ell =
    Env.GetIfArgPrefixFlt("-ell:", 3, "How many edges the temporal motif should have, k-node ell-edge!");
  const TInt whichmotif =
    Env.GetIfArgPrefixFlt("-topo:", 1, "Which motif? 0 wedge, 1 triangles, 2 square");
  
  
  omp_set_num_threads(8);
  double st = get_wall_time(); 
  std::string file_name{temporal_graph_filename.CStr()};
 	std::string motif_name{motif_graph_filename.CStr()};
  TempMotifSamplerStatic tmc(temporal_graph_filename, motif_graph_filename);

  Env.PrepArgs(TStr::Fmt("Temporalmotifs. build: %s, %s. Time: %s",
       __TIME__, __DATE__, TExeTm::GetCurTm()));  
	printf("Time to load the dataset: %lfs\n\n", get_wall_time() - st);
    	srand(time(nullptr));
  
  double init;
  double end;

  double frac = 0;
  int samples = 0;
  double a=0;
	
	
	
	if(type == 1)
	{
		int seeds[10] {10, 254, 4321, 91283, 39438, 10239102, 11864, 65638, 891237, 7684};
		int alphas[10] {1,2,3,4,5,6,7, 8, 9, 10};
		for(int i{ 0 }; i < execs; i++)
		{	
			init = get_wall_time();
			std::clock_t c_start = std::clock();
			double myresult {tmc.EnumerateOnStaticGraph(delta, samp, seeds[i], timeLimit, ell, whichmotif)};
			std::clock_t c_end = std::clock();
			std::cout.precision(2);
			std::cout << "--SSTATIC-- static graphlet count is: " << std::fixed << myresult << '\n';
			std::cout << "--STATICTIME-- Time to enumerate static instances: " << get_wall_time()-init << '\n';
			std::cout << std::fixed << std::setprecision(2) << "CPU Time used: " << 1.* (c_end-c_start) / CLOCKS_PER_SEC << "s\n";
		}
	}

	if(type == 2)
	{
		int seeds[10] {10, 254, 4321, 91283, 39438, 10239102, 11864, 65638, 891237, 7684};
		for(int i{ 0 }; i < execs; i++)
		{	
			init = get_wall_time();
			std::clock_t c_start = std::clock();
			double myresult {tmc.EnumerateOnStaticGraphParallel(delta, samp, seeds[i], timeLimit, ell, whichmotif, num_threads)};
			std::clock_t c_end = std::clock();
			std::cout.precision(2);
			std::cout << "--SSTATIC-- static graphlet count is: " << std::fixed << myresult << '\n';
			std::cout << "--STATICTIME-- Time to enumerate static instances: " << get_wall_time()-init << '\n';
			std::cout << std::fixed << std::setprecision(2) << "CPU Time used: " << 1.* (c_end-c_start) / CLOCKS_PER_SEC << "s\n";
		}
	}
	else if(type == 3)
	{
		int seeds[10] {10, 254, 4321, 91283, 39438, 10239102, 11864, 65638, 891237, 7684};
		int alphas[10] {1,2,3,4,5,6,7, 8, 9, 10};
		for(int i{ 0 }; i < execs; i++)
		{	
			init = get_wall_time();
			std::clock_t c_start = std::clock();
			double myresult {tmc.EnumerateOnStaticGraphWithLemon(delta, samp, seeds[i], timeLimit, ell, whichmotif)};
			std::clock_t c_end = std::clock();
			std::cout.precision(2);
			std::cout << "--SSTATICLEMON-- static graphlet count is: " << std::fixed << myresult << '\n';
			std::cout << "--STATICTIMELEMON-- Time to enumerate static instances: " << get_wall_time()-init << '\n';
			std::cout << std::fixed << std::setprecision(2) << "CPU Time used: " << 1.* (c_end-c_start) / CLOCKS_PER_SEC << "s\n";
		}

	}


  return 0;
}
