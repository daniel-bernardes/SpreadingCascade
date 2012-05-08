/*
  SIMPLE EPIDEMIC CASCADE SIMULATION:
  SIR process such that infected nodes become recovered in one time step
  Output: the complete trace of the spreading -- ie, including the spread
  attempts to removed individuals.

  Daniel.Bernardes@lip6.fr, 2011
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <dirent.h>
#include <omp.h>

#include "prelim.c"
#include "queue.c"

// misc defs and utils
#define VERBOSE 1
#define PARALLEL 1
#define MAX_PATH_LENGTH 4096

inline char *tstamp() {
  time_t now = time(NULL);
  char *str = asctime(localtime(&now));
  str[strlen(str)-1]=' ';
  return str;
}

inline void techo(char *str) {
#if VERBOSE > 1
  #if PARALLEL
    fprintf(stdout, "%s -- Thread %d -- %s\n",
	    tstamp(), omp_get_thread_num(), str);
  #else
    fprintf(stdout, "%s -- %s\n", tstamp(), str);
  #endif
  fflush(stdout);
#endif
}

// Epidemic management
typedef enum _Stop_criterion {MaxTime, NumInfected} Stopc;
const char *stopc_description[] = {"maxdepth","maxsize"};

typedef struct _InitialCondition {
  int id;                 // epidemic id
  int num_infected;       // number of infected nodes
  int *infected;          // list of infected nodes' id
  int bound;              // bounds on epidemic evolution in terms of a ...
  Stopc stop_criterion;   // ... e.g., max time or max num infected)
} InitialCondition;

typedef struct _Epidemic {
  int id;                 // epidemic id
  int t;                  // time steps elapsed
  int num_infected;       // number of currently infected nodes
  int cascade_links;      // number of arcs in the infection cascade
  int bound;              // bounds on epidemic evolution in terms of ...
  Stopc stop_criterion;   // ... e.g., max time or max num infected
  double p;               // neighbor infection probability
  graph *g;               // underlying graph (network)
  FILE *output;           // trace output
  int *infected;          // set of all infected nodes
  Queue *active;          // list of active infected nodes
} Epidemic;

Epidemic *epidemic_new(double p, graph *g, InitialCondition *ic, FILE *output) {
  int i;
  Epidemic *epidemic = (Epidemic *) malloc(sizeof(Epidemic));
  assert(epidemic != NULL);
  assert(ic != NULL);
  epidemic->id             = ic->id;
  epidemic->t              = 1;
  epidemic->num_infected   = ic->num_infected;
  epidemic->cascade_links  = 0;
  epidemic->bound          = ic->bound;
  epidemic->stop_criterion = ic->stop_criterion;
  epidemic->p              = p;
  epidemic->g              = g;
  epidemic->output         = output;
  epidemic->active         = queue_new(g->n);
  epidemic->infected       = (int *) calloc(g->n, sizeof(int));
  assert(epidemic->infected != NULL);
  for (i = 0; i < ic->num_infected; i++) {
    queue_add(epidemic->active, ic->infected[i]);
    epidemic->infected[ic->infected[i]] = 1; // the initial time;
  }
  return epidemic;
}

void epidemic_destroy(Epidemic *epidemic) {
  assert(epidemic != NULL);
  epidemic->g = NULL; // don't destroy the graph, since it's shared a structure generally
  free(epidemic->infected);
  queue_destroy(epidemic->active);
  free(epidemic);
  epidemic = NULL;
}

/**
   Run epidemic spreading until the bound condition (on time or size) is met
 */
void epidemic_run(Epidemic *epidemic) {
  int i, u, v, t;
  
  while (!queue_empty(epidemic->active)) {
    u = queue_get(epidemic->active); // provider
    t = epidemic->infected[u];       // current time
    if (epidemic->stop_criterion == MaxTime && epidemic->bound < t)
      return;
    for (i = 0; i < epidemic->g->degrees[u]; i++) {
      v = epidemic->g->links[u][i];  // client
      if ( (double)rand() <= (double)RAND_MAX * epidemic->p ) {
	if ( !epidemic->infected[v] ) {
	  epidemic->infected[v] = t+1;
	  queue_add(epidemic->active, v);
	  epidemic->num_infected++;
	  epidemic->cascade_links++;
	  epidemic->t = t;
	  if (epidemic->stop_criterion == NumInfected && epidemic->bound == epidemic->num_infected) {
	    if (epidemic->output) // print output: t P C F
	      fprintf(epidemic->output, "%d %d %d %d\n", t, u, v, epidemic->id);
	    return;
	  }
	} else if (epidemic->infected[v] == t+1)
	  epidemic->cascade_links++;
	if (epidemic->output) // print output: t P C F
	  fprintf(epidemic->output, "%d %d %d %d\n", t, u, v, epidemic->id);
      }
    }
  }
}

/**
   Allocates a set of n infected nodes' id
*/
inline void ic_init(InitialCondition *ic, int n) {
  ic->num_infected = n;
  ic->infected = (int *) calloc(n, sizeof(int));
  assert(ic->infected != NULL);
}

/**
   De-allocates a set of n infected nodes' id
*/
inline void ic_clean(InitialCondition *ic) {
  if(ic) {
    free(ic->infected);
    ic->infected = NULL;
    ic->num_infected = 0;
  }
}

/**
   Returns the address of a new initial condition with one infected node (id = 0)
*/
inline InitialCondition *ic_trivial() {
  InitialCondition *ic = (InitialCondition *) calloc(1,sizeof(InitialCondition));
  assert(ic != NULL);
  ic_init(ic, 1);
  ic->infected[0] = 0;
  return ic;
}

/**
   Picks ic->num_infected distinct infected nodes from 0, ..., total_nodes
   and stores their ids into ic->infected.
*/
void ic_infect_randomly(InitialCondition *ic, int total_nodes) {
  int num_infected = ic->num_infected;
  int i, v, k, *infected = (int *) calloc(total_nodes, sizeof(int));
  assert(infected != NULL);
  assert(num_infected < total_nodes);

  // if num_infected > total_nodes/2, pick the non-infected nodes
  k = (num_infected <= total_nodes/2) ? num_infected : total_nodes-num_infected;
  for (i = 0; i < k; i++) {
    do
      v = rand() % total_nodes;
    while (infected[v]);
    infected[v] = !infected[v];
  }

  i = 0;
  for (v = 0; v < total_nodes; v++)
    if ((num_infected < total_nodes/2 &&  infected[v]) ||
	(num_infected > total_nodes/2 && !infected[v]))
      ic->infected[i++] = v;

  free(infected);
}

/**
   Import initial conditions from file into the array *ic. If 'total_nodes'
   is non-zero, pick infected nodes' id randomly from 0, ..., 'total_nodes',
   otherwise read the nodes' id from the corresponding line
   File format:
   <number of epidemics>
   <epidemic id> <N, number of infected nodes> [<node 1> ... <node N>]
   ...
*/
int ic_import(InitialCondition **ic, FILE *input, int total_nodes) {
  int i, j, id, num_infected, tokens_read, epidemics = 0;
  assert(input != NULL);
  tokens_read = fscanf(input, "%d\n", &epidemics);
  assert(tokens_read == 1);
  assert(epidemics > 0);
  assert(ic != NULL);
  *ic = (InitialCondition *) calloc(epidemics, sizeof(InitialCondition));
  assert(*ic != NULL);

  for (i = 0; i < epidemics; i++) {
    tokens_read = fscanf(input, "%d %d", &id, &num_infected);
    assert(tokens_read == 2 && num_infected > 0);
    ic_init(*ic+i, num_infected);
    (*ic+i)->id = id;
    if(total_nodes)
      ic_infect_randomly(*ic+i, total_nodes);
    else
      for (j = 0; j < num_infected; j++) {
	tokens_read = fscanf(input, "%d", &(*ic+i)->infected[j]);
	assert(tokens_read == 1);
      }
  }
  return epidemics;
}

/**
   Import stop bounds for each epidemic in the array *ic from file
   composed of a collection of lines with: <id> <bound>
*/
void ic_import_bounds(InitialCondition *ic, int n, Stopc stop_criterion, FILE *input) {
  int i, id, bound, tokens_read;
  assert(n > 0);
  assert(ic != NULL);
  assert(input != NULL);
  
  for (i = 0; i < n; i++) {
    tokens_read = fscanf(input, "%d %d\n", &id, &bound);
    assert(tokens_read == 2);
    assert(id == ic[i].id);
    ic[i].bound = bound;
    ic[i].stop_criterion = stop_criterion;
  }
}

/**
   Main
*/
int main(int argc, char **argv) {
  int i, j, epidemics, tid = 0;
  char epidemic_output_path[MAX_PATH_LENGTH] = "";
  FILE *graph_input, *ic_list_input, *bounds_list_input, 	\
    *data_output = NULL, *epidemic_output = NULL;
  graph *g;
  InitialCondition *ic;
  Epidemic *epidemic;
  Stopc stop_criterion;

  // default parameters
  double p               = 0;    // neighbor infection probability
  int maxtime            = 0;    // global maximum epidemic simulation time
  int sample_epidemics   = 1;    // number of sample epidemics
  int threads            = 1;    // number of threads
  char *graph_path       = NULL; // input path for graph (network) file
  char *ic_list_path     = NULL; // input path for list of epidemic initial parameters
  char *bounds_list_path = NULL; // input path for list of epidemic bounds
  char *trace_output_path= NULL; // output path for trace

  // parameter parsing
  char syntax[] = "\n Required parameters:\n\t -p SPREADING_PROBABILITY\n\t -g GRAPH_PATH\n \
Required parameters (one choice among the options):\n\t -t GLOBAL_MAX_TIME\n\t -a MAX_TIME_LIST_PATH\n\t -b MAX_INFECTED_LIST_PATH\n \
Optional parameters:\n\t -s NUM_SAMPLE_EPIDEMICS\n\t -i INITIAL_CONDITIONS_DATA_PATH \n\t -h NUM_THREADS\n \t -e [STATUS_OUTPUT_PATH]\n\t -o SPREADING_OUTPUT\n\n";
  fprintf(stderr, "SIMPLE EPIDEMIC CASCADE SIMULATION:\n\n");
  while ((i = getopt(argc, argv, "e::o:p:s:g:i:t:a:b:h:")) != -1)
    switch (i) {
    case 'p':
      p = atof(optarg);
      break;
    case 'e':
      if (optarg)
	data_output = fopen(optarg, "w");
      else
	data_output = stdout;
      break;
    case 'o':
      trace_output_path = optarg;
      break;
    case 's':
      sample_epidemics = atoi(optarg);
      break;
    case 'g':
      graph_path = optarg;
      break;
    case 'i':
      ic_list_path = optarg;
      break;
    case 't':
      assert(bounds_list_path == NULL);
      maxtime = atoi(optarg);
      stop_criterion = MaxTime;
      break;
    case 'a':
    case 'b':
      assert(maxtime <= 0);
      bounds_list_path = optarg;
      stop_criterion = i - 'a';
      break;
    case 'h':
      threads = atoi(optarg);
      break;
    case '?':
      fputs(syntax, stderr);
    default:
      abort();
    }
  assert(p > 0.0 && p <= 1.0);
  assert(sample_epidemics > 0);
  assert(graph_path || ic_list_path);
  assert(bounds_list_path || maxtime > 0);
  assert(threads > 0);

  // preliminaires
  srand((unsigned)time(NULL));
  #if PARALLEL
    omp_set_num_threads(threads);
  #else
    threads = 1;
  #endif
  fprintf(stderr,"Number of threads: %d %s %s\n\n", threads,
	  !trace_output_path? "" : ", with trace output", !trace_output_path? "" : trace_output_path);
  fflush(stderr);

  // set list of initial conditions
  fprintf(stderr,"%s\n Loading list of epidemics %s...\n", tstamp(), ic_list_path? ic_list_path : "");
  fflush(stderr);
  if(!ic_list_path) {
    fprintf(stderr,"  No list of initial conditions given; using 1 epidemic with 1 infected node...\n");
    epidemics = 1;
    ic = ic_trivial();
  } else {
    ic_list_input = fopen(ic_list_path, "r");
    epidemics = ic_import(&ic, ic_list_input, 0);
    fclose(ic_list_input);
  }

  // set bounds and stop criterion
  fprintf(stderr,"Setting bounds (%s) for epidemics...\n", 
	  bounds_list_path? bounds_list_path : ":global:");
  fflush(stderr);
  if(maxtime)
    for(i = 0; i < epidemics; i++) {
      ic[i].bound = maxtime;
      ic[i].stop_criterion = MaxTime;
    }
  else {
    bounds_list_input = fopen(bounds_list_path, "r");
    ic_import_bounds(ic, epidemics, stop_criterion, bounds_list_input);
    fclose(bounds_list_input);
  }
  fprintf(stderr,"  Loaded %d epidemics.\n\n", epidemics);
  fflush(stderr);

  // load underlying graph
  fprintf(stderr,"%s\nLoading the graph %s...\n", tstamp(), graph_path? graph_path : "");
  fflush(stderr);
  if (!graph_path)
    graph_input = stdin;
  else 
    graph_input = fopen(graph_path, "r");
  assert(graph_input != NULL);
  g = graph_from_file(graph_input);
  if (graph_input != stdin)
    fclose(graph_input);
  fprintf(stderr,"  Loaded graph with %d nodes, %d links.\n\n", g->n, g->m);
  fflush(stderr);

  // set global epidemic_output
  assert(sample_epidemics == 1);
  if (trace_output_path && strlen(trace_output_path) > 0) {
    sprintf(epidemic_output_path,"%s-%s.trace",trace_output_path,stopc_description[stop_criterion]);
    epidemic_output = fopen(epidemic_output_path, "w");
    assert(epidemic_output != NULL);
  } else
    epidemic_output = NULL;

  #if PARALLEL
  #pragma omp parallel default(none)					\
  private(tid,epidemic,i,j)						\
  shared(stderr,stopc_description,p,g,ic,epidemics,sample_epidemics,data_output,\
	 stop_criterion,trace_output_path,  epidemic_output,epidemic_output_path)
  #endif
  {
  #if PARALLEL
    tid = omp_get_thread_num();
    #pragma omp for schedule(guided)
  #endif
    for (j = 0; j < epidemics; j++) {
      fprintf(stderr,"%s- thread %d: running epidemic %d with p = %f upto %s = %d %s%s ...\n",
	      tstamp(), tid, ic[j].id, p, stopc_description[stop_criterion], ic[j].bound,
	      !trace_output_path? "" : ", output: ", !trace_output_path? "" : trace_output_path);
      fflush(stderr);
      
      for (i = 1; i <= sample_epidemics; i++) {
	epidemic = epidemic_new(p, g, ic+j, epidemic_output);
	
	if (data_output) {
	  fprintf(data_output,
		  "Epidemic %d #%d: started at t = %d with %d / %d ( %.2f%% ) infected nodes\n",
		  epidemic->id,i, epidemic->t, epidemic->num_infected,
		  epidemic->g->n, 100.0*(float)epidemic->num_infected/(float)epidemic->g->n);
	  fflush(data_output);
	}

	epidemic_run(epidemic);

	if (epidemic_output)
	  fflush(epidemic_output);

	if (data_output) {
	  fprintf(data_output, 
		  "Epidemic %d #%d: stopped at t = %d with %d / %d ( %.2f%% ) infected nodes and %d links\n",
		  epidemic->id,i, epidemic->t, epidemic->num_infected,
		  epidemic->g->n, 100.0*(float)epidemic->num_infected/(float)epidemic->g->n,
		  epidemic->cascade_links);
	  fflush(data_output);
	}

	epidemic_destroy(epidemic);
      }
      ic_clean(ic+j);
    }
  }
  // close global epidemic_output
  if (epidemic_output)
    fclose(epidemic_output);

  // clean up and exit
  if (data_output && data_output != stdout)
    fclose(data_output);

  fputc('\n', stderr);
  fprintf(stderr,"%s\nDone.\n", tstamp());
  fflush(stderr);
  free_graph(g);
  free(ic);
  return 0;
}
