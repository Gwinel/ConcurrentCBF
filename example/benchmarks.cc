#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <pthread.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <sched.h>
#include <inttypes.h>
#include <sys/time.h>
#include <unistd.h>
#include <malloc.h>

#include "cuckoofilter.h"

using cuckoofilter::CuckooFilter;

int init_seq = 0;
size_t num_threads = 1;
size_t num_elements = 2048;
int duration = 1000;
double filling_rate = 0.5;
double update_rate = 0.1;
double put_rate = 0.1;
double get_rate = 0.9;
int print_vals_num = 0;
size_t pf_vals_num = 8191;
size_t density = 50;

__thread unsigned long *seeds;
size_t rand_max;
#define rand_min 1
#define LOAD_FACTOR 2
#define RANGE_EXACT 0

static volatile int stop;

typedef uint64_t ticks;

volatile ticks *putting_succ;
volatile ticks *putting_fail;
volatile ticks *getting_succ;
volatile ticks *getting_fail;
volatile ticks *removing_succ;
volatile ticks *removing_fail;
volatile ticks *putting_count;
volatile ticks *putting_count_succ;
volatile ticks *getting_count;
volatile ticks *getting_count_succ;
volatile ticks *removing_count;
volatile ticks *removing_count_succ;
volatile ticks *total;

#define LFENCE asm volatile ("lfence")

#ifdef COMPUTE_LATENCY
#  define START_TS(s)       \
    asm volatile ("");        \
    start_acq = getticks();     \
    LFENCE;
#  define END_TS(s, i)        \
    asm volatile ("");        \
    end_acq = getticks();     \
    asm volatile ("");
#  define END_TS_ELSE(s, i, inc)    \
  else            \
    {           \
      END_TS(s, i);       \
      ADD_DUR(inc);       \
    }
#  define ADD_DUR(tar) tar += (end_acq - start_acq - correction)
#else
#define START_TS(s)
#define END_TS(s, i)
#define END_TS_ELSE(s, i, inc)
#define ADD_DUR(tar)
#endif

#define MEM_BARRIER __sync_synchronize()

typedef struct barrier
{
  pthread_cond_t complete;
  pthread_mutex_t mutex;
  int count;
  int crossing;
} barrier_t;

void barrier_init(barrier_t *b, int n)
{
  pthread_cond_init(&b->complete, NULL);
  pthread_mutex_init(&b->mutex, NULL);
  b->count = n;
  b->crossing = 0;
}

void barrier_cross(barrier_t *b)
{
  pthread_mutex_lock(&b->mutex);
  /* One more thread through */
  b->crossing++;
  /* If not all here, wait */
  if (b->crossing < b->count) {
    pthread_cond_wait(&b->complete, &b->mutex);
  } else {
    pthread_cond_broadcast(&b->complete);
    /* Reset for next time */
    b->crossing = 0;
  }
  pthread_mutex_unlock(&b->mutex);
}

barrier_t barrier, barrier_global;

// CuckooFilter template parameters
#define TYPE uint64_t
#define BITS 12

typedef struct thread_data
{
  uint32_t id;
  CuckooFilter<TYPE, BITS> *filter;
} thread_data_t;

static inline int is_power_of_two (unsigned int x) 
{
  return ((x != 0) && !(x & (x - 1)));
}

/// Round up to next higher power of 2 (return x if it's already a power
/// of 2) for 32-bit numbers
static inline uint64_t pow2roundup(uint64_t x)
{
  if (x == 0) return 1;
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return x + 1;
}

static inline ticks getticks(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ((unsigned long long)lo) | (((unsigned long long)hi) << 32);
}

// getticks needs to have a correction because the call itself takes a
// significant number of cycles and skewes the measurement
static inline ticks getticks_correction_calc() 
{
#define GETTICKS_CALC_REPS 5000000
  ticks t_dur = 0;
  uint32_t i;
  for (i = 0; i < GETTICKS_CALC_REPS; i++) {
    ticks t_start = getticks();
    ticks t_end = getticks();
    t_dur += t_end - t_start;
  }
  //    printf("corr in float %f\n", (t_dur / (double) GETTICKS_CALC_REPS));
  ticks getticks_correction = (ticks)(t_dur / (double) GETTICKS_CALC_REPS);
  return getticks_correction;
}

static inline unsigned long* seed_rand() 
{
  unsigned long* seeds;
  seeds = (unsigned long*) memalign(64, 64);
  seeds[0] = getticks() % 0x992123E456789LL;
  seeds[1] = getticks() % 0x22136D2436069LL;
  seeds[2] = getticks() % 0x2119F521288629LL;
  return seeds;
}

// Marsaglia's xorshf generator
static inline unsigned long xorshf96(unsigned long* x, unsigned long* y, unsigned long* z) 
{          // period 2^96-1
  unsigned long t;
  (*x) ^= (*x) << 16;
  (*x) ^= (*x) >> 5;
  (*x) ^= (*x) << 1;

  t = *x;
  (*x) = *y;
  (*y) = *z;
  (*z) = t ^ (*x) ^ (*y);

  return *z;
}

#define my_random xorshf96

void*
test(void* thread)
{
  thread_data_t* td = (thread_data_t*) thread;
  uint32_t ID = td->id;
  //phys_id = the_cores[ID];
  //set_cpu(phys_id);

  CuckooFilter<TYPE, BITS> *filter = td->filter;

#if defined(COMPUTE_LATENCY)
  volatile ticks my_putting_succ = 0;
  volatile ticks my_putting_fail = 0;
  volatile ticks my_getting_succ = 0;
  volatile ticks my_getting_fail = 0;
  volatile ticks my_removing_succ = 0;
  volatile ticks my_removing_fail = 0;
#endif

  uint64_t my_putting_count = 0;
  uint64_t my_getting_count = 0;
  uint64_t my_removing_count = 0;

  uint64_t my_putting_count_succ = 0;
  uint64_t my_getting_count_succ = 0;
  uint64_t my_removing_count_succ = 0;

  seeds = seed_rand();

#ifdef COMPUTE_LATENCY
  volatile ticks start_acq, end_acq;
  volatile ticks correction = getticks_correction_calc();
#endif

  barrier_cross(&barrier);

  uint64_t key;
  uint64_t scale_rem = (uint64_t) (update_rate * ((uint64_t) -1));
  uint64_t scale_put = (uint64_t) (put_rate * ((uint64_t) -1));

  size_t i;
  size_t num_slots = filter->TotalSlots();
  size_t num_elements_init = num_slots * filling_rate;
  printf("num_elements_init: %lu\n", num_elements_init);
  if (init_seq) {
    //size_t num_slots = filter->TotalSlots();
    //size_t num_elements_init = num_slots * filling_rate;
    if (ID == 0) {
      printf("** will initialize sequentially, keys 1 .. %zu\n", num_elements_init);
      size_t progress_step = (num_elements_init / 10) - 1;
      size_t progress = progress_step;
      size_t progress_100 = 10;
      printf("Progress: 0%%  "); fflush(stdout);
      for (i = 0; i < num_elements_init; i++) {
        if (unlikely(i == progress)) {
          progress += progress_step;
          printf("%zu%%  ", progress_100); fflush(stdout);
          progress_100 += 10;
        }
        key = i;//rand_max - i - 1;
        filter->Add(key);
      }
      std::cout << filter->Info() << std::endl;
      printf("\n");
    }
  } else {    
    //uint64_t num_elems_thread = (uint64_t) (num_elements * filling_rate / num_threads);
    size_t num_slots = filter->TotalSlots();
    uint64_t num_elems_thread = (uint64_t) (num_slots * filling_rate / num_threads);
    int64_t missing = (uint64_t) (num_elements * filling_rate) - (num_elems_thread * num_threads);
    if (ID < missing) {
      num_elems_thread++;
    }

    for (i = 0; i < num_elems_thread; i++) {
      key = (my_random(&(seeds[0]), &(seeds[1]), &(seeds[2])) & (rand_max)) + rand_min;
      //uint64_t num = rand();
      //key = rand();//(num & (rand_max)) + rand_min;

      //printf("Adding key (%lu) to filter, index: %lu\n", key, i);
      if (filter->Add(key) == cuckoofilter::Retry) {
        //printf("duplicate key.\n");
        i--;
      }
    }
  }
  MEM_BARRIER;
  //std::cout << "filling_rate: " << filling_rate << std::endl;
  //std::cout << "num_elements: " << num_elements << std::endl;
  //std::cout << "rand_max: " << rand_max << std::endl;
  //std::cout << filter->Info() << std::endl;
  //exit(0);

  barrier_cross(&barrier);
  barrier_cross(&barrier_global);

  printf("before while loop in thread %u...\n", ID);
  i = 0;
  while (stop == 0) {
    uint64_t c = (uint64_t)(my_random(&(seeds[0]), &(seeds[1]), &(seeds[2])));
    key = (c & (num_elements_init));//(c & rand_max) + rand_min;

    if (unlikely(c <= scale_put)) {
      cuckoofilter::Status res;
      START_TS(1);
      res = filter->Add(key);
      if (res == cuckoofilter::Ok) {
        END_TS(1, my_putting_count_succ);
        ADD_DUR(my_putting_succ);
        my_putting_count_succ++;
      }
      END_TS_ELSE(4, my_putting_count - my_putting_count_succ,
          my_putting_fail);
      my_putting_count++;
    } else if (unlikely(c <= scale_rem)) {
      cuckoofilter::Status res;
      START_TS(2);
      res = filter->Delete(key);
      if (res == cuckoofilter::Ok) {
        END_TS(2, my_removing_count_succ);
        ADD_DUR(my_removing_succ);
        my_removing_count_succ++;
      }
      END_TS_ELSE(5, my_removing_count - my_removing_count_succ,
          my_removing_fail);
      my_removing_count++;
    } else {
      cuckoofilter::Status res;
      START_TS(0);
      res = filter->Contain(c % (num_elements_init - 1));
      if (res == cuckoofilter::Ok) {
        END_TS(0, my_getting_count_succ);
        ADD_DUR(my_getting_succ);
        my_getting_count_succ++;
      }
      END_TS_ELSE(3, my_getting_count - my_getting_count_succ,
          my_getting_fail);
      my_getting_count++;
    }
  }
  printf("after while loop in thread %u...\n", ID);

  /* printf("gets: %-10llu / succ: %llu\n", num_get, num_get_succ); */
  /* printf("rems: %-10llu / succ: %llu\n", num_rem, num_rem_succ); */
  //barrier_cross(&barrier);

#ifdef COMPUTE_LATENCY
  putting_succ[ID] += my_putting_succ;
  putting_fail[ID] += my_putting_fail;
  getting_succ[ID] += my_getting_succ;
  getting_fail[ID] += my_getting_fail;
  removing_succ[ID] += my_removing_succ;
  removing_fail[ID] += my_removing_fail;
#endif
  putting_count[ID] += my_putting_count;
  getting_count[ID] += my_getting_count;
  removing_count[ID]+= my_removing_count;

  putting_count_succ[ID] += my_putting_count_succ;
  getting_count_succ[ID] += my_getting_count_succ;
  removing_count_succ[ID]+= my_removing_count_succ;

  pthread_exit(NULL);
}

int
main(int argc, char **argv)
{
  struct option long_options[] = {
    // These options don't set a flag
    {"help",                no_argument,       NULL, 'h'},
    {"duration",            required_argument, NULL, 'd'},
    {"initial-size",        required_argument, NULL, 'i'},
    {"init-seq",            required_argument, NULL, 's'},
    {"num-threads",         required_argument, NULL, 'n'},
    {"range",               required_argument, NULL, 'r'},
    {"update-rate",         required_argument, NULL, 'u'},
    {"print-vals",          required_argument, NULL, 'v'},
    {"table-density",       required_argument, NULL, 'f'},
    {NULL, 0, NULL, 0}
  };

  size_t initial = 1024, range = 2048, update = 20, put = 10;
  double load_factor = LOAD_FACTOR;
  int put_explicit = 0;

  int i, c;
  while (1) {
    i = 0;
    c = getopt_long(argc, argv, "hAf:d:i:n:r:su:m:a:l:p:v:f:", long_options, &i);

    if (c == -1)
      break;

    if (c == 0 && long_options[i].flag == 0)
      c = long_options[i].val;

    switch(c) {
    case 0:
      /* Flag is automatically set */
      break;
    case 'h':
      printf("Cuckoo Filter test "
      "\n"
      "Usage:\n"
      "  executable [options...]\n"
      "\n"
      "Options:\n"
      "  -h, --help\n"
      "        Print this message\n"
      "  -d, --duration <int>\n"
      "        Test duration in milliseconds\n"
      "  -i, --initial-size <int>\n"
      "        Number of elements to insert before test\n"
      "  -s, --init-seq\n"
      "        Initialize the hash table sequentially, with keys from 1 .. initial\n"
      "  -n, --num-threads <int>\n"
      "        Number of threads\n"
      "  -r, --range <int>\n"
      "        Range of integer values inserted in set\n"
      "  -u, --update-rate <int>\n"
      "        Percentage of update transactions\n"
      "  -p, --put-rate <int>\n"
      "        Percentage of put update transactions (should be less than percentage of updates)\n"
      "  -l, --load-factor <int>\n"
      "        Number of elements per bucket (initially)\n"
      "  -v, --print-vals <int>\n"
      "        When using detailed profiling, how many values to print.\n"
      "  -f, --table-density <int>\n"
      "        Table density.\n"
      );
      exit(0);
    case 'd':
      duration = atoi(optarg);
      break;
    case 'i':
      initial = atol(optarg);
      break;
    case 's':
      init_seq = 1;
      break;
    case 'n':
      num_threads = atoi(optarg);
      break;
    case 'r':
      range = atol(optarg);
      break;
    case 'u':
      update = atoi(optarg);
      break;
    case 'p':
      put_explicit = 1;
      put = atoi(optarg);
      break;
    case 'l':
      load_factor = atof(optarg);
      break;
    case 'v':
      print_vals_num = atoi(optarg);
      break;
    case 'f':
      density = atoi(optarg);
      break;
    case '?':
      default:
      printf("Use -h or --help for help\n");
      exit(1);
    }
  }

#if RANGE_EXACT != 1
  if (!is_power_of_two(initial)) {
    size_t initial_pow2 = pow2roundup(initial);
    printf("** rounding up initial (to make it power of 2): old: %zu / new: %zu\n", initial, initial_pow2);
    initial = initial_pow2;
  }
#endif

  if (range < initial)
    range = 2 * initial;

#if RANGE_EXACT != 1
  if (!is_power_of_two(range)) {
    size_t range_pow2 = pow2roundup(range);
    printf("** rounding up range (to make it power of 2): old: %zu / new: %zu\n", range, range_pow2);
    range = range_pow2;
  }
#endif

  printf("## Initial: %zu / Range: %zu\n", initial, range);

  if (put > update)
    put = update;

  num_elements = range;
  filling_rate = density / 100.0; //(double) initial / range;
  update_rate = update / 100.0;
  if (put_explicit)
    put_rate = put / 100.0;
  else
    put_rate = update_rate / 2;

  get_rate = 1 - update_rate;

#if RANGE_EXACT == 1
  rand_max = num_elements;
#else
  rand_max = num_elements - 1;
#endif

  struct timeval start, end;
  struct timespec timeout;
  timeout.tv_sec = duration / 1000;
  timeout.tv_nsec = (duration % 1000) * 1000000;

  stop = 0;

  /* Initialize the filter */
  CuckooFilter<TYPE, BITS> *filter = new CuckooFilter<TYPE, BITS>(num_elements);

  /* Initializes the local data */
  putting_succ = (ticks *) calloc(num_threads, sizeof(ticks));
  putting_fail = (ticks *) calloc(num_threads, sizeof(ticks));
  getting_succ = (ticks *) calloc(num_threads, sizeof(ticks));
  getting_fail = (ticks *) calloc(num_threads, sizeof(ticks));
  removing_succ = (ticks *) calloc(num_threads, sizeof(ticks));
  removing_fail = (ticks *) calloc(num_threads, sizeof(ticks));
  putting_count = (ticks *) calloc(num_threads, sizeof(ticks));
  putting_count_succ = (ticks *) calloc(num_threads, sizeof(ticks));
  getting_count = (ticks *) calloc(num_threads, sizeof(ticks));
  getting_count_succ = (ticks *) calloc(num_threads, sizeof(ticks));
  removing_count = (ticks *) calloc(num_threads, sizeof(ticks));
  removing_count_succ = (ticks *) calloc(num_threads, sizeof(ticks));

  pthread_t threads[num_threads];
  pthread_attr_t attr;
  int rc;
  void *status;

  barrier_init(&barrier_global, num_threads + 1);
  barrier_init(&barrier, num_threads);

  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  thread_data_t* tds = (thread_data_t*) malloc(num_threads * sizeof(thread_data_t));

  printf("before creating threads...\n");
  size_t t;
  for (t = 0; t < num_threads; t++) {
    tds[t].id = t;
    tds[t].filter = filter;
    rc = pthread_create(&threads[t], &attr, test, tds + t);
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  /* Free attribute and wait for the other threads */
  pthread_attr_destroy(&attr);

  barrier_cross(&barrier_global);
  gettimeofday(&start, NULL);
  nanosleep(&timeout, NULL);

  stop = 1;
  gettimeofday(&end, NULL);
  duration = (end.tv_sec * 1000 + end.tv_usec / 1000) - (start.tv_sec * 1000 + start.tv_usec / 1000);
  printf("%d\n", duration);

  for (t = 0; t < num_threads; t++) {
    rc = pthread_join(threads[t], &status);
    if (rc) {
      printf("ERROR; return code from pthread_join() is %d\n", rc);
      exit(-1);
    }
  }

  free(tds);

  volatile ticks putting_suc_total = 0;
  volatile ticks putting_fal_total = 0;
  volatile ticks getting_suc_total = 0;
  volatile ticks getting_fal_total = 0;
  volatile ticks removing_suc_total = 0;
  volatile ticks removing_fal_total = 0;
  volatile uint64_t putting_count_total = 0;
  volatile uint64_t putting_count_total_succ = 0;
  volatile uint64_t getting_count_total = 0;
  volatile uint64_t getting_count_total_succ = 0;
  volatile uint64_t removing_count_total = 0;
  volatile uint64_t removing_count_total_succ = 0;

  for (t = 0; t < num_threads; t++) {
    putting_suc_total += putting_succ[t];
    putting_fal_total += putting_fail[t];
    getting_suc_total += getting_succ[t];
    getting_fal_total += getting_fail[t];
    removing_suc_total += removing_succ[t];
    removing_fal_total += removing_fail[t];
    putting_count_total += putting_count[t];
    putting_count_total_succ += putting_count_succ[t];
    getting_count_total += getting_count[t];
    getting_count_total_succ += getting_count_succ[t];
    removing_count_total += removing_count[t];
    removing_count_total_succ += removing_count_succ[t];
  }

#define LLU long long unsigned int

//   int pr = (int) (putting_count_total_succ - removing_count_total_succ);
//   int size_after = clht_size(hashtable->ht);
// #if defined(DEBUG)
//   printf("puts - rems  : %d\n", pr);
// #endif
//   assert(size_after == (initial + pr));

  printf("    : %-10s | %-10s | %-11s | %s\n", "total", "success", "succ %", "total %");
  uint64_t total = putting_count_total + getting_count_total + removing_count_total;
  double putting_perc = 100.0 * (1 - ((double)(total - putting_count_total) / total));
  double getting_perc = 100.0 * (1 - ((double)(total - getting_count_total) / total));
  double removing_perc = 100.0 * (1 - ((double)(total - removing_count_total) / total));
  printf("puts: %-10llu | %-10llu | %10.1f%% | %.1f%%\n", (LLU) putting_count_total, 
   (LLU) putting_count_total_succ,
   (1 - (double) (putting_count_total - putting_count_total_succ) / putting_count_total) * 100,
   putting_perc);
  printf("gets: %-10llu | %-10llu | %10.1f%% | %.1f%%\n", (LLU) getting_count_total, 
   (LLU) getting_count_total_succ,
   (1 - (double) (getting_count_total - getting_count_total_succ) / getting_count_total) * 100,
   getting_perc);
  printf("rems: %-10llu | %-10llu | %10.1f%% | %.1f%%\n", (LLU) removing_count_total, 
   (LLU) removing_count_total_succ,
   (1 - (double) (removing_count_total - removing_count_total_succ) / removing_count_total) * 100,
   removing_perc);

  std::cout << filter->Info() << "\n";

  double throughput = (putting_count_total + getting_count_total + removing_count_total) * 1000.0 / duration;
  printf("#txs %zu\t(%-10.0f)\n", num_threads, throughput);
  printf("#Mops %.3f\n", throughput / 1e6);

  pthread_exit(NULL);
  return 0;
}