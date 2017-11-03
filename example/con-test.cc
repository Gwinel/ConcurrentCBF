#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>


#include <iostream>
#include <vector>

#include <thread>

using cuckoofilter::CuckooFilter;

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

typedef struct thread_data
{
  size_t start;
  size_t end;
} thread_data_t;

size_t total_items = 1000000;
// Create a cuckoo filter where each item is of type size_t and
// use 12 bits for each item:
//    CuckooFilter<size_t, 12> filter(total_items);
// To enable semi-sorting, define the storage of cuckoo filter to be
// PackedTable, accepting keys of size_t type and making 13 bits
// for each key:
//   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
CuckooFilter<size_t, 12> filter(total_items);

void* InsertCon(void *thread)
{

    thread_data_t* td = (thread_data_t*) thread;
    printf("thread: %lu begin at %lu end at %lu\n", pthread_self(), td->start, td->end);
    barrier_cross(&barrier_global);
    for (size_t i = td->start; i < td->end; i++) {
        if (filter.Add(i) != cuckoofilter::Ok) {
            printf("InsertCon failed.\n");
            break;
        }
    }
    barrier_cross(&barrier_global);
    pthread_exit(NULL);
}

void* SearchCon(void *thread)
{
    thread_data_t* td = (thread_data_t*) thread;
    barrier_cross(&barrier_global);
    for (size_t i = td->start; i < td->end; i++) {
        if (filter.Contain(i) != cuckoofilter::Ok) {
            printf("SearchCon failed.\n");
            break;
        }
    }
    barrier_cross(&barrier_global);
    pthread_exit(NULL);
}

int main(int argc, char **argv)
{

  if (argc != 2) {
    printf("usage: %s num_thread\n", argv[0]);
    exit(0);
  }

  int num_threads = atoi(argv[1]);
  if (num_threads < 1) {
    printf("invalid input.\n");
  }

  size_t numPerThread = total_items / num_threads;
  size_t remaining = total_items % num_threads;

  std::thread nthreads[num_threads];
  size_t start = 0, end;
  timeval start_time, end_time;
  double inserts_per_sec;

  size_t vecInt[num_threads][2];
  for (int i = 0; i < num_threads; i++) {
    if (remaining == 0) {
      end = start + numPerThread;
    } else {
      end = start + numPerThread + 1;
      remaining--;
    }
    vecInt[i][0] = start;
    vecInt[i][1] = end;
    start = end;
  }

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

  long t;
  for(t = 0; t < num_threads; t++) {
      tds[t].start = vecInt[t][0];
      tds[t].end = vecInt[t][1];
      rc = pthread_create(&threads[t], &attr, InsertCon, tds + t);
      if (rc) {
        printf("ERROR; return code from pthread_create() is %d\n", rc);
        exit(-1);
      }
  }


  // for (int i = 0; i < threads; i++) {
  //   std::cout << "thread" << i << " begin=>" << vecInt[i][0] << " end=>" << vecInt[i][1] << std::endl;
  //   nthreads[i] = std::thread(&InsertCon, vecInt[i][0], vecInt[i][1]);
  // }

  // for (auto &t : nthreads) {
  //   t.join();
  // }
  /* Free attribute and wait for the other threads */

  barrier_cross(&barrier_global);
  gettimeofday(&start_time, NULL);

  barrier_cross(&barrier_global);

  gettimeofday(&end_time, NULL);
  //in milliseconds
  double elapsed = (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000;

  for(t = 0; t < num_threads; t++) {
      rc = pthread_join(threads[t], &status);
      if (rc) {
        printf("ERROR; return code from pthread_join() is %d\n", rc);
        exit(-1);
      }
  }

  inserts_per_sec = 1000 * total_items / elapsed;

  std::cout << "total_items is: " << total_items << std::endl;
  std::cout << "inserts per sec (M) : " << inserts_per_sec / 1000000 << std::endl;

  for(t = 0; t < num_threads; t++) {
      tds[t].start = vecInt[t][0];
      tds[t].end = vecInt[t][1];
      rc = pthread_create(&threads[t], &attr, SearchCon, tds + t);
      if (rc) {
        printf("ERROR; return code from pthread_create() is %d\n", rc);
        exit(-1);
      }
  }

  barrier_cross(&barrier_global);
  gettimeofday(&start_time, NULL);

  barrier_cross(&barrier_global);

  gettimeofday(&end_time, NULL);
  //in milliseconds
  elapsed = (end_time.tv_sec - start_time.tv_sec) * 1000 + (end_time.tv_usec - start_time.tv_usec) / 1000;

  inserts_per_sec = 1000 * total_items / elapsed;
  std::cout << "search per sec (M) : " << inserts_per_sec / 1000000 << std::endl;


  for(t = 0; t < num_threads; t++) {
      rc = pthread_join(threads[t], &status);
      if (rc) {
        printf("ERROR; return code from pthread_join() is %d\n", rc);
        exit(-1);
      }
  }

  // Check if previously inserted items are in the filter, expected
  // true for all items
  // for (size_t i = 0; i < total_items; i++) {
  //   assert(filter.Contain(i) == cuckoofilter::Ok);
  // }

  // Check non-existing items, a few false positives expected
  size_t total_queries = 0;
  size_t false_queries = 0;
  for (size_t i = total_items; i < 2 * total_items; i++) {
    if (filter.Contain(i) == cuckoofilter::Ok) {
      false_queries++;
    }
    total_queries++;
  }

  pthread_attr_destroy(&attr);
  free(tds);

  // Output the measured false positive rate
  std::cout << "false positive rate is "
            << 100.0 * false_queries / total_queries << "%\n";

  std::cout << filter.Info() << "\n";

  std::cout << "num of items: " << filter.Size() << std::endl;
  std::cout << "size of table: " << filter.SizeInBytes() / (1024 * 1024) << "MB"<< std::endl;

  return 0;
}
