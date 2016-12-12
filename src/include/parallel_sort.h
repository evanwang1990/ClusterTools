#ifndef CLUSTER_TOOLS_PARALLEL_SORT_
#define CLUSTER_TOOLS_PARALLEL_SORT_

#include <thread>
#include <vector>
#include <algorithm>

using namespace std;

template <typename T>
void Sort(T* arr, size_t len, size_t grainsize)
{
  if (len < grainsize) {
    sort(data, data+len)
  } else {
    thread thr(Sort<T>, data, len/2, grainsize);
    Sort(data+len/2, len-len/2, grainsize);
    thr.join();
    inplace_merge(data, data+len/2, data+len);
  }
}

template<T>
void ParallelSort(vector<T>* arr, int num_threads)
{
  size_t grainsize = max(arr->size()/num_threads + 5, (size_t)1024*16);
  Sort(arr->data(), arr->size(), grainsize);
}

