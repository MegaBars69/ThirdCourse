#ifndef HEADER3
#define HEADER3
#include "algorithm.hpp"
#include "initialize_matrix.hpp"

void* thread_func(void *arg);
template <typename T>
void reduce_max(int p, T* a, int n) {
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static T* r = nullptr;
    int i;
    if(p <= 1) return;
    pthread_mutex_lock(&m);
    if(r == nullptr) {
        r = a;
    } else {
        for(i = 0; i < n; ++i) r[i] = std::max(r[i], a[i]);
    }
    ++t_in;
    if(t_in >= p) {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    } else {
        while(t_in < p) {
            pthread_cond_wait(&c_in, &m);
        }
    }
    if(r != a) {
        for(i = 0; i < n; ++i) {
            a[i] = r[i];
        }
    }
    ++t_out;
    if(t_out >= p) {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    } else {
        while (t_out < p) {
            pthread_cond_wait(&c_out, &m);
        }
    }
    pthread_mutex_unlock(&m);
}


inline void barrier(int p) {

    struct closure_barrier {
        pthread_barrier_t barrier;
        closure_barrier(int p) {
            pthread_barrier_init(&barrier, nullptr, p);
        }
        ~closure_barrier() {
            pthread_barrier_destroy(&barrier);
        }
    };

    static closure_barrier bar(p);

    pthread_barrier_wait(&bar.barrier);

}

template <typename T, size_t alignment> 
struct aligned_allocator {
    typedef T value_type;
    typedef T *pointer;

    pointer allocate(size_t n) {
        void *chunk = std::aligned_alloc(alignment,  n * sizeof(value_type));
        return static_cast<pointer>(chunk);
    }

    void deallocate(pointer p, size_t) {
        free(p);
    }

    template <typename U> 
    aligned_allocator(const aligned_allocator<U, alignment> &) {}

    aligned_allocator() {}
    aligned_allocator(const aligned_allocator &) {}

    template <typename U> 
    struct rebind { typedef aligned_allocator<U, alignment> other; };

};
#endif