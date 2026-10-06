#pragma once
#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace crucible {
// A fixed set of threads kept for the life of the object, so a parallel loop costs a wake-up and
// not a thread start (starting them cost about 84 us per call, docs/evidence/TABLE_A.md). The
// calling thread is worker 0. Idle workers spin for 100 us and then sleep until the next job.
// Not reentrant: a job must not call run() on its own pool.
class Pool {
public:
    explicit Pool(int threads=1);
    ~Pool();
    Pool(const Pool&)=delete;
    Pool& operator=(const Pool&)=delete;
    [[nodiscard]] int threads() const { return threads_; }
    // Calls job(worker) once on every worker and returns when all have returned. The first
    // exception thrown is rethrown here, after every worker has finished.
    void run(const std::function<void(int)>& job);
    // Indices [0, n) in blocks of `block`, taken in order from a shared counter by whichever
    // worker is free: body(begin, end, worker). Which worker takes an index varies between calls,
    // so a body whose result for an index depends only on that index gives the same result on any
    // number of threads.
    template<class Body> void blocks(std::size_t n,std::size_t block,Body&& body) {
        if(threads_==1 || n<=block) { if(n) body(std::size_t{0},n,0); return; }
        std::atomic<std::size_t> next{0};
        run([&](int worker) {
            for(std::size_t begin;(begin=next.fetch_add(block,std::memory_order_relaxed))<n;)
                body(begin,std::min(n,begin+block),worker);
        });
    }
private:
    void loop(int worker);
    void execute(int worker);
    int threads_;
    std::vector<std::thread> workers_;
    const std::function<void(int)>* job_{nullptr};
    std::atomic<std::uint64_t> generation_{0};
    std::atomic<int> remaining_{0},sleepers_{0};
    std::atomic<bool> callerSleeps_{false};
    std::mutex mutex_,errorMutex_;
    std::condition_variable wake_,done_;
    std::exception_ptr error_;
    bool stop_{false};
};
} // namespace crucible
