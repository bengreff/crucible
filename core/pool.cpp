#include "core/pool.hpp"
#include <chrono>
#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

namespace crucible {
namespace {
// How long an idle thread spins before it sleeps: longer than the serial gaps between the loops of
// one flow step, so the step's loops do not pay a wake-up each.
constexpr auto kSpin=std::chrono::microseconds(100);
inline void relax() {
#if defined(__x86_64__) || defined(_M_X64)
    _mm_pause();
#elif defined(__aarch64__)
    asm volatile("yield");
#endif
}
// Spins until done() or kSpin has passed; returns done().
template<class Done> bool spin(Done done) {
    const auto until=std::chrono::steady_clock::now()+kSpin;
    for(;;) {
        for(int n=0;n<256;++n) { if(done()) return true; relax(); }
        if(std::chrono::steady_clock::now()>until) return done();
    }
}
}
Pool::Pool(int threads):threads_(std::max(1,threads)) {
    for(int w=1;w<threads_;++w) workers_.emplace_back([this,w]{ loop(w); });
}
Pool::~Pool() {
    { std::lock_guard<std::mutex> lock(mutex_);stop_=true;generation_.fetch_add(1); }
    wake_.notify_all();
    for(auto& t:workers_) t.join();
}
void Pool::execute(int worker) {
    try { (*job_)(worker); }
    catch(...) { std::lock_guard<std::mutex> lock(errorMutex_);if(!error_) error_=std::current_exception(); }
}
// A sleeper counts itself in sleepers_ before it checks for a job under the mutex, and run() bumps
// the generation before it reads sleepers_ (both sequentially consistent), so either the sleeper
// sees the new job or run() takes the mutex and wakes it. The caller's wait for the last worker
// follows the same pattern with done_.
void Pool::loop(int worker) {
    std::uint64_t seen=0;
    for(;;) {
        auto fresh=[&]{ return generation_.load()!=seen; };
        if(!spin(fresh)) {
            std::unique_lock<std::mutex> lock(mutex_);
            sleepers_.fetch_add(1);
            wake_.wait(lock,fresh);
            sleepers_.fetch_sub(1);
        }
        seen=generation_.load();
        if(stop_) return;
        execute(worker);
        if(remaining_.fetch_sub(1)==1 && callerSleeps_.load()) {
            { std::lock_guard<std::mutex> lock(mutex_); }
            done_.notify_one();
        }
    }
}
void Pool::run(const std::function<void(int)>& job) {
    if(threads_==1) { job(0);return; }
    job_=&job;error_=nullptr;
    remaining_.store(threads_-1);
    generation_.fetch_add(1);
    if(sleepers_.load()>0) {
        { std::lock_guard<std::mutex> lock(mutex_); }
        wake_.notify_all();
    }
    execute(0);
    auto finished=[&]{ return remaining_.load()==0; };
    if(!spin(finished)) {
        std::unique_lock<std::mutex> lock(mutex_);
        callerSleeps_.store(true);
        done_.wait(lock,finished);
        callerSleeps_.store(false);
    }
    job_=nullptr;
    if(error_) { auto e=error_;error_=nullptr;std::rethrow_exception(e); }
}
} // namespace crucible
