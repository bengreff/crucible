#include "core/session.hpp"
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace crucible {
Session::Session(Definition d):worker_([this,d]{work(d);}) {}
Session::~Session() {
    { std::lock_guard lock(mutex_);stopRequested_=true; }
    changed_.notify_all();
    if(worker_.joinable()) worker_.join();
}
void Session::run() {
    { std::lock_guard lock(mutex_);runRequested_=true; }
    changed_.notify_all();
}
void Session::pause() {
    { std::lock_guard lock(mutex_);runRequested_=false; }
    auto expected=RunState::Running;
    status_.compare_exchange_strong(expected,RunState::PauseRequested);
    changed_.notify_all();
}
std::uint64_t Session::setTotalPressure(double p) {
    if(!(p>0) || !std::isfinite(p)) throw std::invalid_argument("Invalid reservoir pressure.");
    std::lock_guard lock(mutex_);
    auto sequence=++nextSequence_;
    pressure_={{sequence,p}}; // Coalesce pending, not already-applied, control changes.
    changed_.notify_all();return sequence;
}
std::shared_ptr<const FieldSnapshot> Session::latest() const { std::lock_guard lock(snapshotMutex_); return latest_; }
std::string Session::error() const { std::lock_guard lock(mutex_);return error_; }
std::vector<ControlRecord> Session::controls() const { std::lock_guard lock(mutex_);return controls_; }
void Session::work(Definition d) {
    try {
        Flow flow(d);
        std::uint64_t sequence=0,generation=0;
        auto publish=[&] {
            auto snapshot=std::make_shared<FieldSnapshot>(flow.snapshot());
            snapshot->appliedControlSequence=sequence;snapshot->generation=++generation;
            std::lock_guard lock(snapshotMutex_); latest_=std::move(snapshot);
        };
        publish();status_=RunState::Paused;
        auto lastPublished=std::chrono::steady_clock::now();
        for(;;) {
            bool shouldRun;
            {
                std::unique_lock lock(mutex_);
                if(!runRequested_) {
                    // Publish the final accepted state before acknowledging pause.
                    lock.unlock();publish();status_=RunState::Paused;lock.lock();
                    changed_.wait(lock,[&]{return stopRequested_||runRequested_||pressure_.has_value();});
                }
                if(stopRequested_) break;
                if(pressure_) {
                    sequence=pressure_->first;flow.setTotalPressure(pressure_->second);
                    controls_.push_back({sequence,flow.time(),pressure_->second});pressure_.reset();
                    lock.unlock();publish();lock.lock();
                }
                shouldRun=runRequested_;
            }
            if(!shouldRun) continue;
            status_=RunState::Running;flow.step();
            auto now=std::chrono::steady_clock::now();
            if(now-lastPublished>=std::chrono::milliseconds(50)) {publish();lastPublished=now;}
        }
        publish();status_=RunState::Stopped;
    } catch(const std::exception& e) {
        {std::lock_guard lock(mutex_);error_=e.what();}
        status_=RunState::Failed;
    }
}
}
