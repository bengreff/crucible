#pragma once
#include "core/flow.hpp"
#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>

namespace crucible {
enum class RunState { Initializing, Paused, Running, PauseRequested, Failed, Stopped };
struct ControlRecord { std::uint64_t sequence{}; double time{}, totalPressure{}; };
class Session {
public:
    explicit Session(Definition definition);
    ~Session();
    Session(const Session&)=delete;
    Session& operator=(const Session&)=delete;
    void run();
    void pause();
    std::uint64_t setTotalPressure(double pressure);
    [[nodiscard]] std::shared_ptr<const FieldSnapshot> latest() const;
    [[nodiscard]] RunState status() const { return status_.load(); }
    [[nodiscard]] std::string error() const;
    [[nodiscard]] std::vector<ControlRecord> controls() const;
private:
    void work(Definition definition);
    mutable std::mutex mutex_;
    std::condition_variable changed_;
    bool runRequested_{false},stopRequested_{false};
    std::optional<std::pair<std::uint64_t,double>> pressure_;
    std::uint64_t nextSequence_{};
    std::vector<ControlRecord> controls_;
    std::string error_;
    std::atomic<RunState> status_{RunState::Initializing};
    mutable std::mutex snapshotMutex_;
    std::shared_ptr<const FieldSnapshot> latest_;
    std::thread worker_;
};
}
