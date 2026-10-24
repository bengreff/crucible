#pragma once
// The drivers' restart files (docs/evidence/RESTART.md), around core/checkpoint.hpp's stream.
//   CRUCIBLE_CHECKPOINT=<stem>       write <stem>.ckpt at a step boundary every CRUCIBLE_CHECKPOINT_EVERY wall seconds
//                                    (default 1800) and at the run's end: to <stem>.ckpt.tmp, fsync, the previous
//                                    file kept as <stem>.ckpt.prev (a hard link), renamed over <stem>.ckpt, the
//                                    directory fsynced, so a kill at any point leaves a readable file.
//   CRUCIBLE_RESTART=auto            read <stem>.ckpt, or <stem>.ckpt.prev if the latest fails its check; neither
//                                    there: a fresh start from t = 0. CRUCIBLE_RESTART=<file>: that file (missing: an error).
//   CRUCIBLE_CHECKPOINT_STEPS=<n>    (tests) also every n driver steps.
//   CRUCIBLE_KILL_AT_STEP=<n>        (tests) leave at driver step n as a kill would: _Exit(75), no flush, no done-file.
//   CRUCIBLE_STATE_OUT=<file>        (tests) at the run's end, the engine objects' checkpoint fields alone (no
//                                    header, no driver fields), to compare two runs' complete states byte for byte.
// A file holds a header (magic, version, body length, FNV-1a hash of the body) and the body, which opens with the
// run's signature: its arguments and CRUCIBLE_* environment less these variables and those a restart may change
// (CRUCIBLE_T1_THREADS and a driver's thread argument: the result does not depend on the thread count;
// CRUCIBLE_T1_END: a run may be extended).
#include "core/checkpoint.hpp"

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

extern char** environ;

namespace crucible {

class CheckpointFile {
public:
  // free: the argument positions a restart may change (a thread count).
  CheckpointFile(int argc, char** argv, std::vector<int> free = {}) {
    if (const char* v = std::getenv("CRUCIBLE_CHECKPOINT")) stem_ = v;
    if (const char* v = std::getenv("CRUCIBLE_CHECKPOINT_EVERY")) every_ = std::atof(v);
    if (const char* v = std::getenv("CRUCIBLE_CHECKPOINT_STEPS")) everySteps_ = std::atol(v);
    if (const char* v = std::getenv("CRUCIBLE_KILL_AT_STEP")) killAt_ = std::atol(v);
    if (const char* v = std::getenv("CRUCIBLE_RESTART")) restart_ = v;
    if (!restart_.empty() && restart_ == "auto" && stem_.empty())
      throw std::invalid_argument("CRUCIBLE_RESTART=auto needs CRUCIBLE_CHECKPOINT=<stem>");
    std::vector<std::string> env;
    for (char** e = environ; *e; ++e) {
      const std::string s(*e);
      const std::string name = s.substr(0, s.find('='));
      if (name.rfind("CRUCIBLE_", 0) != 0) continue;
      if (name.rfind("CRUCIBLE_CHECKPOINT", 0) == 0 || name == "CRUCIBLE_RESTART" || name == "CRUCIBLE_KILL_AT_STEP" || name == "CRUCIBLE_STATE_OUT" ||
          name == "CRUCIBLE_T1_THREADS" || name == "CRUCIBLE_T1_END")
        continue;
      env.push_back(s);
    }
    std::sort(env.begin(), env.end());
    for (int a = 0; a < argc; ++a)
      signature_ += std::string(a ? " " : "") + (std::find(free.begin(), free.end(), a) != free.end() ? "*" : argv[a]);
    for (const auto& s : env) signature_ += "\n" + s;
    last_ = std::chrono::steady_clock::now();
  }
  [[nodiscard]] bool enabled() const { return !stem_.empty(); }
  [[nodiscard]] const std::string& signature() const { return signature_; }
  // Whether a checkpoint is due after driver step `step`.
  [[nodiscard]] bool due(long step) const {
    if (!enabled()) return false;
    if (everySteps_ > 0 && step % everySteps_ == 0) return true;
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - last_).count() >= every_;
  }
  void killPoint(long step) const {
    if (killAt_ > 0 && step == killAt_) std::_Exit(75);
  }
  // Writes the checkpoint: the signature, then body's fields.
  void write(const std::function<void(CheckpointWriter&)>& body) {
    std::ostringstream os(std::ios::binary);
    CheckpointWriter w(os);
    w.field(signature_);
    body(w);
    const std::string bytes = os.str();
    const std::uint64_t header[4] = {kMagic, kVersion, bytes.size(), hash(bytes)};
    const std::string latest = stem_ + ".ckpt", previous = latest + ".prev", tmp = latest + ".tmp";
    const int fd = ::open(tmp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) throw std::runtime_error("checkpoint: cannot create " + tmp);
    auto put = [&](const char* p, std::size_t n) {
      while (n) {
        const auto k = ::write(fd, p, n);
        if (k <= 0) { ::close(fd); throw std::runtime_error("checkpoint: write to " + tmp + " failed"); }
        p += k, n -= static_cast<std::size_t>(k);
      }
    };
    put(reinterpret_cast<const char*>(header), sizeof header);
    put(bytes.data(), bytes.size());
    if (::fsync(fd) != 0 || ::close(fd) != 0) throw std::runtime_error("checkpoint: fsync of " + tmp + " failed");
    std::error_code ec;
    if (std::filesystem::exists(latest)) {
      std::filesystem::remove(previous, ec);
      std::filesystem::create_hard_link(latest, previous, ec);
      if (ec) throw std::runtime_error("checkpoint: cannot keep " + previous + ": " + ec.message());
    }
    std::filesystem::rename(tmp, latest);
    std::string dir = std::filesystem::path(latest).parent_path().string();
    if (dir.empty()) dir = ".";
    if (const int d = ::open(dir.c_str(), O_RDONLY); d >= 0) ::fsync(d), ::close(d);
    last_ = std::chrono::steady_clock::now();
    ++written_;
  }
  // Reads the restart file into body; false for a fresh start (no restart asked, or auto with no file).
  bool read(const std::function<void(CheckpointReader&)>& body) {
    if (restart_.empty()) return false;
    std::vector<std::string> files;
    if (restart_ == "auto") files = {stem_ + ".ckpt", stem_ + ".ckpt.prev"};
    else files = {restart_};
    bool any = false;
    std::string why;
    for (const auto& f : files) {
      if (!std::filesystem::exists(f)) continue;
      any = true;
      try {
        std::ifstream in(f, std::ios::binary);
        const std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::uint64_t header[4];
        if (all.size() < sizeof header) throw std::runtime_error("shorter than its header");
        std::memcpy(header, all.data(), sizeof header);
        if (header[0] != kMagic || header[1] != kVersion) throw std::runtime_error("not a checkpoint of this version");
        if (header[2] != all.size() - sizeof header) throw std::runtime_error("truncated");
        const std::string bytes = all.substr(sizeof header);
        if (hash(bytes) != header[3]) throw std::runtime_error("its hash does not match");
        std::istringstream is(bytes, std::ios::binary);
        CheckpointReader r(is);
        r.same(signature_, "the run's arguments and CRUCIBLE_* environment");
        body(r);
        r.end();
        std::printf("restart: from %s\n", f.c_str());
        std::fflush(stdout);
        return true;
      } catch (const std::exception& e) {
        why += "\n  " + f + ": " + e.what();
        std::printf("restart: %s not usable: %s\n", f.c_str(), e.what());
        std::fflush(stdout);
      }
    }
    if (!any && restart_ == "auto") {
      std::printf("restart: no checkpoint, fresh start\n");
      std::fflush(stdout);
      return false;
    }
    if (!any) throw std::runtime_error("restart: no file " + restart_);
    throw std::runtime_error("restart: no usable checkpoint:" + why);
  }
  [[nodiscard]] long written() const { return written_; }

private:
  static constexpr std::uint64_t kMagic = 0x54504b4349435243ull;  // "CRCICKPT"
  static constexpr std::uint64_t kVersion = 1;
  std::string stem_, restart_, signature_;
  double every_{1800};
  long everySteps_{0}, killAt_{0}, written_{0};
  std::chrono::steady_clock::time_point last_;
  static std::uint64_t hash(const std::string& s) {
    std::uint64_t h = 1469598103934665603ull;
    for (const unsigned char c : s) h = (h ^ c) * 1099511628211ull;
    return h;
  }
};

// CRUCIBLE_STATE_OUT's file, if asked for.
inline void writeStateOut(const std::function<void(CheckpointWriter&)>& body) {
  const char* path = std::getenv("CRUCIBLE_STATE_OUT");
  if (!path) return;
  std::ofstream out(path, std::ios::binary);
  CheckpointWriter w(out);
  body(w);
  out.close();
  if (!out) throw std::runtime_error(std::string("cannot write ") + path);
}
// Cuts an output file back to `size` bytes (its length at the checkpoint) and opens it for appending.
inline std::FILE* reopenAt(const std::string& path, std::uint64_t size) {
  std::filesystem::resize_file(path, size);
  std::FILE* f = std::fopen(path.c_str(), "a");
  if (!f) throw std::runtime_error("cannot reopen " + path);
  return f;
}
// An open output file's length after flushing it.
inline std::uint64_t flushedSize(std::FILE* f) {
  std::fflush(f);
  return static_cast<std::uint64_t>(std::ftell(f));
}

}  // namespace crucible
