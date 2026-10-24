#pragma once
#include <cstdint>
#include <cstring>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace crucible {

// A restart checkpoint's stream (docs/evidence/RESTART.md). Each class lists its persistent state once, in a
// template over the stream (Flow::checkpoint, Amr::checkpoint, ...), so the writer and the reader walk the same
// fields in the same order: trivially copyable values and vectors of them as raw bytes, a vector after its
// length, and a named mark before each object's part so a reader out of step stops at once. The bytes are this
// build's: a checkpoint is read back by the same binary, not a portable format.
class CheckpointWriter {
public:
    explicit CheckpointWriter(std::ostream& out):out_(out) {}
    static constexpr bool reading=false;
    template<class T> void field(const T& v) {
        static_assert(std::is_trivially_copyable_v<T>);
        out_.write(reinterpret_cast<const char*>(&v),sizeof(T));
        check();
    }
    template<class T> void field(const std::vector<T>& v) {
        static_assert(std::is_trivially_copyable_v<T>);
        field(static_cast<std::uint64_t>(v.size()));
        if(!v.empty()) out_.write(reinterpret_cast<const char*>(v.data()),static_cast<std::streamsize>(v.size()*sizeof(T)));
        check();
    }
    void field(const std::string& s) {
        field(static_cast<std::uint64_t>(s.size()));
        out_.write(s.data(),static_cast<std::streamsize>(s.size()));
        check();
    }
    // A value the reader must find equal to its own (a mesh size, a species count, a mark's name).
    template<class T> void same(const T& v,const char*) { field(v); }
    void mark(const char* name) { field(std::string(name)); }
    template<class T> void object(const T& o) { o.writeCheckpoint(*this); }
private:
    std::ostream& out_;
    void check() const { if(!out_) throw std::runtime_error("checkpoint: write failed"); }
};

class CheckpointReader {
public:
    explicit CheckpointReader(std::istream& in):in_(in) {}
    static constexpr bool reading=true;
    template<class T> void field(T& v) {
        static_assert(std::is_trivially_copyable_v<T>);
        in_.read(reinterpret_cast<char*>(&v),sizeof(T));
        check();
    }
    template<class T> void field(std::vector<T>& v) {
        static_assert(std::is_trivially_copyable_v<T>);
        std::uint64_t n{};
        field(n);
        if(n>(std::uint64_t{1}<<40)/sizeof(T)) throw std::runtime_error("checkpoint: implausible vector length");
        v.resize(static_cast<std::size_t>(n));
        if(n) in_.read(reinterpret_cast<char*>(v.data()),static_cast<std::streamsize>(n*sizeof(T)));
        check();
    }
    void field(std::string& s) {
        std::uint64_t n{};
        field(n);
        if(n>(std::uint64_t{1}<<30)) throw std::runtime_error("checkpoint: implausible string length");
        s.resize(static_cast<std::size_t>(n));
        if(n) in_.read(s.data(),static_cast<std::streamsize>(n));
        check();
    }
    template<class T> void same(const T& v,const char* what) {
        T read{};
        field(read);
        if constexpr(std::is_trivially_copyable_v<T>) {
            if(std::memcmp(&read,&v,sizeof(T))!=0) throw std::runtime_error(std::string("checkpoint: ")+what+" differs from this run's");
        } else if(read!=v) throw std::runtime_error(std::string("checkpoint: ")+what+" differs from this run's");
    }
    void mark(const char* name) { same(std::string(name),"a section mark"); }
    template<class T> void object(T& o) { o.readCheckpoint(*this); }
    // The whole stream consumed (a reader that stops short has read a different layout).
    void end() {
        in_.peek();
        if(!in_.eof()) throw std::runtime_error("checkpoint: bytes left after the last field");
    }
private:
    std::istream& in_;
    void check() const { if(!in_) throw std::runtime_error("checkpoint: truncated or unreadable"); }
};

}  // namespace crucible
