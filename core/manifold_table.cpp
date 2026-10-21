#include "core/manifold_table.hpp"

#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace crucible {
namespace {
constexpr char kMagic[] = "CRUCIBLE manifold table 2\n";  // 2: per-point c, the reference pressure
template <class T>
void put(std::ofstream& o, const T& v) { o.write(reinterpret_cast<const char*>(&v), sizeof v); }
template <class T>
void putVector(std::ofstream& o, const std::vector<T>& v) {
  put(o, static_cast<std::uint64_t>(v.size()));
  o.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)));
}
void putString(std::ofstream& o, const std::string& s) {
  put(o, static_cast<std::uint64_t>(s.size()));
  o.write(s.data(), static_cast<std::streamsize>(s.size()));
}
template <class T>
T get(std::ifstream& in) {
  T v{};
  in.read(reinterpret_cast<char*>(&v), sizeof v);
  return v;
}
template <class T>
std::vector<T> getVector(std::ifstream& in) {
  std::vector<T> v(get<std::uint64_t>(in));
  in.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)));
  return v;
}
std::string getString(std::ifstream& in) {
  std::string s(get<std::uint64_t>(in), '\0');
  in.read(s.data(), static_cast<std::streamsize>(s.size()));
  return s;
}
}  // namespace

void ManifoldTable::save(const std::string& path) const {
  std::ofstream o(path, std::ios::binary);
  if (!o) throw std::runtime_error("Cannot write " + path);
  o.write(kMagic, sizeof kMagic - 1);
  putString(o, provenance);
  putVector(o, fH);
  putVector(o, zN);
  put(o, nEta);
  put(o, nLnRho);
  put(o, lnRhoMin);
  put(o, lnRhoMax);
  put(o, referencePressure);
  put(o, static_cast<std::uint64_t>(reactingSpecies.size()));
  for (const auto& s : reactingSpecies) putString(o, s);
  put(o, lowerGeometric);
  put(o, lowerUniform);
  put(o, upperUniform);
  putVector(o, eLo);
  putVector(o, eHi);
  putVector(o, nodes);
  putVector(o, points);
  if (!o) throw std::runtime_error("Write failed: " + path);
}

ManifoldTable ManifoldTable::load(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot read manifold table " + path);
  std::string magic(sizeof kMagic - 1, '\0');
  in.read(magic.data(), static_cast<std::streamsize>(magic.size()));
  if (magic != kMagic) throw std::runtime_error(path + " is not a manifold table of layout version 2.");
  ManifoldTable t;
  t.provenance = getString(in);
  t.fH = getVector<double>(in);
  t.zN = getVector<double>(in);
  t.nEta = get<int>(in);
  t.nLnRho = get<int>(in);
  t.lnRhoMin = get<double>(in);
  t.lnRhoMax = get<double>(in);
  t.referencePressure = get<double>(in);
  t.reactingSpecies.resize(get<std::uint64_t>(in));
  for (auto& s : t.reactingSpecies) s = getString(in);
  t.lowerGeometric = get<int>(in);
  t.lowerUniform = get<int>(in);
  t.upperUniform = get<int>(in);
  t.eLo = getVector<double>(in);
  t.eHi = getVector<double>(in);
  t.nodes = getVector<Node>(in);
  t.points = getVector<double>(in);
  const std::size_t totalNodes = t.fH.size() * t.zN.size() * static_cast<std::size_t>(t.nEta) *
                                  static_cast<std::size_t>(t.nLnRho);
  if (!in || t.fH.empty() || t.zN.empty() || t.nEta < 1 || t.nLnRho < 1 || t.reactingSpecies.empty() ||
      t.eLo.size() != t.fH.size() * t.zN.size() || t.eHi.size() != t.eLo.size() || t.nodes.size() != totalNodes ||
      t.points.size() != totalNodes * static_cast<std::size_t>(t.pointsPerNode()) * t.pointWidth())
    throw std::runtime_error(path + " is truncated or inconsistent.");
  return t;
}

}  // namespace crucible
