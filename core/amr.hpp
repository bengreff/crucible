#pragma once
#include "core/flow.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace crucible {

// Block-structured adaptive resolution in the index space of the body-fitted mapping
// (docs/evidence/AMR.md): levels 0..L, level l a Flow on (nz 2^l) x (nr 2^l) cells spanning the whole
// domain, ratio 2 in space and time, Berger and Colella's recursive subcycling with refluxing. The
// finest mesh is built from the definition and every coarser one by Mesh::coarsened, so a coarse cell
// is the exact union of its children.
// Step B1 (AMR.md, Deviations): every level's kernels run over the whole level and every cell that a
// level does not update is filled from the level below before each stage; inviscid only.
class Amr {
public:
    // base: the level-0 definition; finer: the number of levels above it (1 to 3).
    Amr(const Definition& base, int finer);
    [[nodiscard]] int levels() const { return static_cast<int>(levels_.size()); }
    Flow& level(int l) { return *levels_[static_cast<std::size_t>(l)].flow; }
    [[nodiscard]] const Flow& level(int l) const { return *levels_[static_cast<std::size_t>(l)].flow; }
    void setThreads(int threads);
    // Level l's (1..L) active cells in its own index space (1 active, 0 not), whole parent cells. The
    // nesting (each active cell's ghost-width box over active parents) is checked at the next step.
    void setActive(int l, std::vector<unsigned char> active);
    [[nodiscard]] const std::vector<unsigned char>& active(int l) const { return levels_[static_cast<std::size_t>(l)].flow->active_; }
    // Whether level-l cell q is covered by level l + 1's active cells.
    [[nodiscard]] bool covered(int l, std::size_t q) const;
    // Ghost width in cells of the finer level (2, the widest stencil). In B1 every cell that is not
    // active is filled, so the width enters only the nesting check.
    void setGhostWidth(int width) { ghost_=width;checked_=false; }
    // Replace each covered cell by the volume-weighted mean of its children, finest level first (after
    // setting the levels' initial states).
    void averageDown();
    // One base step: the smallest of 2^l times each level's stable step over its active cells, ended
    // on the base's schedule events; on a non-admissible stage anywhere the whole base step is redone
    // from the saved levels at half the size.
    double step(double maxDt = std::numeric_limits<double>::infinity());
    void advanceTo(double time);
    [[nodiscard]] double time() const { return time_; }
    [[nodiscard]] std::uint64_t steps() const { return steps_; }
    [[nodiscard]] std::uint64_t rejectedSteps() const { return rejected_; }

private:
    // A level's state saved for a redone base step.
    struct Saved {
        std::vector<Conserved> state;
        std::vector<double> species, turbulence, temperature;
        double time{}, dt{}, heat{}, mass{}, energy{}, correction{}, momentum{}, gross{}, clipped{}, prescribed{};
        std::uint64_t steps{};
    };
    struct Level {
        std::unique_ptr<Flow> flow;
        // The state at the start of the level's current step, for the ghosts of the level above.
        std::vector<Conserved> oldState;
        std::vector<double> oldSpecies, oldTurbulence;
        // The level's values at the fill's time, per cell in stride 4 + ns + nt (scratch of fill()).
        std::vector<double> values;
        // Levels above the base: the face integrals of its steps since the level below last stepped.
        std::vector<double> fineSum;
        Saved saved;
    };
    std::vector<Level> levels_;
    int ghost_{2}, top_{0};
    bool checked_{false};
    double time_{};
    std::uint64_t steps_{}, rejected_{};
    void check();
    bool advance(int l, double t, double dt, int half);
    // Fill every cell of level l that is not active from level l - 1 at the fraction theta of l - 1's
    // current step: linear in time between its old and new states, conservative linear in space.
    void fill(int l, std::vector<Conserved>& state, std::vector<double>& species, std::vector<double>& turbulence, double theta);
    void averageDown(int l);
    bool reflux(int l);
    void save();
    void restore();
};

}  // namespace crucible
