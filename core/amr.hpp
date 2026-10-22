#pragma once
#include "core/flow.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

namespace crucible {

// Block-structured adaptive resolution in the index space of the body-fitted mapping
// (docs/evidence/AMR.md): levels 0..L, level l a Flow on (nz 2^l) x (nr 2^l) cells spanning the whole
// domain, ratio 2 in space and time, Berger and Colella's recursive subcycling with refluxing. The
// finest mesh is built from the definition and every coarser one by Mesh::coarsened, so a coarse cell
// is the exact union of its children.
// Step B2 (AMR.md, Step B2: setup): a level above the base runs its kernels over the cells within the
// ghost width of its active cells (Flow::setRegions), fills the children of the parents of its ghosts
// from the level below before each stage, and keeps its other cells (off cells) fixed. Step B1, every
// kernel over the whole level and every non-active cell filled, stays selectable (setRestricted) and
// gives the same active cells bit for bit. The registers carry the convective and viscous fluxes; the
// wall functions and the SST property corrections are refused. Step C (setTags) moves the levels'
// active cells with the flow; without it the pattern is setActive's.
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
    // Ghost width in cells of the finer level: 3 with turbulence (the Strang-split half-source at a
    // distance-2 ghost needs that ghost's gradients), else 2 (the flux stencils). It sets the nesting
    // check and, in B2, the regions; B2 refuses less than 3 with turbulence.
    void setGhostWidth(int width) { ghost_=width;checked_=false; }
    // B2 (true, the default) or B1 (false).
    void setRestricted(bool restricted) { restricted_=restricted;checked_=false; }
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
    // Called after each step of a level (its index and dt), before the levels above it step: e.g. a
    // reaction step, which must keep each cell's mass and energy.
    void setAfterStep(std::function<void(int, double)> afterStep) { afterStep_=std::move(afterStep); }
    // Step C (AMR.md, Step C: setup): tags on each level below the top decide the level above's active
    // cells, before the first base step (one pass per level above the base, so the initial levels reach the
    // top) and after every interval-th accepted base step. A threshold of zero or less turns its tag off;
    // empty hydrogen weights turn the f_H tag off.
    struct Tags {
        int interval{2}, buffer{2};
        double throatHalfWidth{0.5};  // in throat radii, about the wall's narrowest interior station
        int throatLevel{2};           // the throat and the wall strip refine to this level
        double wallYPlus{300};        // the wall row's target first-cell y+ (laminar estimate)
        double mixture{0.02}, temperature{0.05}, pressure{0.10}, vorticity{0.05};
        std::vector<double> hydrogen;  // per species, the H element's mass over the species' mass
        std::function<bool(int, int, int)> custom;  // a caller's own tag on level l's cell (i, j), e.g. a test's
    };
    void setTags(Tags tags);
    // Regrids that changed the pattern, and the largest content error of the children against their parent in
    // the last such regrid's fill (mass, momentum, energy and partial densities): over the parent's content
    // (momentum's magnitude for both components), and over the children's gross content (the sum of |content|).
    [[nodiscard]] std::uint64_t regrids() const { return regrids_; }
    [[nodiscard]] double fillError() const { return fillError_; }
    [[nodiscard]] double fillGross() const { return fillGross_; }
    // The filled parents' mass change (children's content less the parent's), signed and summed over every regrid.
    [[nodiscard]] double fillMassChange() const { return fillMassChange_+fillMassCompensation_; }
    // Level l's integrated mass exchange over its owned boundaries (compensated), as the budget takes it.
    [[nodiscard]] double massExchange(int l) const;
    // The composite budget (AMR.md, criterion 2 (b)): the base's content, which holds the composite
    // after average-down and reflux, less its content at the first step and every level's exchanges
    // booked at its owned cells since (boundary faces, igniter, corrections, spray). The residuals are
    // absolute; momentumGross and igniterEnergy are the levels' sums.
    struct Budget {
        double mass{}, energy{}, momentum{};
        double massResidual{}, energyResidual{}, momentumResidual{}, momentumGross{}, igniterEnergy{};
    };
    [[nodiscard]] Budget budget() const;

private:
    // A level's state saved for a redone base step.
    struct Saved {
        std::vector<Conserved> state;
        std::vector<double> species, turbulence, temperature;
        double time{}, dt{}, heat{}, heatCompensation{}, mass{}, energy{}, correction{}, momentum{}, gross{}, clipped{}, prescribed{};
        double massCompensation{}, energyCompensation{}, correctionCompensation{}, momentumCompensation{};
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
        // B2, levels above the base, in the level below's cells: the parents of the cells within the
        // ghost width that are not active, and those with their face neighbours (whose values the
        // fill's gradients read).
        std::vector<std::size_t> parents, around;
        Saved saved;
    };
    std::vector<Level> levels_;
    int ghost_{2}, top_{0};
    bool checked_{false}, started_{false}, restricted_{true};
    std::function<void(int, double)> afterStep_;
    double time_{};
    std::uint64_t steps_{}, rejected_{};
    void check();
    bool advance(int l, double t, double dt, int half);
    // Fill the children of level l's parents (B2) or every cell that is not active (B1) from level l - 1
    // at the fraction theta of l - 1's current step: linear in time between its old and new states,
    // conservative linear in space; each filled cell's temperature guess becomes its parent's.
    void fill(int l, std::vector<Conserved>& state, std::vector<double>& species, std::vector<double>& turbulence, double theta);
    // The fill's spatial part: level l - 1's parent p's values (Level::values, at p and its face
    // neighbours) onto its four children on level l.
    void spread(int l, std::size_t p, std::vector<Conserved>& state, std::vector<double>& species, std::vector<double>& turbulence);
    // Step C: new active sets from the tags, top down; the cells that become active filled, bottom up.
    // initial: before the first base step, when the whole wall row is tagged.
    void regrid(bool initial);
    Tags tags_;
    bool tagged_{false}, regridded_{false};
    std::uint64_t regrids_{};
    double fillError_{}, fillGross_{}, fillMassChange_{}, fillMassCompensation_{};
    void averageDown(int l);
    bool reflux(int l);
    void save();
    void restore();
    Budget initial_{};
    [[nodiscard]] Budget content() const;
};

}  // namespace crucible
