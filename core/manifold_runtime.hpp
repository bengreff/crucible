#pragma once
#include <array>
#include <cstddef>
#include <vector>

#include "core/manifold_table.hpp"
#include "core/medium.hpp"

namespace crucible {

// Table B's runtime for a table built with the specific entropy (candidate 3): the cell placed on the 16
// corners' blended curve by its own moles per kilogram, and its reaction advanced along that curve
// (docs/evidence/TABLE_B.md, "The cell placed by its moles per kilogram" and "The runtime's step on the placed
// point"). No Cantera here: the species' element fractions come from the caller. A cell the table does not
// handle returns its cause, and the caller integrates it directly (CVODES), so the table's fallback is exact.
class ManifoldRuntime {
 public:
  // The element mass fractions of one medium species (kg of H, O and N per kg of the species).
  struct Elements { double h{}, o{}, n{}; };

  // Why a cell is not on the table (Table: it is). The placement's causes, then the step's.
  enum class Cause : int {
    Table,
    NoHydrogenOrOxygen,
    OutsideAxes,
    Argon,             // AR above 1e-12 (the table's nodes carry none)
    CornerRefused,     // a corner the builder did not store
    CornerFlagged,     // a corner flagged outside the thermo range, or degenerate (f_H 0 or 1)
    NoFallInN,
    FrozenReach,       // N between the two branches' ends: a frozen corner's gap
    BelowNCp,          // N below the blended complete-combustion N
    OffCurveAboveSEq,  // w < 1 and s* above the blended s_eq
    OffCurveBelowStart,
    DYAboveMax,        // the placed point's max |dY| above 0.01
    Induction,         // on the lower branch below phi_N 0.3: its induction part, where a mixed cell's own
                       // intermediates set the rate and the curve does not carry them
    FrozenCorner,      // a corner of nonzero weight whose placed branch is frozen
    NoRate,            // a corner of nonzero span whose lambda at the cell's u is not positive and finite
    Count
  };
  static const char* causeName(Cause c);

  struct Placement {
    Cause cause{Cause::Table};
    int branch{-1};       // 0 lower, 1 upper
    double phi{}, w{}, tau{}, u{}, dY{};
    int iterations{};     // Illinois iterations for tau_N
    int substeps{};       // the step's sub-steps (advance only)
    // FrozenCorner: the heaviest frozen corner of nonzero weight on the placed branch, its weight and node.
    double frozenW{};
    std::size_t frozenNode{};
  };

  // Per-thread scratch.
  struct Work {
    std::vector<double> yc, yb, yr, products;
  };

  // medium: the flow's species in order; elements: one per medium species. Every species must be one of the
  // table's reacting species, N2, or AR.
  ManifoldRuntime(const ManifoldTable& table, const Medium& medium, std::vector<Elements> elements);

  [[nodiscard]] Work work() const;

  // Places a cell (its internal energy e [J/kg], rho, T and every medium species' Y) on the blended curve.
  Placement place(double e, double rho, double t, const double* y, Work& work) const;
  // Places it and, if it is on the table, advances its reaction over dt at fixed e and rho: y returns the new
  // composition (the reacting species from the blend at the new u; N2 and AR unchanged). y is untouched
  // otherwise.
  Placement advance(double e, double rho, double t, double* y, double dt, Work& work) const;

  // The blended curve's reacting-species Y at u on a branch (0 lower, 1 upper), at the cell's corners; false if the
  // cell has no usable corners. For checks.
  bool curve(double e, double rho, const double* y, int branch, double u, double* out, Work& work) const;

  // The mixture's specific entropy [J/(kg K)] on the medium's NASA-7 data at reference pressure pRef, with
  // the constant-cp continuation above each species' tHigh: s = sum_k Y_k/W_k (R s°_k - R ln(X_k p / p°)).
  static double entropy(const Medium& medium, double t, double rho, const double* y, double pRef);

  // The cell's Table B coordinates: f_H, Z_N and eta (against its own complete-combustion products at 200 K
  // and 6000 K). f_H is NaN with no H or O.
  struct Coordinates { double fH{}, zN{}, eta{}; };
  Coordinates coordinates(double e, const double* y, Work& work) const;

  [[nodiscard]] const ManifoldTable& table() const { return table_; }

 private:
  const ManifoldTable& table_;
  const Medium& medium_;
  std::vector<Elements> elements_;
  std::vector<std::size_t> reacting_;  // medium index of each table reacting species
  std::vector<double> reactingMass_;   // their molar masses
  std::size_t n2_, ar_, h2o_, h2_, o2_;
  double hInWater_{};
  std::vector<double> pointN_;         // N of every stored point's reacting species
  std::vector<double> baseU_[2];       // each branch's base u, decreasing from 1 to 0
};

}  // namespace crucible
