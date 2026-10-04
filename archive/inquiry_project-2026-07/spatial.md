# Spatial Material System

The unified geometry and material representation shared by the radiation transport engine and the flight simulator. One source of truth for all physics modules.

---

## OVERVIEW

The ship is represented at two levels:

1. **Component Graph** -- the human-editable definition. A list of components (engines, tanks, shields, cabins) with shapes, materials, positions, and connections. This is what you write to define a ship.

2. **Voxel Octree** -- the derived spatial representation. Generated from the component graph by rasterizing shapes into a 3D grid with adaptive resolution. This is what physics modules query: "what material is at point X?"

The component graph is the source of truth. The octree is rebuilt when the graph changes (ablation, failure, reconfiguration). All physics modules read from the same octree, guaranteeing consistency.

---

## COMPONENT GRAPH

### What a component is

```rust
struct Component {
    id: ComponentId,
    shape: Shape,              // geometric primitive
    position: Vec3,            // meters, ship reference frame
    orientation: Quat,         // rotation from shape's local frame
    material: MaterialId,      // wall/shell material
    wall_thickness: f64,       // meters (for hollow shapes)
    contents: Contents,        // what fills the interior
    
    connections: Vec<Connection>,
    sub_components: Vec<SubComponent>,
    
    // Mutable state
    state: ComponentState,
}

enum Shape {
    Cylinder { radius: f64, length: f64 },
    Cone { radius_base: f64, radius_top: f64, length: f64 },
    Sphere { radius: f64 },
    HemiSphere { radius: f64 },
    Box { width: f64, height: f64, depth: f64 },
    Disk { radius: f64, thickness: f64 },
    Torus { major_radius: f64, minor_radius: f64 },
    Panel { width: f64, height: f64, thickness: f64 },
    Capsule { radius: f64, length: f64 },  // cylinder + hemispheres
}

enum Contents {
    Vacuum,
    Propellant { fuel: MaterialId, fill_fraction: f64 },
    Atmosphere { composition: AtmosphereId, pressure: f64 },
    Solid { material: MaterialId },  // for solid components like truss members
}

struct ComponentState {
    wall_thickness_current: f64,   // may differ from initial (ablation)
    fill_fraction: f64,            // for tanks
    temperature: f64,              // K, from thermal module
    cumulative_dpa: f64,           // average displacement damage
    failed: bool,                  // structural failure flag
}
```

### Connections

```rust
struct Connection {
    target: ComponentId,
    link_type: LinkType,
}

enum LinkType {
    Structural {
        joint: JointType,          // bolted, welded, hinged, etc.
        cross_section: f64,        // m², load-bearing area
    },
    Thermal {
        contact_area: f64,         // m², heat transfer area
    },
    Resource {
        resource: ResourceType,    // LOX, RP1, coolant, air, etc.
        max_flow_rate: f64,        // kg/s
        valve_state: ValveState,   // open, closed, throttled
    },
    Electrical {
        bus: BusId,
        max_power: f64,            // watts
    },
}
```

### Sub-components

Locations within a component that need specific tracking (electronics for SEE, crew for dose, sensors for degradation):

```rust
struct SubComponent {
    id: SubComponentId,
    position: Vec3,              // relative to parent component
    sc_type: SubComponentType,
}

enum SubComponentType {
    Electronics {
        device: DeviceId,        // links to SEE cross-section data
        power_draw: f64,         // watts
        critical: bool,          // failure = cascading consequence
    },
    CrewStation {
        max_occupants: u32,
    },
    Sensor {
        sensor_type: SensorType,
        degradation_model: DegradationId,
    },
    Valve {
        resource: ResourceType,
        state: ValveState,
    },
    Engine {
        propulsion: PropulsionId,  // links to propulsion module
    },
    Reactor {
        reactor_type: ReactorId,   // links to nuclear module
    },
    BlackHoleDrive {
        bh_mass: f64,             // kg, links to hawking module
    },
}
```

---

## MATERIAL DATABASE

### Structure

```rust
struct MaterialDatabase {
    materials: Vec<Material>,
}

struct Material {
    id: MaterialId,
    name: String,
    
    // Composition (for radiation cross-sections)
    elements: Vec<(Element, f64)>,    // (Z, atom fraction)
    nominal_density: f64,              // g/cm³
    
    // Thermal properties
    thermal_conductivity: PropertyFn,  // f(T, dpa) → W/(m·K)
    specific_heat: PropertyFn,         // f(T) → J/(kg·K)
    emissivity: PropertyFn,            // f(T) → dimensionless
    melting_point: f64,                // K
    heat_of_vaporization: f64,         // J/kg
    
    // Structural properties
    youngs_modulus: PropertyFn,        // f(T, dpa) → Pa
    yield_strength: PropertyFn,        // f(T, dpa) → Pa
    ultimate_strength: PropertyFn,     // f(T, dpa) → Pa
    fracture_toughness: PropertyFn,    // f(T, dpa) → Pa·√m
    poisson_ratio: f64,
    thermal_expansion: PropertyFn,     // f(T) → 1/K
    
    // Radiation damage response
    swelling: PropertyFn,              // f(dpa, T) → fractional volume change
    nrt_displacement_energy: f64,      // eV
    sputter_yield_table: TableId,      // from SRIM data
    
    // Radiation interaction tables (from Geant4)
    interaction_tables: TableId,
}

// PropertyFn: either a lookup table or a closure over (T, dpa)
// Stored as polynomial coefficients or piecewise-linear tables
```

### Materials catalog (initial set)

**Structural metals:**
- Al-2219-T87 (primary spacecraft aluminum)
- Al-6061-T6 (general purpose aluminum)
- Ti-6Al-4V (titanium alloy)
- Inconel 718 (high-temperature nickel alloy)
- 304L stainless steel
- Tungsten (shield, high-Z)
- Tungsten-rhenium alloy (high-temperature shield)

**Shielding materials:**
- Polyethylene (HDPE, neutron moderator)
- Borated polyethylene (neutron absorber)
- Lithium hydride (LiH, reactor shadow shield)
- Lead (gamma shield)
- Water (neutron moderator + consumable)

**Composites and thermal protection:**
- Carbon fiber reinforced polymer (CFRP)
- Carbon-carbon composite (high-temperature TPS)
- PICA (phenolic impregnated carbon ablator)
- Silica aerogel (thermal insulation)

**Propellants and fluids:**
- LOX (liquid oxygen)
- RP-1 (kerosene)
- LH2 (liquid hydrogen)
- Hydrazine
- NTO (nitrogen tetroxide)

**Other:**
- Kapton (insulation film)
- Fused silica (optics, windows)
- Silicon (electronics substrate)
- Lunar regolith (for surface habitats)

Each material needs: elemental composition, density, thermal properties, structural properties, and Geant4 interaction tables. Thermal and structural properties come from engineering handbooks (ASM, MPDB). Radiation damage correlations from the irradiation effects literature.

---

## VOXEL OCTREE

### Data structure

```rust
struct Octree {
    root: OctreeNode,
    bounds: AABB,           // bounding box of the entire ship
    max_depth: u32,         // limits minimum voxel size
    min_voxel_size: f64,    // meters (typically ~1mm)
}

enum OctreeNode {
    // Leaf: uniform material in this region
    Leaf {
        voxel: Voxel,
        bounds: AABB,
    },
    // Branch: subdivided into 8 children
    Branch {
        children: Box<[OctreeNode; 8]>,
        bounds: AABB,
    },
}

struct Voxel {
    material_id: u16,      // 0 = vacuum
    density_frac: f16,     // fraction of nominal density
    damage_dpa: f16,       // cumulative displacement damage
}
// 6 bytes per leaf voxel
```

### Construction from component graph

```rust
fn build_octree(components: &[Component], max_depth: u32) -> Octree {
    let bounds = compute_bounding_box(components);
    let root = build_node(bounds, components, 0, max_depth);
    Octree { root, bounds, max_depth, min_voxel_size: bounds.size() / 2^max_depth }
}

fn build_node(bounds: AABB, components: &[Component], depth: u32, max_depth: u32) -> OctreeNode {
    // Which components intersect this region?
    let overlapping: Vec<_> = components.iter()
        .filter(|c| c.shape_bounds().intersects(&bounds))
        .collect();
    
    if overlapping.is_empty() {
        return Leaf { voxel: Voxel::vacuum(), bounds };
    }
    
    // Is this region uniform? (entirely inside one component's wall, or entirely vacuum)
    if let Some(uniform_material) = check_uniform(&bounds, &overlapping) {
        return Leaf { voxel: Voxel::new(uniform_material, 1.0, 0.0), bounds };
    }
    
    // Need to subdivide (material boundary crosses this region)
    if depth >= max_depth {
        // At max depth: pick the dominant material
        return Leaf { voxel: dominant_material(&bounds, &overlapping), bounds };
    }
    
    // Subdivide into 8 children
    let children = bounds.octants().map(|child_bounds| {
        build_node(child_bounds, &overlapping, depth + 1, max_depth)
    });
    
    Branch { children: Box::new(children), bounds }
}
```

### Key operations

**Point query: what material is at position P?**
```rust
fn lookup(&self, point: Vec3) -> &Voxel  // O(log N), tree descent
```

**Ray march: step through voxels along a ray**
```rust
fn ray_march(&self, origin: Vec3, direction: Vec3) -> RayMarchIterator
// Yields: (voxel, entry_point, exit_point, distance_through_voxel)
```

**Distance to exit: how far until the ray leaves this voxel?**
```rust
fn distance_to_exit(&self, position: Vec3, direction: Vec3) -> f64  // O(1), ray-AABB
```

**Mutate: modify a single voxel**
```rust
fn set_voxel(&mut self, position: Vec3, voxel: Voxel)  // O(log N)
// May need to split a leaf node if the modification creates non-uniformity
```

**Local rebuild: regenerate a region from the component graph**
```rust
fn rebuild_region(&mut self, bounds: AABB, components: &[Component])
// Used after component state changes (ablation, failure)
// Only rebuilds the affected region, not the entire tree
```

**Volume query: find all voxels of a given material**
```rust
fn find_material(&self, material_id: MaterialId) -> Vec<(Vec3, &Voxel)>
// For tallying total mass, computing center of mass, etc.
```

**Enclosed volume detection: find sealed cavities**
```rust
fn find_enclosed_volumes(&self) -> Vec<EnclosedVolume>
// Flood-fill from vacuum voxels to identify sealed pressurized volumes
// Used by atmosphere module to determine cabin volume and detect leaks
```

### Memory estimates

| Ship | Dimensions | Estimated leaves | Memory |
|---|---|---|---|
| Falcon 9 | 70m × 3.7m | ~20M | ~120 MB |
| NTP Mars ship | 50m × 8m | ~15M | ~90 MB |
| BH drive starship | 100m × 15m | ~30M | ~180 MB |

These are estimates. The actual count depends on geometric complexity (more complex = more boundary refinement = more leaves). The octree compresses well: 99%+ of a rocket's volume is either vacuum or uniform propellant, represented as large leaf nodes at coarse depth.

---

## SURFACE MESH (for aerodynamics)

The octree's voxelized surface has staircase artifacts. For aerodynamic computation, the external surface is extracted directly from the component graph's geometric primitives (cylinders, cones, etc.), NOT from the octree. These are smooth analytic surfaces.

```rust
struct SurfaceMesh {
    panels: Vec<SurfacePanel>,
}

struct SurfacePanel {
    component_id: ComponentId,
    surface_type: SurfaceType,  // cylinder, cone, sphere, flat
    area: f64,                  // m²
    normal: Vec3,               // outward surface normal (or function for curved surfaces)
    material: MaterialId,       // surface material (for emissivity, ablation, aerodynamic heating)
}
```

Aerodynamic coefficients are precomputed by the characterization tool from these smooth shapes, NOT from the octree. The octree is for internal physics (radiation, thermal conduction, structural loads). The surface mesh is for external physics (drag, heating, radiation to space).

---

## SHIP DEFINITIONS

Each ship is a file defining its component graph. Examples:

### Falcon 9 / Dragon (simplified)

```
components:
  # First stage
  - id: s1_tank_lox
    shape: Cylinder(r=1.83, l=16.0)
    position: (0, 0, 20.0)
    material: al-2219-t87
    wall_thickness: 4mm
    contents: LOX(fill=0.95)
    
  - id: s1_tank_rp1
    shape: Cylinder(r=1.83, l=8.5)
    position: (0, 0, 4.25)
    material: al-2219-t87
    wall_thickness: 4mm
    contents: RP1(fill=0.95)
    
  - id: s1_engine_section
    shape: Cone(r_base=1.83, r_top=2.6, l=3.0)
    position: (0, 0, -1.5)
    material: inconel-718
    wall_thickness: 3mm
    sub_components:
      - Engine(id: merlin_1, type: merlin_1d, position: (0, 0, -2.5))
      # ... 8 more Merlins
    
  - id: interstage
    shape: Cylinder(r=1.83, l=4.0)
    position: (0, 0, 30.0)
    material: cfrp
    wall_thickness: 2mm
    contents: Vacuum

  # Second stage (similar pattern)
  # ...

  # Dragon capsule
  - id: dragon_cabin
    shape: Capsule(r=1.85, l=3.0)
    position: (0, 0, 42.0)
    material: al-2219-t87
    wall_thickness: 8mm
    contents: Atmosphere(cabin_air, 101.3kPa)
    sub_components:
      - CrewStation(id: seat_1, pos: (0.3, 0, 42.5))
      - CrewStation(id: seat_2, pos: (-0.3, 0, 42.5))
      - Electronics(id: flight_comp, pos: (0, 0.5, 43.0), device: dragon_fc)
      
  - id: dragon_heatshield
    shape: Disk(r=1.85, thickness=0.08)
    position: (0, 0, 40.4)
    material: pica
    
connections:
  - Structural(s1_tank_lox, interstage, bolted_ring)
  - Structural(s1_tank_rp1, s1_engine_section, welded)
  - Resource(s1_tank_lox, s1_engine_section, LOX, 2500 kg/s)
  - Resource(s1_tank_rp1, s1_engine_section, RP1, 1100 kg/s)
  # ...
```

This definition is human-readable and editable. The octree, thermal network, structural model, and resource graph are all derived from it automatically.

---

## USAGE BY PHYSICS MODULES

| Module | Uses | How |
|---|---|---|
| **Radiation transport** | Octree + material database | Ray-march through octree, sample interactions from material tables |
| **Thermal** | Component graph + material database | Thermal network from components + connections. k(T, dpa) from material database. |
| **Structural** | Component graph + material database | Stiffness from component shapes + materials. Loads from connections. σ_y(T, dpa) from database. |
| **Aerodynamics** | Surface mesh | Drag and heating from smooth analytic surfaces (characterization tables). |
| **Resources** | Component graph | Tank volumes, fill fractions, pipe connections, valve states. |
| **Atmosphere** | Octree (enclosed volumes) | Flood-fill to find sealed cavities. Volume from voxel count. Leaks from breached walls. |
| **Pressure** | Component graph + octree | Tank pressure from fill fraction + temperature. Cabin pressure from atmosphere module. Leak rate from hole size (detected via octree damage). |
| **Electrical** | Component graph | Power bus topology. Generator → bus → load connections. |

All modules read from the same material database for property evaluation. When any module changes the ship state (thermal changes T, radiation changes dpa, structural changes failed, resources change fill), the affected derived representations update.
