# Terminal Interface

The simulator's user interface. Terminal-based with a command line, telemetry display, window system, and manual flight mode. Scriptable.

---

## MODES

### Command mode (default)

Text prompt. Type commands, hit enter. The simulation runs in the background (or paused, player's choice).

```
> throttle 0.8
Throttle set to 80%.

> target mars
Target set: Mars. Current distance: 78.3 million km. Hohmann ΔV: 3.6 km/s.

> timewarp 1000
Time acceleration: 1000x. Coordinate time advancing at 1000 s/s.

> query stars 50pc filter G-type
Querying G-type stars within 50 pc... 847 results.
[list of stars with ID, type, distance, luminosity]

> inspect GCS-R8.247-P000.12-Z+019-0003
Star: GCS-R8.247-P000.12-Z+019-0003
  Type: G2V, Mass: 1.02 M_sun, Age: 4.8 Gyr, Distance: 12.3 pc
  [generating system...]
  Planets: 4 (1 gas giant, 2 sub-Neptunes, 1 super-Earth)
  Planet c: 1.4 R_earth, 3.1 M_earth, in habitable zone
  [inspect planet c for detail]
```

### Flight mode (toggle with `flight` command)

Real-time control via keypresses. Terminal switches to raw input mode (crossterm/ncurses).

```
W/S: increase/decrease throttle
A/D: yaw left/right
↑/↓: pitch up/down
Q/E: roll left/right
Space: kill rotation (point-and-hold)
Tab: cycle target
Esc: return to command mode
```

During flight mode, a HUD is displayed with critical telemetry (velocity, altitude, fuel, attitude indicators).

### Script mode

Execute a sequence of commands from a file:

```
> run scripts/hohmann_mars.script
```

Scripts are text files with one command per line. Support variables, conditionals, and wait-for-condition:

```
# hohmann_mars.script
set target mars
wait altitude > 200km
attitude prograde
throttle 1.0
wait deltav > 3600
throttle 0
timewarp 10000
wait distance_to target < 1000000km
timewarp 1
attitude retrograde
throttle 1.0
wait velocity_relative target < 100
throttle 0
print "Mars orbit insertion complete"
```

---

## TELEMETRY DISPLAY

The default view shows critical ship data, updated every display frame (~30 Hz):

```
╔══════════════════════════════════════════════════════════════════╗
║ INQUIRY FLIGHT SIMULATOR          T: 2000-03-15 14:32:07 TDB   ║
║ Ship: Falcon 9 / Dragon           τ: 0d 00h 14m 32.07s        ║
╠══════════════════════════════════════════════════════════════════╣
║ TRAJECTORY                        │ PROPULSION                  ║
║  Alt:  187.3 km (LEO)             │  Throttle: 0% (coast)       ║
║  Vel:  7,784 m/s                  │  Isp: -- (engine off)       ║
║  Acc:  0.00 m/s² (freefall)       │  Fuel: 1,247 kg (S2 RP-1)  ║
║  Inc:  28.5°                      │  Oxid: 3,124 kg (S2 LOX)   ║
║  Ecc:  0.0012                     │  ΔV remaining: 1,847 m/s    ║
║  Per:  186.1 km                   │                             ║
║  Apo:  189.4 km                   │                             ║
╠══════════════════════════════════════════════════════════════════╣
║ ENVIRONMENT                       │ SHIP STATUS                 ║
║  Gravity: 9.34 m/s² (Earth)       │  Hull temp: 293 K (nominal) ║
║  Atm density: 2.1e-10 kg/m³      │  Cabin: 101.3 kPa, 22°C    ║
║  Solar flux: 1361 W/m²            │  O2: 20.9%  CO2: 0.04%     ║
║  Radiation: 0.82 μSv/hr          │  Power: 1.2 kW / 1.5 kW    ║
║  Nearest body: Earth (187 km)     │  Battery: 94%               ║
║  Timewarp: 1x                     │  Crew dose: 0.12 mSv total  ║
╠══════════════════════════════════════════════════════════════════╣
║ > _                                                             ║
╚══════════════════════════════════════════════════════════════════╝
```

## WINDOWS

Open additional views alongside the main telemetry:

```
> window map             # Galaxy/system map with ship position
> window orbit           # Orbital elements and trajectory prediction
> window thermal         # Component temperature map
> window radiation       # Dose rates and shielding status
> window systems         # All subsystem status (resources, electrical, atmosphere)
> window structural      # Load/stress indicators per joint
> window bodies          # List of nearby bodies with data
```

Windows are terminal panes (split view via crossterm). Each window updates at its own rate (map: 1 Hz, telemetry: 30 Hz, thermal: 1 Hz).

## CORE COMMANDS

### Navigation
- `target <body/star/ID>` -- set navigation target
- `attitude <prograde/retrograde/normal/antinormal/radial/target/sun/custom(x,y,z)>` -- set attitude hold
- `maneuver <prograde/retrograde/normal/...> <dv> [at <time>]` -- plan or execute a maneuver

(Control laws, actuator allocation, the burn executor, and the display-path trajectory
predictor behind these commands are specified in ship.md Section 9. Scripted open-loop
guidance profiles — pitch programs, entry schedules — are ship.md Section 9 + script mode.)

### Propulsion
- `throttle <0-1 | percent>` -- set engine throttle
- `engine <on/off> [engine_id]` -- toggle specific engine
- `reactor <startup/shutdown/set_power percent>` -- reactor control
- `bh_drive feed_rate <kg/s>` -- BH drive feed control (there is NO off switch: Hawking
  emission cannot be stopped; thrust changes only via BH mass drift — see hawking.md)

### Time
- `timewarp <multiplier>` -- set time acceleration (1, 10, 100, 1000, ...)
- `timewarp off` -- return to 1x
- `pause` / `resume` -- pause/resume simulation
- `advance <duration>` -- advance by specific time (e.g., `advance 6h`)

### Query
- `query stars <radius> [filter ...]` -- find stars within radius
- `query bodies <radius> [filter ...]` -- find all bodies (stars, planets, asteroids)
- `query system [star_id]` -- list bodies in a star system
- `inspect <body_id>` -- detailed information about a body
- `bookmark <name>` -- save current target
- `goto <bookmark>` -- recall a bookmark

Filters: `filter G-type`, `filter habitable`, `filter has-life`, `filter mass>10`, `filter binary`, `filter neutron-star`

### Ship management
- `status` -- full ship status
- `status <subsystem>` -- specific subsystem (thermal, structural, electrical, ...)
- `valve <valve_id> <open/close/throttle percent>` -- control valves
- `systems <system> <on/off>` -- toggle subsystem (heaters, sensors, comms)

### Data
- `log start [filename]` -- start recording telemetry to file
- `log stop` -- stop recording
- `export trajectory [filename]` -- export trajectory data
- `export system [star_id] [filename]` -- export system data

### Meta
- `save [filename]` -- save simulation state
- `load <filename>` -- load simulation state
- `ship <ship_definition>` -- switch to a different ship (loads new tables)
- `help [command]` -- help
- `quit` -- exit

---

## ALERTS AND WARNINGS

Critical events produce alerts that appear in the telemetry header:

```
⚠ FUEL LOW: S2 RP-1 at 5% (124 kg remaining)
⚠ O2 DROPPING: 18.2% (crew at risk below 16%)
🔴 STRUCTURAL: Joint S1-interstage stress at 94% of yield
🔴 RADIATION: SPE detected. Dose rate: 15.2 mSv/hr. Seek shelter.
⚠ TIMEWARP PAUSED: Ship state at lookup table boundary
```

Alerts are logged and can be reviewed with `alerts` command.

---

## IMPLEMENTATION

| Component | Rust crate | Purpose |
|---|---|---|
| **crossterm** | Terminal raw mode, colors, cursor control | Input handling, display rendering |
| **tui-rs / ratatui** | Terminal UI framework | Window layout, widgets, charts |
| Custom | Script parser + interpreter | Scripting engine |
| Custom | Command parser | Command-line parsing with tab completion |
