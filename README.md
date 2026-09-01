# Kryonomos.jl

*The law of the cold: coupled ice-sheet modeling.*

Worked examples for running the [Yelmo.jl](https://github.com/fesmc/Yelmo.jl) ice-sheet
model coupled to the systems that drive it — solid-Earth deformation and subglacial
hydrology — together with the glue code those couplings need.

The ice sheet is the core. Each example wires one more system to it and shows the whole
loop: build the model, push the coupled fields in, step, pull the response back out.

Both Yelmo backends appear throughout:

- **`YelmoModel`** — the native Julia solver.
- **`YelmoMirror`** — the same physics, namelists and input data, `ccall`ing into the
  Fortran `libyelmo_c_api.so` instead of running Julia code.

Coupled fields cross the boundary the same way in both cases, so an example written
against one backend reads almost identically against the other.

## Before you run

The examples resolve input data and Yelmo's Fortran namelists relative to `examples/`.
Create both links yourself — neither is repository content:

```bash
cd examples
ln -s /path/to/ice_data    ice_data
ln -s /path/to/yelmo/input input
```

Then run any example **from `examples/`**, not from its own folder:

```bash
julia --project=. fasthydrology/Greenland_yelmo-fasthydrology.jl
```

Each script `cd`s to `examples/` itself, so it works from anywhere — but the project
environment and the Fortran side's relative paths both live there.

## Solid Earth — `examples/fastisostasy/`

Bedrock responding to the ice load through
[FastIsostasy.jl](https://github.com/JanJereczek/FastIsostasy.jl).

| Script | What it shows |
|---|---|
| `Greenland_yelmo-fastisostasy.jl` | Greenland via the Fortran-backed `YelmoMirror` |
| `Greenland_yelmo-fastisostasy-native.jl` | The same coupling against the native Julia backend |

The two are deliberately close to each other: reading them side by side is the quickest
way to see what changes between backends and what doesn't.

## Subglacial hydrology — `examples/fasthydrology/`

Effective pressure and basal water from
[FastHydrology.jl](https://github.com/fesmc/FastHydrology.jl), fed back into Yelmo's basal
friction each step.

| Script | What it shows |
|---|---|
| `Greenland_yelmo-fasthydrology.jl` | Greenland, 106×181 @ 16 km |
| `Antarctica_yelmo-fasthydrology.jl` | Antarctica, 191×191 @ 32 km |

Both select the hydrology model at the top — Kazmierczak et al. (2024), height-above-
buoyancy, or Shakti — and are structured identically so the two domains diff cleanly.

They are **coupling-mechanics demonstrations, not physically forced runs**: no `smb`,
`T_srf` or `Q_geo` forcing is applied, and `Q_geo` falls back to a constant. Add real
forcing before reading anything physical into the output.

## Layout

```
examples/
  Project.toml                  shared environment for every example
  fastisostasy/                 Yelmo + FastIsostasy
  fasthydrology/                Yelmo + FastHydrology
  run01/                        namelists written by the FastIsostasy examples
```
