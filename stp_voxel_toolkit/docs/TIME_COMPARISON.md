# Time Comparison Report

Benchmark input:

`D:\File\Project\Section\体素化代码\test1.stp`

Benchmark grid:

`247 x 400 x 97 = 9,583,600` total cells

STL:

`40,574` triangles, reused for both runs.

## C++ Voxelizer

| Stage | Seconds | Percent | Cumulative | Detail |
|---|---:|---:|---:|---|
| STEP -> STL | 0.000 | 0.0% | 0.000 | reused test1.stl |
| Read STL / Grid | 0.327 | 22.0% | 0.327 | 40,574 triangles, grid 247x400x97=9,583,600 |
| Voxelize | 1.140 | 76.6% | 1.467 | C++, occupied 110,801, fill 0.01156152 |
| Write Occupancy VTK | 0.021 | 1.4% | 1.488 | 9.14 MB |
| Total | 1.488 | 100.0% | 1.488 | all stages |

## Python Voxelizer

| Stage | Seconds | Percent | Cumulative | Detail |
|---|---:|---:|---:|---|
| STEP -> STL | 0.000 | 0.0% | 0.000 | reused test1.stl |
| Read STL / Grid | 0.334 | 2.8% | 0.334 | 40,574 triangles, grid 247x400x97=9,583,600 |
| Voxelize | 11.682 | 97.0% | 12.016 | Python, occupied 110,801, fill 0.01156152 |
| Write Occupancy VTK | 0.022 | 0.2% | 12.038 | 9.14 MB |
| Total | 12.039 | 100.0% | 12.039 | all stages |

## Summary

The voxelization kernel is the bottleneck.

- C++ voxelization: `1.140 s`
- Python voxelization: `11.682 s`
- Speedup: `10.25x`

For high-resolution runs, use the default C++ path.

## Full Pipeline Example

A cold run that includes STEP -> STL at `resolution=160` produced:

| Stage | Seconds | Percent | Detail |
|---|---:|---:|---|
| STEP -> STL | 7.635 | 92.2% | 40,574 triangles, 8.78 MB |
| Read STL / Grid | 0.331 | 4.0% | 99x160x39=617,760 |
| Voxelize | 0.300 | 3.6% | C++, occupied 7,602 |
| Write Occupancy VTK | 0.001 | 0.0% | 0.59 MB |
| Extract Surface VTK | 0.014 | 0.2% | 29,672 quad faces |
| Total | 8.282 | 100.0% | all stages |

When STL is not reused, CAD triangulation dominates runtime. When STL is reused, voxelization dominates runtime.
