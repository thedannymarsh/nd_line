# nd_line

[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/thedannymarsh/nd_line)
[![PyPI](https://img.shields.io/pypi/v/nd_line.svg)](https://pypi.org/project/nd_line/)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/nd_line?period=month&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads/month)](https://pepy.tech/projects/nd_line)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/nd_line?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/nd_line)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/nd_line.svg)](https://pypi.org/project/nd_line/)

Interpolate points on an n-dimensional line by Euclidean arc length.

### Installation

`pip install nd_line`

#### Classes

- `nd_line`: Linear representation of ordered n-dimensional points.
- `nd_spline`: B-Spline representation of ordered n-dimensional points. Child class of nd_line.

#### Methods

- `ln.interp(dist)`: returns the point `dist` along the arc of the line

- `ln.interp_rat(ratio)`: `ratio` is between 0 and 1; returns the point at `ratio * length`

- `ln.to_spline(samples)`: returns a **new** `nd_spline` sampled from a spline fit (does not mutate in place); `samples` is the number of points on the new line (defaults to the original count)

- `ln.dist_from(point)`: returns the distance from all line points to the given point

- `ln.closest_idx(point)`: returns the index of the closest point

- `ln.within_idx(point, dist)`: returns the indices of all points within `dist` of the given point

- `ln.contig_within_idx(point, dist)`: same as `within_idx`, but indices are consecutive and include the closest point

- `ln.resample(dists)`: returns a new `nd_line` with points interpolated at the given distances

- `nd_spline.recursive_upsample(tol)`: returns a new `nd_spline` upsampled until the Euclidean arc length is within `tol`

- `nd_spline.resample(dists)`: same as `ln.resample`, evaluated on the spline

#### Attributes

- `ln.points`: the points of the line
- `ln.length`: the length of the line
- `ln.lengths`: Euclidean length of each segment (array of len(points) - 1)
- `ln.cumul`: cumulative arc length at each point (starts at 0.0, ends at ln.length)
- `ln.type`: `'linear'` unless created by `to_spline` or `nd_spline`, then `'spline'` (read-only)
- `nd_spline.tck`: the knots, coefficients, and degree of the spline
- `nd_spline.u`: the underlying parameter of the spline from 0 (first point) to 1 (last point)

#### Example

```python
from nd_line.nd_line import nd_line

# 3-4-5 polyline
ln = nd_line([[0, 0], [3, 0], [3, 4]])
ln.type  # 'linear'
ln.length  # 7.0
ln.lengths  # array([3., 4.])
ln.cumul  # array([0., 3., 7.])
ln.interp(1.5)  # array([1.5, 0. ])
ln.interp([0.5, 1.5])  # one point per distance
ln.interp(5.0)  # array([3., 2.])
ln.interp_rat(0.5)  # midpoint by arc length

# spline fit needs at least 4 points (cubic); original line is unchanged
ln = nd_line([[0, 0], [1, 0.5], [2, 0], [3, -0.5], [4, 0]])
spline = ln.to_spline(samples=20)
spline.type  # 'spline'
spline.points.shape  # (20, 2)
ln.type  # still 'linear'
```

#### nd_spline Example

```python
from nd_line.nd_line import nd_spline
import numpy as np

# Make a semi-circle
uu = np.linspace(0, np.pi, 8)
# input points with shape (n_points, n_dimensions)
points = np.array([np.cos(uu), np.sin(uu)]).T

nds = nd_spline(points)
halfway_point = nds.interp_rat(0.5)  # [0, 1]

# Length of a semicircle with radius=1 should be pi
nds.length  # Not close to pi because points are sparse
upsampled_spline = (
    nds.recursive_upsample()
)  # upsample to estimate spline arc length with Euclidean arc length
upsampled_spline.length  # Close to pi
```

### Notes

`interp` and `interp_rat` accept a single distance or ratio, or a sequence of them.

See the [changelog](CHANGELOG.md) for release history.
