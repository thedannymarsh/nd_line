# nd_line

[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?logo=github)](https://github.com/thedannymarsh/nd_line)
[![PyPI](https://img.shields.io/pypi/v/nd_line.svg)](https://pypi.org/project/nd_line/)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/nd_line?period=month&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads/month)](https://pepy.tech/projects/nd_line)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/nd_line?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/nd_line)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/nd_line.svg)](https://pypi.org/project/nd_line/)

Interpolate points on an n-dimensional line by euclidean arc length.

### Installation

`pip install nd_line`

#### Methods

- `ln.interp(dist)`: returns a point `dist` along the arc of the line
- `ln.interp_rat(ratio)`: `ratio` is between 0 and 1; returns the point at `ratio * length`
- `ln.to_spline(samples)`: returns a new line sampled from a spline fit; `samples` is the number of points on the new line (defaults to the original count)

#### Attributes

- `ln.points`: the points of the line
- `ln.length`: the length of the line
- `ln.lengths`: Euclidean length of each segment
- `ln.cumul`: cumulative distance at each point
- `ln.type`: `'linear'` unless created by `to_spline`, then `'spline'`

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
ln.interp(5.0)  # array([3., 2.])
ln.interp_rat(0.5)  # midpoint by arc length

# spline fit needs at least 4 points (cubic); original line is unchanged
ln = nd_line([[0, 0], [1, 0.5], [2, 0], [3, -0.5], [4, 0]])
spline = ln.to_spline(samples=20)
spline.type  # 'spline'
spline.points.shape  # (20, 2)
ln.type  # still 'linear'
```

### Notes

Currently points must be sampled one at a time, future version will allow interpolation of a list of distances along the line
