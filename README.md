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

- `ln.interp(dist)`: returns a point dist length along the arc of the line

- `ln.interp_rat(ratio)`: ratio should be a value between 0 and 1, returns a value ratio*length along the line

- `ln.splineify(samples)`: returns a new line sampled from a spline fit; `samples` is the number of points on the new line (defaults to the original count)

#### Attributes

- `ln.points`: the points of the line
- `ln.length`: the length of the line
- `ln.type`: linear if not spline approximated, spline otherwise

#### Example

```python
from nd_line.nd_line import nd_line
import numpy as np

ln = nd_line(np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]]))
interpolated_point = ln.interp(1.5)
line_length = ln.length
halfway_point = ln.interp_rat(0.5)
```

### Notes

`interp` and `interp_rat` currently accept a single distance or ratio at a time.
Batch interpolation of a list of distances is planned for a future version.

See the [changelog](CHANGELOG.md) for release history.
