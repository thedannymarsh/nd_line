"""Tests for the nd_line module."""

import math
import random
import sys

import numpy as np
import pytest

sys.path.append('.')

from src.nd_line.nd_line import nd_line  # noqa E402


class Test2D:
    """Test the nd_line class for a two-dimensional line."""

    def setup_method(self):
        """Set up the 2D line."""
        random.seed(a=123)
        pts = np.array([[random.random() for _ in range(100)] for _ in range(2)])
        pts = np.transpose(pts)
        self.line = nd_line(pts)

    def test_zero(self):
        """Test interpolation at zero distance."""
        np.testing.assert_allclose(self.line.interp(0), np.array([0.0523636, 0.75081494]))

    def test_length(self):
        """Test the line length calculation."""
        np.testing.assert_allclose(self.line.length, 51.784153048659896)

    def test_interp(self):
        """Test interpolation at half distance."""
        np.testing.assert_allclose(self.line.interp(self.line.length / 2), np.array([0.11157182, 0.28764942]))

    def test_end(self):
        """Test interpolation at end distance."""
        np.testing.assert_allclose(self.line.interp(self.line.length), np.array([0.47251074, 0.41472736]))


class Test10D:
    """Test the nd_line class for a ten-dimensional line."""

    def setup_method(self):
        """Set up the 10D line."""
        random.seed(a=123)
        pts = np.array([[random.random() for _ in range(100)] for _ in range(10)])
        pts = np.transpose(pts)
        self.line = nd_line(pts)

    def test_zero(self):
        """Test interpolation at zero distance."""
        np.testing.assert_allclose(
            self.line.interp(0),
            np.array(
                [
                    0.0523636,
                    0.75081494,
                    0.5004748,
                    0.67094985,
                    0.20534254,
                    0.15618528,
                    0.15576653,
                    0.19688572,
                    0.66846312,
                    0.99964834,
                ]
            ),
        )

    def test_length(self):
        """Test the line length calculation."""
        np.testing.assert_allclose(self.line.length, 126.77686142601377)

    def test_interp(self):
        """Test interpolation at half distance."""
        np.testing.assert_allclose(
            self.line.interp(self.line.length / 2),
            np.array(
                [
                    0.66056648,
                    0.45677264,
                    0.58577795,
                    0.20230632,
                    0.0345588,
                    0.61912402,
                    0.59964352,
                    0.1829004,
                    0.26046928,
                    0.68669252,
                ]
            ),
        )

    def test_end(self):
        """Test interpolation at end distance."""
        np.testing.assert_allclose(
            self.line.interp(self.line.length),
            np.array(
                [
                    0.47251074,
                    0.41472736,
                    0.11271949,
                    0.07060848,
                    0.67520735,
                    0.00524097,
                    0.77656087,
                    0.6270458,
                    0.76898746,
                    0.92127103,
                ]
            ),
            rtol=1e-06,
        )


class TestGeometry:
    """Tests against a polyline with known segment lengths."""

    def setup_method(self):
        """3-4-5 right triangle: (0,0) -> (3,0) -> (3,4)."""
        self.points = np.array([[0.0, 0.0], [3.0, 0.0], [3.0, 4.0]])
        self.line = nd_line(self.points)

    def test_constructor_accepts_list(self):
        """Points may be passed as a nested list."""
        line = nd_line(self.points.tolist())
        np.testing.assert_array_equal(line.points, self.points)
        assert line.length == 7.0

    def test_type_is_linear(self):
        """A new line is linear, not a spline."""
        assert self.line.type == 'linear'

    def test_points_shape(self):
        """Stored points keep (n_points, n_dims)."""
        assert self.line.points.shape == (3, 2)

    def test_length_is_sum_of_segments(self):
        """Total length is the sum of Euclidean segment lengths."""
        assert self.line.length == 7.0
        segments = [nd_line.e_dist(self.points[i], self.points[i + 1]) for i in range(len(self.points) - 1)]
        np.testing.assert_allclose(self.line.length, sum(segments))

    def test_cumul(self):
        """Cumulative distance starts at 0 and ends at length."""
        np.testing.assert_array_equal(self.line.cumul, np.array([0.0, 3.0, 7.0]))
        assert self.line.cumul[0] == 0.0
        assert self.line.cumul[-1] == self.line.length

    def test_interp_at_vertices(self):
        """Interpolation at a vertex distance returns that vertex."""
        np.testing.assert_array_equal(self.line.interp(0.0), self.points[0])
        np.testing.assert_array_equal(self.line.interp(3.0), self.points[1])
        np.testing.assert_array_equal(self.line.interp(7.0), self.points[2])

    def test_interp_mid_segment(self):
        """Interpolation inside a segment is linear."""
        np.testing.assert_allclose(self.line.interp(1.5), np.array([1.5, 0.0]))
        np.testing.assert_allclose(self.line.interp(5.0), np.array([3.0, 2.0]))

    def test_interp_rat(self):
        """Ratio interpolation matches distance interpolation."""
        np.testing.assert_array_equal(self.line.interp_rat(0.0), self.points[0])
        np.testing.assert_array_equal(self.line.interp_rat(1.0), self.points[2])
        np.testing.assert_allclose(self.line.interp_rat(3.0 / 7.0), self.points[1])
        np.testing.assert_allclose(self.line.interp_rat(0.5), self.line.interp(self.line.length / 2))

    def test_interp_rejects_out_of_range(self):
        """Distance outside [0, length] raises AssertionError."""
        with pytest.raises(AssertionError):
            self.line.interp(-1e-9)
        with pytest.raises(AssertionError):
            self.line.interp(self.line.length + 1e-9)

    def test_interp_rat_rejects_above_one(self):
        """Ratio above 1 raises AssertionError."""
        with pytest.raises(AssertionError):
            self.line.interp_rat(1.0 + 1e-9)

    def test_interp_rat_negative_is_rejected(self):
        """Negative ratio fails when forwarded to interp."""
        with pytest.raises(AssertionError):
            self.line.interp_rat(-0.1)


class Test3D:
    """Test a three-dimensional polyline."""

    def setup_method(self):
        """Unit steps along x, then y, then z."""
        self.points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]])
        self.line = nd_line(self.points)

    def test_length(self):
        """Three unit segments sum to 3."""
        assert self.line.length == 3.0

    def test_interp_along_each_axis(self):
        """Midpoints of each unit segment."""
        np.testing.assert_allclose(self.line.interp(0.5), np.array([0.5, 0.0, 0.0]))
        np.testing.assert_allclose(self.line.interp(1.5), np.array([1.0, 0.5, 0.0]))
        np.testing.assert_allclose(self.line.interp(2.5), np.array([1.0, 1.0, 0.5]))


class TestEDist:
    """Tests for the Euclidean distance helper."""

    def test_zero(self):
        """Distance from a point to itself is 0."""
        pt = np.array([1.0, 2.0, 3.0])
        assert nd_line.e_dist(pt, pt) == 0.0

    def test_3_4_5(self):
        """Classic 3-4-5 triangle."""
        assert nd_line.e_dist(np.array([0.0, 0.0]), np.array([3.0, 4.0])) == 5.0

    def test_matches_math_hypot(self):
        """Agrees with math.hypot in 2D."""
        a = np.array([1.25, -0.5])
        b = np.array([4.0, 2.5])
        np.testing.assert_allclose(nd_line.e_dist(a, b), math.hypot(*(b - a)))


class TestSplineify:
    """Tests for in-place spline approximation."""

    def setup_method(self):
        """Enough points for a cubic spline (k=3)."""
        self.points = np.array(
            [
                [0.0, 0.0],
                [1.0, 0.5],
                [2.0, 0.0],
                [3.0, -0.5],
                [4.0, 0.0],
            ]
        )
        self.line = nd_line(self.points)

    def test_default_keeps_point_count(self):
        """Default samples equals the original number of points."""
        self.line.splineify()
        assert self.line.points.shape[0] == len(self.points)
        assert self.line.type == 'spline'

    def test_custom_sample_count(self):
        """Samples sets the number of points on the new polyline."""
        self.line.splineify(samples=12)
        assert self.line.points.shape == (12, 2)
        assert self.line.type == 'spline'

    def test_endpoints_preserved(self):
        """S=0 spline interpolation keeps the first and last points."""
        self.line.splineify()
        np.testing.assert_allclose(self.line.points[0], self.points[0], atol=1e-15)
        np.testing.assert_allclose(self.line.points[-1], self.points[-1], atol=1e-15)

    def test_length_is_refreshed(self):
        """Length is recomputed from the new points."""
        self.line.splineify(samples=20)
        recomputed = sum(
            nd_line.e_dist(self.line.points[i], self.line.points[i + 1]) for i in range(len(self.line.points) - 1)
        )
        np.testing.assert_allclose(self.line.length, recomputed)
        np.testing.assert_allclose(self.line.interp(0.0), self.line.points[0])
        np.testing.assert_allclose(self.line.interp(self.line.length), self.line.points[-1])
