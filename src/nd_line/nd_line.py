"""Module for creating an n-dimensional line.

Copyright Daniel Marshall
"""

import math
from functools import cached_property

import numpy as np
from numpy import ndarray
from numpy.typing import ArrayLike
from scipy.interpolate import splev, splprep


class nd_line:
    """Class for n-dimensional line."""

    def __init__(self, points: ArrayLike, name: str | None = None) -> None:
        """Create a line from a list of points.

        :param points: list of points
        :param name: line name
        """
        self._points = np.array([tuple(x) for x in points])
        self._type = 'linear'
        self.name = name  # editable

    @property
    def points(self) -> ndarray:
        """Input points from which the line was constructed."""
        return self._points  # noqa: DAR201

    @property
    def type(self) -> str:
        """Whether the line is linear or sampled from a spline."""
        return self._type  # noqa: DAR201

    @cached_property
    def lengths(self) -> ndarray:
        """Euclidean distance between each point along the line.

        Returns vector of length nd_line.points - 1.
        """  # noqa: DAR201
        return np.array([self.e_dist(self.points[i], self.points[i + 1]) for i in range(len(self.points) - 1)])

    @cached_property
    def length(self) -> float:
        """Sum of the Euclidean distance between each point along the line.

        Identical to nd_line.cumul[-1].
        """  # noqa: DAR201
        return sum(self.lengths)

    @cached_property
    def cumul(self) -> ndarray:
        """Cumulative Euclidean distance between each point along the line.

        Same length as nd_line.points.
        """  # noqa: DAR201
        return np.concatenate(([0.0], np.cumsum(self.lengths)))

    def dist_from(self, point: ArrayLike) -> ndarray:
        """Calculate the distance between a given point and the points in the nd_line.

        :param point: numpy array of a point
        :return: distance of each nd_line.points from point
        """
        # translate points relative to point of interest
        translated_points = self._points - point
        # calculate norm from new origin at point (Euclidean distance)
        return np.linalg.norm(translated_points, axis=1)

    def closest_idx(self, point: ndarray) -> int:
        """Get the index of the closest point in nd_line to a point.

        :param point: numpy array of point
        :return: index of closest point
        """
        return np.argmin(self.dist_from(point))

    def within_idx(self, point: ndarray, distance: float) -> ndarray:
        """Get the index of all points in nd_line within distance from point.

        :param point: numpy array of point
        :param distance: radius within which points will be returned
        :return: indices of nd_line points
        """
        dists = self.dist_from(point)
        return np.flatnonzero(dists[dists < distance])

    def contig_within_idx(self, point: ndarray, distance: float) -> ndarray:
        """Get the index of contiguous points in nd_line within distance from point around the closest point on nd_line.

        :param point: numpy array of point
        :param distance: radius within which points will be returned
        :return: indices of nd_line points
        """
        dists = self.dist_from(point)
        closest_idx = np.argmin(dists)
        # start at the closest point
        start_idx = closest_idx
        while dists[start_idx] < distance and start_idx > 0:
            # increment backwards in line until line point is outside boundary
            start_idx -= 1
        end_idx = closest_idx
        while dists[end_idx] < distance and end_idx < self._points.shape[0] - 1:
            # increment forward in line until line point is outside boundary
            end_idx += 1
        # return start to end
        within_idx = np.arange(start_idx, end_idx + 1)
        return within_idx

    def interp(self, dist: float | ArrayLike) -> ndarray:
        """Return a point a specified distance along the line.

        :param dist: distance along the line
        :type dist: float or 1-d array of distances
        :return: numpy array of the point coordinates
        """
        flat = np.asarray(dist, dtype=float).reshape(-1)
        if flat.size != 1:
            assert np.all(
                np.logical_and(0.0 <= flat, flat <= self.length)
            ), "All new_lengths must between 0 and nd_line length"
            return np.vstack([np.interp(flat, self.cumul, self.points[:, n]) for n in range(self.points.shape[1])]).T
        dist = float(flat[0])
        assert dist <= self.length, 'length cannot be greater than line length'
        assert dist >= 0, 'length cannot be less than zero'
        if dist == 0:
            return self.points[0]
        if dist == self.length:
            return self.points[-1]
        index = np.where(self.cumul < dist)[0][-1]
        d = self.cumul[index]
        vector = (self.points[index + 1] - self.points[index]) / self.e_dist(self.points[index], self.points[index + 1])
        remdist = dist - d
        final_point = remdist * vector + self.points[index]
        return final_point

    def interp_rat(self, ratio: float | ArrayLike) -> ndarray:
        """Return a point a specified ratio along the line.

        :param ratio: ratio along the line
        :return: numpy array of the point coordinates
        """
        flat = np.asarray(ratio, dtype=float).reshape(-1)
        if flat.size != 1:
            assert np.all(
                np.logical_and(0.0 <= flat, flat <= 1.0)
            ), "Ratio for interp_rat() must be a value from 0 to 1"
            return self.interp(flat * self.length)
        ratio = float(flat[0])
        assert ratio <= 1, "Ratio for interp_rat() must be a value from 0 to 1"
        return self.interp(ratio * self.length)

    def resample(self, new_dists: ArrayLike) -> 'nd_line':
        """Resample the line from new lengths along the line.

        :param new_dists: vector of distances along the line at which to resample
        :return: new nd_line from resampled points
        """
        new_points = self.interp(new_dists)
        return nd_line(new_points, name=self.name)

    def to_spline(self, samples: int | None = None, s: float | None = 0, **kwargs) -> 'nd_spline':  # noqa: ANN003
        """Return a new line from a spline approximation of this line.

        :param samples: number of samples to use for spline approximation
        :param s: smoothing factor for spline approximation
        :param kwargs: keyword arguments forwarded to nd_spline
        :return: new nd_spline sampled from the spline
        """
        nds = nd_spline(self.points, name=self.name, s=s, **kwargs)
        if samples is not None:
            new_dists = np.linspace(0, nds.length, samples)
            nds = nds.resample(new_dists)
        return nds

    def splinify(self, samples: int | None = None, s: float | None = 0, **kwargs) -> 'nd_spline':  # noqa: ANN003
        """Alias of nd_line.to_spline().

        :param samples: number of samples to use for spline approximation
        :param s: smoothing factor for spline approximation
        :param kwargs: keyword arguments forwarded to nd_spline
        :return: new nd_spline sampled from the spline
        """
        return self.to_spline(samples, s=s, **kwargs)

    @staticmethod
    def e_dist(a: ndarray, b: ndarray) -> float:
        """Calculate the euclidean distance between two points.

        :param a: numpy array of point a
        :param b: numpy array of point b
        :return: euclidean distance between a and b
        """
        return math.sqrt(sum((a - b) ** 2))


class nd_spline(nd_line):
    """Class for n-dimensional spline."""

    def __init__(
        self, points: ArrayLike, name: str | None = None, s: float | None = 0, **kwargs  # noqa: ANN003
    ) -> None:
        """Create a spline from a list of points.

        :param points: list of points
        :param name: line name
        :param s: smoothing parameter for splprep. Default is 0 (no smoothing). Use None for splprep default.
        :param kwargs: keyword arguments for splprep.
            See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splprep.html.
            Argument ``u`` is replaced by the normalized line lengths from 0 to 1.
        """
        super().__init__(points, name=name)
        self._type = 'spline'

        # protected attributes
        # Create u, parameter from 0 (line startpoint) to 1 (line endpoint)
        self._u = self.cumul / self.length
        # Add kwargs to the class definition for reproducability
        self.__dict__.update(kwargs)

        # u is overridden if passed in kwargs
        if 'u' in kwargs:
            del kwargs['u']

        tck, _ = splprep(self._points.T, u=self._u, s=s, **kwargs)
        self._tck = tck

    @property
    def tck(self) -> tuple:
        """Knots, coefficients, and degree from splprep.

        :return: (t, c, k) tuple.
            See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splprep.html.
        """
        return self._tck

    @property
    def u(self) -> ndarray:
        """Parameter from 0 to 1 at each point.

        :return: normalized cumulative length used by splprep.
            See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splprep.html.
        """
        return self._u

    def resample(self, new_dists: ArrayLike, **kwargs) -> 'nd_spline':  # noqa: ANN003
        """Resample the spline from new lengths along the line.

        The returned nd_spline may be dissimilar to the original nd_spline, especially if new_dists is sparse
        compared to the original nd_spline.cumul.
        The overall length of the nd_spline will be different from the original nd_spline because the length is
        approximated by Euclidean distance between points.

        :param new_dists: vector of distances along the line at which to resample
        :param kwargs: keyword arguments forwarded to the new nd_spline
        :return: new nd_spline from resampled points
        """
        new_dists = np.array(new_dists).reshape(-1)
        assert new_dists.ndim == 1, "new_lengths must be a 1-D vector of lengths"
        assert np.all(
            np.logical_and(0.0 <= new_dists, new_dists <= self.length)
        ), "All new_lengths must between 0 and nd_line length"
        new_points = self.interp(new_dists)
        return nd_spline(new_points, name=self.name, **kwargs)

    @staticmethod
    def _recursive_bisection(
        tck: tuple, u1: float, u3: float, pt1: ndarray, pt3: ndarray, d13: float, r_tol: float = 1e-4
    ) -> list[float]:
        """Recursively upsample one spline segment until it is within tolerance.

        Upsample a spline parameterized by u between u1 and u3 until the Euclidean distance between
        each point represents the spline arc length within tolerance.

        :param tck: (t, c, k) tuple from splprep.
            See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splprep.html.
        :param u1: parameter at segment start in range [0, 1]
        :param u3: parameter at segment end in range [0, 1]
        :param pt1: nd point at segment start
        :param pt3: nd point at segment end
        :param d13: distance between pt1 and pt3
        :param r_tol: relative tolerance at which to stop upsampling the spline between points.
            Upsampling stops when the proportion of dist(pt1 to pt2) + dist(pt2 to pt3) to dist(pt1 to pt3)
            is less than r_tol.
        :return: list of new parameterization values u without the first point u1
        """
        # get parameter between the boundary points
        u2 = (u1 + u3) / 2
        # evaluate the nd point between the boundaries
        pt2 = np.array(splev(u2, tck)).T
        # get distances between boundaries and new point
        d12 = nd_line.e_dist(pt1, pt2)
        d23 = nd_line.e_dist(pt2, pt3)
        # calculate proportion difference between distances pt1->pt3 and pt1->pt2->pt3
        err = (d12 + d23 - d13) / d13
        if err > r_tol:
            # if the new distances are large, bisect each new segment and return in a single list
            return nd_spline._recursive_bisection(tck, u1, u2, pt1, pt2, d12, r_tol) + nd_spline._recursive_bisection(
                tck, u2, u3, pt2, pt3, d23, r_tol
            )
        else:
            # return only the last point
            return [u3]

    def recursive_upsample(self, tol: float = 1e-4) -> 'nd_spline':
        """Add points until the polyline represents the curve.

        Add a point between existing points until the curve is reasonably represented by nd_line.points.
        This will ensure that nd_line.lengths are representative of the length along the curve and that existing points
        are preserved.

        :param tol: incremental relative tolerance at which to stop incrementing upsample.
        :return: new nd_spline from resampled points
        """
        new_u = [self._u[0]]
        for i in range(0, self._u.shape[0] - 1):
            new_u.extend(
                self._recursive_bisection(
                    self._tck, self._u[i], self._u[i + 1], self._points[i], self._points[i + 1], self.lengths[i], tol
                )
            )

        return self.resample(np.array(new_u) * self.length)

    def interp(self, dist: ArrayLike) -> ndarray:
        """Return a point a specified distance along the spline.

        :param dist: distance along the spline, a float or 1-d array of distances
        :return: numpy array of the point coordinates
        """
        new_dist = np.asarray(dist, dtype=float).reshape(-1)
        assert np.all(new_dist <= self.length), 'length cannot be greater than line length'
        assert np.all(new_dist >= 0), 'length cannot be less than zero'
        if new_dist.size == 1:
            value = float(new_dist[0])
            if value == 0.0:
                return self.points[0]
            if value == self.length:
                return self.points[-1]
        u = np.asarray(dist, dtype=float) / self.length
        return self.interp_rat(u)

    def interp_rat(self, ratio: ArrayLike) -> ndarray:
        """Return a point a specified ratio along the spline.

        :param ratio: ratio along the spline
        :return: numpy array of the point coordinates
        """
        new_ratio = np.array(ratio).reshape(-1)
        assert np.all(
            np.logical_and(0 <= new_ratio, new_ratio <= 1)
        ), "Ratio for nd_spline.interp_rat() must be a value from 0 to 1"
        return np.array(splev(new_ratio, self._tck)).T.squeeze()

    def to_line(self) -> nd_line:
        """Return a nd_line representation of the nd_spline.

        :return: nd_line of nd_spline
        """
        return nd_line(self.points, name=self.name)
