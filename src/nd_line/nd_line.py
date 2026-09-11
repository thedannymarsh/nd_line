"""Module for creating an n-dimensional line.

Copyright Daniel Marshall
"""

import math
from typing import Optional

import numpy as np
from numpy import ndarray
from numpy.typing import ArrayLike
from scipy.interpolate import splev, splprep


class nd_line:
    """Class for n-dimensional line."""

    def __init__(self, points: ArrayLike) -> None:
        """Create a line from a list of points.

        :param points: list of points
        """
        self._points = np.array([tuple(x) for x in points])
        self.type = 'linear'

    @property
    def points(self) -> ndarray:
        """Input points from which the line was constructed."""
        return self._points  # noqa: DAR201

    @property
    def lengths(self) -> ndarray:
        """Euclidean distance between each point along the line.

        Returns vector of length nd_line.points - 1.
        """  # noqa: DAR201
        return np.array([self.e_dist(self.points[i], self.points[i + 1]) for i in range(len(self.points) - 1)])

    @property
    def length(self) -> float:
        """Sum of the Euclidean distance between each point along the line.

        Identical to nd_line.cumul[-1].
        """  # noqa: DAR201
        return sum(self.lengths)

    @property
    def cumul(self) -> ndarray:
        """Cumulative Euclidean distance between each point along the line.

        Same length as nd_line.points.
        """  # noqa: DAR201
        return np.concatenate(([0.0], np.cumsum(self.lengths)))

    def interp(self, dist: float) -> ndarray:
        """Return a point a specified distance along the line.

        :param dist: distance along the line
        :type dist: float
        :return: numpy array of the point coordinates
        """
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

    def interp_rat(self, ratio: float) -> ndarray:
        """Return a point a specified ratio along the line.

        :param ratio: ratio along the line
        :return: numpy array of the point coordinates
        """
        assert ratio <= 1, "Ratio for interp_rat() must be a value from 0 to 1"
        return self.interp(ratio * self.length)

    def splineify(self, samples: Optional[int] = None, s: float = 0) -> None:
        """Turn line into a spline approximation, currently occurs in place.

        :param samples: number of samples to use for spline approximation
        :param s: smoothing factor for spline approximation
        """
        if samples is None:
            samples = len(self.points)
        tck, u, _, _, _ = splprep([self.points[:, i] for i in range(self.points.shape[1])], s=s)
        self._points = np.transpose(splev(np.linspace(0, 1, num=samples), tck))
        self.type = 'spline'

    @staticmethod
    def e_dist(a: ndarray, b: ndarray) -> float:
        """Calculate the euclidean distance between two points.

        :param a: numpy array of point a
        :param b: numpy array of point b
        :return: euclidean distance between a and b
        """
        return math.sqrt(sum((a - b) ** 2))
