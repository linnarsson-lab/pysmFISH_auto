from typing import Union
import pandas as pd
import numpy as np
from KDEpy import FFTKDE
from scipy.special import logsumexp
from scipy.interpolate import interp1d


class FasterKDE:
	def __init__(self, kernel: str = "epa", bandwidth: Union[str, float] = "ISJ"):
		self.bandwidth = bandwidth
		self.kernel = kernel
		
		self._fftkde = FFTKDE(bw=self.bandwidth, kernel=self.kernel)
	
	def fit(self, x: np.ndarray) -> None:
		self._fftkde.fit(x)
		self.x, self.y = self._fftkde.evaluate()

	def transform(self, z: np.ndarray) -> np.ndarray:
		return interp1d(self.x, self.y, kind="linear", assume_sorted=True, bounds_error=False, fill_value=10**-6)(z)
	

class BayesEEL:
	def __init__(self, kernel: str = "epa", bandwidth: Union[str, float] = "ISJ") -> None:
		self.kernel = kernel
		self.bandwidth = bandwidth
	
	def fit(self, known_x: np.ndarray, known_barcodes: np.ndarray) -> "BayesEEL":
		"""
		Fit the reference spots
		
		Args:
			known_x: the spot intensity data, of shape (n_reference_spots, n_cycles) for reference spots
			known_barcodes: a bool mask array of shape (n_reference_spots, n_cycles) indicating the correct barcodes for the reference spots
		"""
		self.known_x = known_x
		self.known_barcodes = known_barcodes
		n_spots, n_cycles = known_x.shape
		
		# Estimate the unconditional intensity probability distribution, per cycle
		self.m = []
		for i in range(n_cycles):
			x0 = known_x[:, i]
			kde = FasterKDE(bandwidth=self.bandwidth, kernel=self.kernel)
			kde.fit(x0[:, None])
			self.m.append(kde)

		# Estimate the conditional intensity probability distributions, per cycle
		self.p0 = []
		self.p1 = []
		for i in range(n_cycles):
			x0 = known_x[:, i][~known_barcodes[:, i]]
			kde = FasterKDE(bandwidth=self.bandwidth, kernel=self.kernel)
			kde.fit(x0[:, None])
			self.p0.append(kde)

			x1 = known_x[:, i][known_barcodes[:, i]]
			kde = FasterKDE(bandwidth=self.bandwidth, kernel=self.kernel)
			kde.fit(x1[:, None])
			self.p1.append(kde)
		return self
		
	def transform(self, x: np.ndarray, valid_barcodes: np.ndarray, prior_p1: float = 0.5) -> np.ndarray:
		"""
		Compute the log posterior probability for each of a set of valid candidate barcodes
		
		Args:
			x: The spot intensity data, of shape (n_spots, n_cycles)
			valid_barcodes: The candidate barcodes given as a bool mask array of shape (n_barcodes, n_cycles)
			prior_p1: the prior probability of a positive barcode (per cycle), default 0.5 (which seems to work best in practice)
		Returns:
			calls: indices of the called barcodes
			probs: posterior probabilities for the called barcodes
		"""
		n_spots, n_cycles = x.shape
		n_barcodes = valid_barcodes.shape[0]
		
		log_posterior = np.zeros(shape=(n_spots, n_barcodes))
		for b in range(n_barcodes):
			logp = np.zeros((n_spots,))
			for j in range(n_cycles):
				if valid_barcodes[b, j]:
					kde = self.p1[j]
					prior = prior_p1
				else:
					kde = self.p0[j]
					prior = 1 - prior_p1
				logp += np.log(kde.transform(x[:, j])) + np.log(prior) - np.log(self.m[j].transform(x[:, j]))
			log_posterior[:, b] = logp
		# Normalize so that the total probability across all barcodes is 1
		log_posterior = (log_posterior.T - logsumexp(log_posterior, axis=1)).T # Use logsumexp to avoid underflow
		# return log_posterior #np.minimum(log_posterior, 0)
		return np.argmax(log_posterior, axis=1), np.exp(np.max(log_posterior, axis=1))



