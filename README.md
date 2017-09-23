# distcorr

A python (cython) implementation of the [O(nlog(n)) algorithm](https://arxiv.org/abs/1410.1503) for computing the [distance correlation measure](https://projecteuclid.org/euclid.aos/1201012979).

# Installation

`pip install git+https://github.com/hoihui/distcorr`

# Example

```python
import distcorr
import numpy as np
N=1000
print distcorr.distcorr(np.random.random(N),np.random.random(N))
```

# Functions

* `distcorr.distcov(x,y)`: distance covariance between array-like variables `x` and `y`
* `distcorr.distcorr(x,y)`: distance correlation between array-like variables `x` and `y`
