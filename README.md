# distcorr

A python (cython) implementation of the [O(nlog(n)) algorithm](https://arxiv.org/abs/1410.1503) for computing the [distance correlation measure](https://projecteuclid.org/euclid.aos/1201012979).

# Installation

`pip install git+https://github.com/hoihui/distcorr`

# Example

```python
import distcorr
import numpy as np
v=np.linspace(-1,1,10001)
print np.corrcoef(v,np.abs(v))[0,1], distcorr.distcorr(v,np.abs(v))
```
outputs `1.53431196195e-16 0.499812417943`, demonstrating the superiority over Pearson's correlation.


# Functions

* `distcorr.distcov(x,y)`: distance covariance between array-like variables `x` and `y`
* `distcorr.distcorr(x,y)`: distance correlation between array-like variables `x` and `y`
