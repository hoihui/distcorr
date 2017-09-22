import numpy as np
cimport numpy as np
from libc.math cimport sqrt
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def DyadUpdate(np.ndarray[np.intp_t] y, np.ndarray[double] c):
    cdef np.intp_t n,rowsh,k,partialsum2,kj
    n=len(y)
    L=int(np.ceil(np.log(n)/np.log(2)))
    cdef np.ndarray[double] gamma=np.zeros(n,dtype=float)
    cdef np.ndarray[double] s=np.zeros(1<<(L+1),dtype=float)
    for i in range(1,n):
        rowsh=0
        for l in reversed(range(L)):
            k=(y[i-1]-1)>>l
            s[rowsh+k]+=c[i-1]
            rowsh+=1<<L-l
        partialsum2=0
        for bitpos in reversed(range(L)):
            if ((y[i]-1)&(1<<bitpos))!=0:
                kj=partialsum2>>bitpos
                gamma[i]+=s[(1<<L-bitpos)-2+kj]
                partialsum2+=1<<bitpos
    return gamma

@cython.boundscheck(False)
@cython.wraparound(False)
def PartialSum2D(np.ndarray[double] x,np.ndarray[double] y, np.ndarray[double] c):
    cdef np.intp_t n
    n = len(x)
    cdef np.ndarray[np.intp_t] Ix0 = np.argsort(x)
    cdef np.ndarray[np.intp_t] Ix = np.empty(len(x), int)
    Ix[Ix0] = np.arange(n)+1
    x = x[Ix0]
    y = y[Ix0]
    c = c[Ix0]
    cdef np.ndarray[np.intp_t] Iy0 = np.argsort(y)
    cdef np.ndarray[np.intp_t] Iy = np.empty(len(y), int)
    Iy[Iy0] = np.arange(n)+1

    cdef np.ndarray[double] sy = np.cumsum(c[Iy0])-c[Iy0]
    cdef np.ndarray[double] sx = np.cumsum(c) - c
    cdef double cdot = np.sum(c)

    cdef np.ndarray[double] gamma1 = DyadUpdate(Iy,c)
    cdef np.ndarray[double] gamma = cdot - c - 2*sy[Iy-1] - 2*sx + 4*gamma1
    gamma = gamma[Ix-1]
    return gamma

@cython.boundscheck(False)
@cython.wraparound(False)
def distcov(x,y):
    x=np.asarray(x,dtype=np.float64)
    y=np.asarray(y,dtype=np.float64)
    cdef np.intp_t n=len(x)
    cdef np.ndarray[np.intp_t] Ix0=np.argsort(x);
    cdef np.ndarray[np.intp_t] alphax = np.empty(len(x), int)
    alphax[Ix0] = np.arange(len(x))
    cdef np.ndarray[np.intp_t] Iy0=np.argsort(y);
    cdef np.ndarray[np.intp_t] alphay = np.empty(len(y), int)
    alphay[Iy0] = np.arange(len(y))
    cdef np.ndarray[double] sx = np.append(np.cumsum(x[Ix0]),0)
    cdef np.ndarray[double] sy = np.append(np.cumsum(y[Iy0]),0)
    cdef np.ndarray[double] betax = sx[alphax-1]
    cdef np.ndarray[double] betay = sy[alphay-1]
    cdef double xdot=np.sum(x)
    cdef double ydot=np.sum(y)

    cdef np.ndarray[double] aidot=xdot+(2*alphax-n)*x-2*betax
    cdef np.ndarray[double] bidot=ydot+(2*alphay-n)*y-2*betay
    cdef double Sab = np.sum(aidot*bidot); 

    cdef double adotdot=2*np.sum(alphax*x)-2*np.sum(betax)
    cdef double bdotdot=2*np.sum(alphay*y)-2*np.sum(betay)

    cdef np.ndarray[double] gamma_1  = PartialSum2D(x,y, np.ones(n))
    cdef np.ndarray[double] gamma_x  = PartialSum2D(x,y, x)
    cdef np.ndarray[double] gamma_y  = PartialSum2D(x,y, y)
    cdef np.ndarray[double] gamma_xy = PartialSum2D(x,y, x*y)

    cdef double aijbij = np.sum(x*y*gamma_1 + gamma_xy - x*gamma_y - y*gamma_x)
    return aijbij/n/(n-3) - 2*Sab/n/(n-2)/(n-3) + adotdot*bdotdot/n/(n-1)/(n-2)/(n-3)

@cython.boundscheck(False)
@cython.wraparound(False)
def distcorr(x,y):
    x=np.asarray(x,dtype=np.float64)
    y=np.asarray(y,dtype=np.float64)
    cdef double dcovXY = distcov(x,y)
    cdef double dcovX  = distcov(x,x)
    cdef double dcovY  = distcov(y,y)
    if np.abs(dcovX * dcovY)<1e-10:
        return 0; 
    else:
        return np.sign(dcovXY)*sqrt(abs(dcovXY)/sqrt(dcovX*dcovY))