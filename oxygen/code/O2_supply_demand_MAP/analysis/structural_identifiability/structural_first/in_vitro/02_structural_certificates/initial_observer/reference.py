"""Independent analytical references, without importing the adapter."""
def finite_family(sp,grid,mean,width):
    # Exponential-family form cancels the common mean-square factor separately.
    raw=[sp.exp(z*mean/width**2-z**2/(2*width**2)) for z in grid]
    p=sp.Matrix([v/sum(raw) for v in raw])
    ez=sum(z*v for z,v in zip(grid,p))
    esq=sum((z-mean)**2*v for z,v in zip(grid,p))
    dm=sp.Matrix([v*(z-ez)/width**2 for z,v in zip(grid,p)])
    ds=sp.Matrix([v*((z-mean)**2-esq)/width**3 for z,v in zip(grid,p)])
    return p,dm,ds

def manual_kernel(sp,width):
    u=sp.exp(-1/(2*width**2))
    raw=sp.Matrix([[1,u,u**4],[u,1,u],[u**4,u,1]])
    derivatives=sp.Matrix([[0,u/width**3,4*u**4/width**3],
                           [u/width**3,0,u/width**3],
                           [4*u**4/width**3,u/width**3,0]])
    values=sp.zeros(3);dw=sp.zeros(3)
    for j in range(3):
        c=sum(raw[:,j]);dc=sum(derivatives[:,j])
        for i in range(3):
            values[i,j]=raw[i,j]/c
            dw[i,j]=(derivatives[i,j]*c-raw[i,j]*dc)/c**2
    return values,dw
