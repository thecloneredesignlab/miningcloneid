"""Finite toy-interface operations using an injected existing symbolic library."""
def positive_width(width):
    if width.is_positive is not True: raise ValueError('A declared positive width is required')

def initial_weights(sp,grid,mean,width):
    positive_width(width)
    w=sp.Matrix([sp.exp(-(z-mean)**2/(2*width**2)) for z in grid])
    return w/sum(w)

def output_jets(sp,states,dynamics,outputs,order):
    states=sp.Matrix(states);dynamics=sp.Matrix(dynamics);h=sp.Matrix(outputs)
    if dynamics.shape!=states.shape or states.cols!=1 or h.cols!=1:
        raise ValueError('State, dynamics and output must be column arrays')
    if order not in (0,1,2): raise ValueError('Jet order outside the registered fixture scope')
    result=[h]
    for k in range(order): result.append(result[-1].jacobian(states)*dynamics)
    return result

def initial_pullback(sp,jets,states,initial,unknowns):
    if len(states)!=len(initial) or len(unknowns)!=len(set(unknowns)):
        raise ValueError('Initial shape or unknown coordinate mismatch')
    substitution=dict(zip(states,initial))
    values=[j.subs(substitution,simultaneous=True) for j in jets]
    return values,[j.jacobian(unknowns) for j in values]

def kary_kernel(sp,observed,true,width):
    positive_width(width)
    raw=sp.Matrix([[sp.exp(-(o-t)**2/(2*width**2)) for t in true] for o in observed])
    result=raw.copy()
    for j in range(raw.cols):
        total=sum(raw[:,j])
        for i in range(raw.rows):result[i,j]=raw[i,j]/total
    return result

def flow_density(sp,observed,chromosomes,weights,width,n_unit,step):
    positive_width(width)
    if n_unit.is_positive is not True or step.is_positive is not True:
        raise ValueError('Positive known conversion and grid step are required')
    raw=sp.Matrix([sum(p*sp.exp(-(o-n/n_unit)**2/(2*width**2))/(sp.sqrt(2*sp.pi)*width)
                       for n,p in zip(chromosomes,weights)) for o in observed])
    return raw/(sum(raw)*step)
