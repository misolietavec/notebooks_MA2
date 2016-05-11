# not IPython notebook
import matplotlib
matplotlib.rcParams['font.size'] = 18
matplotlib.rcParams['axes.titlesize'] = 18
from sympy import symbols, Heaviside as H, factor, apart, latex, Function, diff, expand_mul
from sympy import laplace_transform as LT,inverse_laplace_transform as ILT, plot, simplify
from sympy import Basic, Rational, collect, solve, numer, denom, init_printing
from numpy import poly, array, linspace

y,Y = symbols("y,Y", cls=Function)
t = symbols("t", positive=True)
s = symbols("s")
y1,y2 = symbols("y_1,y_2", cls=Function)

def to_str(n,rad):
    rs = "y'(t)" if rad==1 else "y(t)" 
    if n == 0:
        return ""
    if abs(n) !=1:
        return "%+d %s" %(n,rs)
    else:
        return "+ %s" %rs if n==1 else "- %s" %rs 
    
def odelatex(drc):
    a2,a1,a0 = drc
    s="%d y''(t)" %a2 if abs(a2)!=1 else ("y''(t)" if a2==1 else "-y''(t)")
    s += "%s %s" %(to_str(a1,1),to_str(a0,0))
    return s

def ilt_pfe(Lpf):
    Ife = []
    for pf in Lpf:
        Ife.append(ILT(pf,s,t))
    return sum(Ife)


def apart_fact(Fs):
    S=0
    Ls = apart(Fs).as_ordered_terms()
    for l in Ls:
        S += numer(l)/factor(denom(l))
    return expand_mul(S)


def plot_riesenie(t1,yt,yder,yder2):
    Py = plot((yt,(t,0,t1)),depth=12,title=r"$y(t)$",xlabel="",ylabel="")
    Pdy  = plot((yder,(t,0,t1)),title=r"$y'(t)$",xlabel="",ylabel="")
    Pddy = plot((yder2,(t,0,t1)),sings=[t1,t2],title=r"$y''(t)$",xlabel="",ylabel="")

