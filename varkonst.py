from sympy import latex, symbols, Function, exp, cos, sin, solve, I, arg, Abs, Matrix
from sympy import diff, factor, collect, expand, roots, simplify

y = symbols("y", cls=Function)
t = symbols("t",real=True)

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

def formulacia(drc,f):
    LS = odelatex(drc)
    PS = latex(f)
    Zstr=r"$%s = %s$" %(LS,PS)
    return Zstr

def fund_ries(drc):
    a2,a1,a0 = drc
    r = roots(a2*t**2 + a1*t +a0,t)
    kl = list(r.keys())
    if len(kl) == 2:  # jednoduche korene
        r1,r2 = kl        
        re_k,im_k = r1.as_real_imag()
        if im_k == 0: # jednoduche realne korene
            return exp(r1*t), exp(r2*t)
        else:         # komplexne zdruz.korene
            return exp(re_k*t)*sin(im_k*t), exp(re_k*t)*cos(im_k*t)
    else: # dvojnasobny koren
        return exp(kl[0]*t), t*exp(kl[0]*t)

def sucin_na_sucet(f,g,u,v):
    assert(f in (sin,cos)) 
    assert(g in (sin,cos))
    u1, u2 = u-v,u+v
    if f==sin and g==sin:
        return (cos(u1) - cos(u2))/2
    elif f==cos and g==cos:
        return (cos(u1) + cos(u2))/2
    elif f==sin and g==cos:
        return (sin(u1) + sin(u2))/2
    else:  # cos, sin
        return (sin(u2) - sin(u1))/2
    