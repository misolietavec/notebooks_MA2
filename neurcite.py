from sympy import latex, symbols, Function, exp, cos, sin, solve, I 
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

def polycoeffs(Pm,m):
    Pmc = [Pm.subs(t,0)] + [Pm.coeff(t**k) for k in range(1,m+1)]
    return Pmc[::-1]

def formulacia(drc,al,Pm,be=0,Qn=[0]):
    Pm,Qn = Pm[::-1], Qn[::-1]
    LS = odelatex(drc)
    als = str(al) if al != 1 else ""
    expart = r"e^{%s t}" %als if al != 0 else ""
    Pmp,Qnp =latex(make_poly(Pm)), latex(make_poly(Qn))
    if be > 0:
        PS = r"%s\left({(%s) \cos %s t + (%s) \sin %s t}\right)" %(expart, Pmp, be, Qnp, be)
    if be == 1:
        PS = r"%s\left({(%s) \cos t + (%s) \sin t}\right)" %(expart, Pmp, Qnp)
    if be < 0:
        PS = r"%s\left({(%s) \cos %s t - (%s) \sin %s t}\right)" %(expart, Pmp, -be, Qnp, -be)
    if be == -1:
        PS = r"%s\left({(%s) \cos t - (%s) \sin t}\right)" %(expart, Pmp, Qnp)
    if be == 0:
        PS = r"%s (%s)" %(expart, Pmp)
    Zstr=r"$%s = %s.$" %(LS,PS)
    return Zstr

def nas_kor(albe, drc):
    a2,a1,a0 = drc
    kor = roots(a2*t**2 + a1*t + a0,t)
    if albe in kor.keys():
        return kor[albe]
    return 0    

def vypis_kor(drc):
    a2,a1,a0 = drc
    kor = roots(a2*t**2 + a1*t + a0,t)
    Rs = ""
    for r in kor:
        nas = kor[r]
        nstr = ', ' if nas==1 else ' (násobnosť: %d), ' %nas
        Rs += latex(r)+nstr
    return Rs[:-2]    
    
def make_poly(pc):
    n = len(pc)
    return sum([pc[k]*t**k for k in range(n)])

def make_respoly(n,k=0):
    Rc = list(symbols("r:%d" %(n+k)))
    return Rc,sum([Rc[k]*t**k for k in range(n+k)])

def zeropad(P,Q):
    np,nq = len(P),len(Q)
    if np < nq:
        P = P + [0]*(nq-np)
    elif np > nq:
        Q = Q + [0]*(np-nq)
    return [P[k] - I*Q[k] for k in range(max(np,nq))],max(np,nq) 

def sustava(LS,f,s):
    sust =[LS.subs(t,0)-f.subs(t,0)]  # absolut. clen
    for i in range(1,s):
        sust.append(expand(LS).coeff(t**i) - f.coeff(t**i))
    return sust

def riesenie(drc,al,Pm,be=0,Qn=[0],Rm=None):
    Pm,Qn = Pm[::-1], Qn[::-1]
    if Rm:
        Rm = Rm[::-1]
    a2,a1,a0 = drc
    A = al + I*be
    k = nas_kor(A,drc)
    Sc,s = zeropad(Pm,Qn)
    Ss = make_poly(Sc)
    f = Ss  # podelene exp(A*t)
    if k>0 and Rm:
        Rc,R = make_respoly(s,k)
    else:
        Rc,R = make_respoly(s)
    yp = (t**k)*exp(A*t)*R
    LS = factor((a2*diff(yp,t,2) + a1*diff(yp,t) + a0*yp)/exp(A*t))
    sust = sustava(LS,f,s)
    ries = solve(sust,Rc)
    if k>0 and Rm:
        for i in range(s):
            Rc[i+k] = ries[Rc[i]]
        for i in range(k):
            Rc[i] = Rm[i]
    else:
        for i in range(s):
            Rc[i] = ries[Rc[i]]        
    R_res = make_poly(Rc)
    if not (k>0 and Rm):
        yp_res = (t**k)*exp(A*t)*R_res
    else:    
        yp_res = exp(A*t)*R_res
    re_c,im_c = yp_res.as_real_imag()
    return re_c

def skuska(drc,yp, albe=0):
    a2,a1,a0 = drc
    LS = factor((a2*diff(yp,t,2) + a1*diff(yp,t) + a0*yp))
    if albe != 0:
        al,be = albe
        LS = collect(expand(LS/exp(al*t)),[cos(abs(be)*t),sin(abs(be)*t)])
        LS = exp(al*t)*LS          
    return LS

def vymysli(drc,al,Pm,be=0,Qn=[0],by_coeff=False):
    a2,a1,a0 = drc
    Pm,Qn = Pm[::-1], Qn[::-1]
    m,n = len(Pm), len(Qn)
    yp = exp(al*t)*(make_poly(Pm)*cos(be*t)+make_poly(Qn)*sin(be*t))
    LS = a2*diff(yp,t,2) + a1*diff(yp,t) + a0*yp
    if be != 0:
        CS = collect(expand(LS/exp(al*t)),[cos(be*t),sin(be*t)])
        if not by_coeff:
            return exp(al*t)*CS
        else:
            Pmp, Qnp = CS.coeff(cos(be*t)),CS.coeff(sin(be*t))
            s = max(m,n)
            Pm_f = [Pmp.subs(t,0)]+[Pmp.coeff(t**i) for i in range(1,s)]
            Qn_f = [Qnp.subs(t,0)]+[Qnp.coeff(t**i) for i in range(1,s)]
            return exp(al*t)*CS, Pm_f[::-1], Qn_f[::-1]
    if not by_coeff:    
        return factor(LS)
    else:
        Pmp = expand(LS/exp(al*t))
        Pm_fraw = [Pmp.coeff(t**k) for k in range(m-1,0,-1)]+[Pmp.subs(t,0)]
        Pm_f = [c for c in Pm_fraw if c != 0]
        return factor(LS), Pm_f
    