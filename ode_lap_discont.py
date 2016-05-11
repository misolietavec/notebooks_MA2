# not IPython notebook
import matplotlib
matplotlib.use('PDF')
matplotlib.rcParams['font.size'] = 18
matplotlib.rcParams['axes.titlesize'] = 18
from sympy import symbols, Heaviside as H, factor, apart, latex, Function, diff, expand_mul
from sympy import laplace_transform as LT,inverse_laplace_transform as ILT, plot, simplify
from sympy import Basic, Rational, collect, solve, numer, denom, init_printing
from random import choice
from numpy import poly, array, linspace
from tempfile import mkstemp

y,Y = symbols("y,Y", cls=Function)
t = symbols("t", positive=True)
s = symbols("s")
y1,y2 = symbols("y_1,y_2", cls=Function)


def choice_except(a,b,exvals=None):
    R = range(a,b+1)
    if exvals:
        R = set(R) - set(exvals) 
    return choice(list(R))

def to_str(n):
    if abs(n) !=1:
        return "%+d \," %n
    else:
        return "+" if n==1 else "-" 

def odelatex(a2,a1,a0):
    s="%d \, y''(t)" %a2 if abs(a2)!=1 else ("y''(t)" if a2==1 else "-y''(t)")
    if a1:
        s = s + "%s y'(t)" %to_str(a1)
    if a0:    
        s = s + "%s y(t)" %to_str(a0)
    return s

def params(realroots=True,continuous=False):
    r1,r2=choice_except(-4,-1),choice_except(-3,0)
    if not realroots:
        if r2 == 0:
            r2=choice_except(-3,-1)
        kk1=complex(r1,r2); kk2=complex(r1,-r2)
        r1,r2=kk1,kk2
    cp=poly((r1,r2))
    brk1,brk2 = choice([1,2,3]),choice([4,5,6])
    k1 = choice_except(1,3)
    if not continuous:
        k2 = choice_except(-3,-1)
    else:
        k2 = Rational(-k1*brk1,brk2-brk1)
    y0 = choice_except(-2,2)
    yd0 = choice_except(-2,2,(0,)) if (y0==0) else choice_except(-2,2)
    return cp,(brk1,brk2),(k1,k2),(y0,yd0)


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


def solution(pars):
    (a2,a1,a0),(t1,t2),(k1,k2),(y0,yd0)=pars
    a2,a1,a0 = map(int,[a2,a1,a0])
    den = a2*s**2+a1*s+a0
    V1 = apart(1/s/den).as_ordered_terms()            # const.
    V2 = apart(1/s**2/den).as_ordered_terms()         # lin.
    Iv1,Iv2 = ilt_pfe(V1),ilt_pfe(V2)
    LIC  = (a2*(s*y0+yd0) + a1*y0)/den
    y1 = ILT(LIC, s,t)
    y2 = k1*Iv2                                       #        k1,  t
    y3 = (k2-k1)*Iv2.subs(t,t-t1)*H(t-t1)             #  k2-k1,(t-t1)*H(t-t1)
    y4 = (t1*(k2-k1)-k2*t2)*Iv1.subs(t,t-t1)*H(t-t1)  # t1*(k2-k1)-k2*t2,H(t-t1), nulove pre spojite
    y5 = -k2*Iv2.subs(t,t-t2)*H(t-t2)                 #  -k2, (t-t2)*H(t-t2)
    y = expand_mul(simplify(y1+y2+y3+y4+y5))
    return y


def formulation(pars):
    (a2,a1,a0),(t1,t2),(k1,k2),(y0,yd0)=pars
    # a2,a1,a0 = map(int,[a2,a1,a0])
    # LS = latex(a2*diff(y(t),t,2) + a1*diff(y(t),t) + a0*y(t),order='rev-lex')
    LS = odelatex(a2,a1,a0)
    f1, f2 = latex(k1*t), latex(-k2*(t2-t))
    Zstr=r"""Riešme rovnicu $%s = f(t)$ s počiatočnými podmienkami
    $$y(0)=%d,\, y'(0)=%d.$$ 
    Pravá strana $f(t)$ je daná rovnicami:
    $$f(t)=\left\{{
    \begin{array}{ll}
    %s & \mbox{ pre } 0\le t \le %d,\\
    \displaystyle{%s} & \mbox{ pre } %d < t \le %d,\\
    %s & \mbox{ pre } t \ge %d.
    \end{array}
    }\right. $$
    """ %(LS,y0,yd0,f1,t1,f2,t1,t2,latex(0),t2)
    return Zstr


def plotexprt(ft,sings=[],depth=10,xlabel='t',ylabel='y(t)'):
    P=plot(ft,singularities=sings,depth=depth,xlabel=xlabel,ylabel=ylabel)

def plotnp(ft,t1,t2, xlab=r'$t$',ylab=r'$y(t)$',npts=200):
    tn = linspace(t1,t2,npts)
    yn = array([ft(ta) for ta in tn],dtype='float')
    plot(tn,yn)
    xlabel(xlab); ylabel(ylab)

def explanation(pars,for_tex=False):
    (a2,a1,a0),(t1,t2),(k1,k2),(y0,yd0)=pars
    a2,a1,a0 = map(int,[a2,a1,a0])
    Chpol = a2*s**2+a1*s+a0
    r12 = solve(Chpol,s)
    if len(r12)==2:
        r1,r2=r12
    else:
        r1=r12[0]; r2=r1
    # nezavisle od parametrov    
    realroots = (complex(r1).imag == 0)
    kps = t1*(k2-k1)-k2*t2
    continuous = (kps == 0)
    
    if not for_tex:
        beg_eqn = r"$$\begin{eqnarray*}"
        end_eqn = r"\end{eqnarray*}$$"
    else:
        beg_eqn = r"\begin{eqnarray*}"
        end_eqn = r"\end{eqnarray*}"        
    Lls = Chpol * Y(s)
    IC = a1*y0+s*y0+yd0
    ict = IC/Chpol
    y0t = simplify(ILT(ict,s,t))
    if realroots:
        y0t = H(t)*expand_mul(simplify(y0t/H(t)))
    if not continuous:
        V1 = apart(1/s/Chpol).as_ordered_terms()
        Iv1 = ilt_pfe(V1)
        y1t = expand_mul(simplify(Iv1/H(t)))
    V2 = apart(1/s**2/Chpol).as_ordered_terms(); Iv2 = ilt_pfe(V2)
    y2t = expand_mul(simplify(Iv2/H(t)))
    if realroots:
        icpfe=latex(apart_fact(ict))
    else:
        icpfe=latex(ict)
    if not continuous:
        Ykstr=r"""%s
        Y_0(s) &=& %s,\\
        Y_1(s) &=& \frac{1}{s (%s)},\\ 
        Y_2(s) &=& \frac{1}{s^2 (%s)}.\\
        %s""" %(beg_eqn,latex(ict),latex(Chpol),latex(Chpol),end_eqn)
        ytstr=r"""%s
        y_0(t) &=& %s,\\ 
        y_1(t) &=& \theta(t) \left({%s}\right),\\
        y_2(t) &=& \theta(t) \left({%s}\right).\\
        %s""" %(beg_eqn,latex(y0t,fold_func_brackets=True), 
                               latex(y1t,fold_func_brackets=True),
                               latex(y2t,fold_func_brackets=True),end_eqn)
        pfestr=r"""%s
        Y_0(s) &=& %s,\\
        Y_1(s) &=& %s,\\ 
        Y_2(s) &=& %s.\\
        %s""" %(beg_eqn,icpfe,latex(apart(1/s/Chpol)),latex(apart(1/s**2/Chpol)),end_eqn) 
    else:
        Ykstr=r"$$Y_0(s) = %s,\ Y_1(s) = \frac{1}{s^2 (%s)}.$$" %(latex(ict),latex(Chpol))
        ytstr=r"""%s
        y_0(t) &=& %s,\\
        y_1(t) &=& \theta(t) \left({%s}\right).\\
        %s""" %(beg_eqn,latex(y0t,fold_func_brackets=True),
                               latex(y2t,fold_func_brackets=True),end_eqn)
        pfestr=r"$$Y_0(s) =%s,\ Y_1(s) = %s.$$" %(icpfe,latex(apart_fact(1/s**2/Chpol)))
    ic1,ic2 = s*Y(s)- y0, s**2*Y(s) - y0*s - yd0
    f = k1*t*(H(t)-H(t-t1))+k2*(t-t2)*(H(t-t1)-H(t-t2))
    ct1,ct2=factor((k2-k1)*(t-t1)),factor(-k2*(t-t2))
    fh =k1*t*H(t)+ct1*H(t-t1)+ct2*H(t-t2)
    if not continuous:
        fh = fh + kps*H(t-t1)
    Lfh = expand_mul(LT(fh,t,s)[0])
    y2t1 = y2t.subs(t,t-t1); y2t2 = y2t.subs(t,t-t2)
    if not continuous:
        y1t1 = y1t.subs(t,t-t1)
    yt =  y0t + k1*y2t + (k2-k1)*y2t1*H(t-t1) - k2*y2t2*H(t-t2)
    if not continuous:
        yt = yt + kps*y1t1*H(t-t1)
    yder = diff(y0t/H(t),t)*H(t) + k1*diff(y2t,t)*H(t)+(k2-k1)*diff(y2t1,t)*H(t-t1)-k2*diff(y2t2,t)*H(t-t2)
    if not continuous:
        yder = yder + kps*diff(y1t1,t)*H(t-t1)
        
    yder2 = diff(y0t/H(t),t,2)*H(t) + k1*diff(y2t,t,2)*H(t)+(k2-k1)*diff(y2t1,t,2)*H(t-t1)-k2*diff(y2t2,t,2)*H(t-t2)
    if not continuous:
        yder2 = yder2 + kps*diff(y1t1,t,2)*H(t-t1)
    if continuous:    
        sol2 = latex(k1 * y1(t) + (k2-k1)*y1(t-t1) - k2 * y1(t-t2))
    else:
        sol2 = latex(k1 * y2(t) + (k2-k1)*y2(t-t1) - k2 * y2(t-t2) + kps*y1(t-t1))
    sol2 = sol2.replace('\\operatorname','')
    if kps < 0:
        SolStr = "y(t) = y_0(t) " + sol2
    else:    
        SolStr = "y(t) = y_0(t) + " + sol2
        
    Fcheck=a2*yder2+a1*yder+a0*yt
    Is_solution=simplify(Fcheck-f)
    Estr=[r"""Nech $y(t) \rightarrow Y(s)$. Pre transformáciu ľavej strany použijeme vzorce:
    %s
    y'(t)  &\rightarrow & s Y(s)- y(0)=%s,\\ 
    y''(t) &\rightarrow & s^2 Y(s) - s y(0) -y'(0)= %s.\\
    %s""" %(beg_eqn,latex(ic1), latex(ic2),end_eqn),
    r"""Pravú stranu prepíšeme pomocou Heavisideovho skoku $\theta(t)\,$ takto:
    $$f(t)=%s.$$
    Keď chceme mať rovnaké posunutia vo funkcii aj v Heavisidovom skoku, upravíme to
    $$f(t)=%s.$$
    Je to treba kvôli použitiu pravidla o posunutí originálu
    $$\mbox{ ak } f(t) \rightarrow F(s), 
    \mbox{ potom } f(t-a) \rightarrow \mathrm{e}^{-a s} F(s),$$
    pre transformáciu pravej strany.
    Teraz môžeme napísať L-transformovanú rovnicu (nezlomkové členy vznikli z počiatočných podmienok, ktoré sme zahrnuli do pravej strany):
    $$%s = %s.$$
    """ %(latex(f),latex(fh),latex(Lls),latex(Lfh +a1*y0 +s* y0 + yd0)),
    r"""Korene charakteristickej rovnice sú $r_1=%s,\, r_2=%s.$
    Vzhľadom na horeuvedené pravidlo o posunutí originálu a linearitu L-transformácie 
    stačí hľadať originály len pre výrazy %s
    Urobili by sme to rozkladom na elementárne zlomky, ale robotu nám ušetrí softvér
    ($\textit{sympy}$) a dostaneme %s
    Pre originály máme %s
    Nakoniec, použitím pravidla o posunutí a linearity L-transformácie dostaneme celkové riešenie
    $$%s.$$ 
    Pomocou $\textit{sympy}$ overíme, či je to riešenie. Je, lebo rozdiel ľavej a pravej strany je $%s$.
    """ %(latex(r1),latex(r2),Ykstr, pfestr, ytstr, SolStr, latex(Is_solution))
    ]
    return Estr[0]+Estr[1]+Estr[2], t1, t2, f, yt, yder, yder2 

def save_plots(t1,t2,f,yt,yder,yder2):
    Pf = plot((f,(t,0,t2+1)),title=r"$f(t)$",xlabel="",ylabel="",show=False)  # vypis formulacie
    Py = plot((yt,(t,0,t2+1)),depth=12,title=r"$y(t)$",xlabel="",ylabel="",show=False)
    Pdy  = plot((yder,(t,0,t2+1)),title=r"$y'(t)$",xlabel="",ylabel="",show=False)
    Pddy = plot((yder2,(t,0,t2+1)),sings=[t1,t2],title=r"$y''(t)$",xlabel="",ylabel="",show=False)
    ffn = mkstemp(suffix=".pdf",dir=".")[1]
    yfn = mkstemp(suffix=".pdf",dir=".")[1]
    dyfn  = mkstemp(suffix=".pdf",dir=".")[1]
    ddyfn = mkstemp(suffix=".pdf",dir=".")[1]
    Pf.save(ffn)
    Py.save(yfn)
    Pdy.save(dyfn)
    Pddy.save(ddyfn)
    return ffn, yfn, dyfn, ddyfn

def plot_riesenie(t1,t2,yt,yder,yder2):
    Py = plot((yt,(t,0,t2+1)),depth=12,title=r"$y(t)$",xlabel="",ylabel="")
    Pdy  = plot((yder,(t,0,t2+1)),title=r"$y'(t)$",xlabel="",ylabel="")
    Pddy = plot((yder2,(t,0,t2+1)),sings=[t1,t2],title=r"$y''(t)$",xlabel="",ylabel="")

def to_latex(pars,fname="odediscont.tex"):
    Z = formulation(pars)
    Estr, t1, t2, f, yt, yder, yder2 = explanation(pars,for_tex=True)
    ffn, yfn, dyfn, ddyfn = save_plots(t1, t2, f, yt, yder, yder2)
    header = r"""\documentclass[a4paper,10pt]{article}
    \usepackage[utf8]{inputenc}
    \usepackage[T1]{fontenc}
    \usepackage[slovak]{babel}
    \usepackage{a4wide}
    \usepackage{graphicx}
    \begin{document}
    
    \noindent
    """
    LtxStr=header + Z + Estr
    of = open(fname,'w')
    of.write(LtxStr)
    of.flush()
    PlotStr=r"""\begin{figure}
    \includegraphics[width=0.45\textwidth]{%s}
    \includegraphics[width=0.45\textwidth]{%s}
    \caption{Pravá strana a riešenie.}
    \end{figure}
    \begin{figure}
    \includegraphics[width=0.45\textwidth]{%s}
    \includegraphics[width=0.45\textwidth]{%s}
    \caption{Prvá a druhá derivácia riešenia.}
    \end{figure}    
    """ % (ffn, yfn, dyfn, ddyfn)
    of.write(PlotStr + r"\end{document}") 
    of.close()
    