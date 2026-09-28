"""New covariant Dirac trace and two-state algebra; SymPy verifies identities.

No previous project code is imported. C++ evaluates physical polarization
vectors independently. The contractions below use only the Clifford algebra.
"""
from itertools import product
from pathlib import Path
import sympy as S

def matchings(indices):
    if not indices:
        yield 1, []
        return
    first = indices[0]
    for position in range(1, len(indices)):
        rest = indices[1:position] + indices[position+1:]
        for sign, pairs in matchings(rest):
            yield (-1)**(position-1)*sign, [(first, indices[position])] + pairs

def trace_word(word, gram):
    """Contract two repeated Lorentz labels using connected components."""
    total = 0
    for sign, edges in matchings(list(range(len(word)))):
        parent = list(range(len(word)))
        def root(i):
            while parent[i] != i:
                i = parent[i]
            return i
        contractions = edges.copy()
        for label in ('a', 'b'):
            positions = [i for i, x in enumerate(word) if x == label]
            if not positions:
                continue
            assert len(positions) == 2
            contractions.append(tuple(positions))
        for i,j in contractions:
            parent[root(i)] = root(j)
        groups = {}
        for i, label in enumerate(word):
            groups.setdefault(root(i), []).append(label)
        term = 4*sign
        for labels in groups.values():
            ends = [x for x in labels if x not in ('a','b')]
            assert len(ends) in (0,2)
            term *= 4 if not ends else gram[ends[0],ends[1]]
        total += term
    return S.expand(total)

def derive_trace():
    a,b,c,d = S.symbols('a b c d', nonzero=True)
    # p,k,q with p^2=c, k^2=0, q^2=d, a=2p.k, b=d-2p.q.
    basis = S.Matrix([[c,a/2,(d-b)/2],[a/2,0,(a+b)/2],[(d-b)/2,(a+b)/2,d]])
    components = {'p':S.Matrix([1,0,0]), 'f':S.Matrix([1,1,-1]), 'r':S.Matrix([1,1,0]), 'h':S.Matrix([1,0,-1])}
    gram = {(x,y):(vx.T*basis*vy)[0] for x,vx in components.items() for y,vy in components.items()}
    vertices = [('b','r','a'), ('a','h','b')]
    answer = 0
    for i,j in product(range(2),repeat=2):
        full = ['f',*vertices[i],'p',*vertices[j][::-1]]
        slots = [k for k,x in enumerate(full) if x in components]
        subtotal = 0
        for retain in product([False,True],repeat=4):
            omit = {k for k,keep in zip(slots,retain) if not keep}
            word = [x for k,x in enumerate(full) if k not in omit]
            if len(word)%2 == 0:
                subtotal += c**(len(omit)//2)*trace_word(word,gram)
        answer += subtotal/([a,b][i]*[a,b][j]*4)
    compact = -2*(a/b+b/a)+4*(2*c+d)*(1/a+1/b+c*(1/a+1/b)**2-d/(a*b))
    assert S.factor(answer-compact) == 0
    E,w,m = S.symbols('E w m', positive=True)
    cost = 1-m*(1/w-1/E)
    assert S.simplify(compact.subs({a:2*m*E,b:-2*m*w,c:m*m,d:0})-2*(E/w+w/E-1+cost**2)) == 0
    assert S.simplify(compact-compact.xreplace({a:b,b:a})) == 0
    print('PASS: fresh covariant trace equals compact F; Klein-Nishina differential identity; s/u symmetry')
    return str(compact)

def derive_mixing():
    eps,M2,g2,E,G,L = S.symbols('epsilon M2 mg2 E Gamma L', positive=True)
    matrix = S.Matrix([[g2+eps**2*M2,-eps*M2],[-eps*M2,M2]])
    discriminant = S.expand(S.trace(matrix)**2-4*matrix.det())
    expected = (M2-g2)**2+2*eps**2*M2*(M2+g2)+eps**4*M2**2
    assert S.expand(discriminant-expected)==0
    delta,beta = S.symbols('delta beta', real=True)
    q = S.I*delta-G/2
    integral = beta*(S.exp(q*L)-1)/q
    assert S.simplify(S.diff(integral,L)-beta*S.exp(q*L))==0
    # The thick-source amplitude modulus squared: beta=epsilon*M2/(2E).
    assert S.simplify((beta**2/(delta**2+G**2/4)).subs({beta:eps*M2/(2*E),delta:(g2-M2)/(2*E)})-eps**2*M2**2/((g2-M2)**2+E**2*G**2))==0
    print('PASS: Hamiltonian discriminant including epsilon^4; damped-amplitude primitive; thick-source probability')

def derive_absorption_polarization():
    a,b,m,d = S.symbols('a b m d', nonzero=True)
    E=(a-d)/(2*m); w=-b/(2*m); p2=E**2-d
    # Basis p,q,k,R with R=M*|q|*epsilon_L. R^2=-M^2 |q|^2.
    basis=S.Matrix([[m*m,(a-d)/2,-b/2,m*p2],[(a-d)/2,d,(a+b)/2,0],[-b/2,(a+b)/2,0,-d*w+E*(a+b)/2],[m*p2,0,-d*w+E*(a+b)/2,-d*p2]])
    components={'p':S.Matrix([1,0,0,0]),'f':S.Matrix([1,1,-1,0]),'r':S.Matrix([1,1,0,0]),'h':S.Matrix([1,0,-1,0]),'R':S.Matrix([0,0,0,1])}
    gram={(x,y):(vx.T*basis*vy)[0] for x,vx in components.items() for y,vy in components.items()}
    vertices=[('a','r','R'),('R','h','a')]
    answer=0
    for i,j in product(range(2),repeat=2):
        full=['f',*vertices[i],'p',*vertices[j][::-1]]
        slots=[k for k,x in enumerate(full) if x in ('f','p','r','h')]
        subtotal=0
        for retain in product([False,True],repeat=4):
            omit={k for k,keep in zip(slots,retain) if not keep}
            word=[x for k,x in enumerate(full) if k not in omit]
            if len(word)%2==0:
                subtotal+=m**len(omit)*trace_word(word,gram)
        answer-=subtotal/(2*[a,b][i]*[a,b][j]*d*p2)
    answer=S.factor(answer)
    assert S.simplify(S.limit(answer,d,0))==0
    print('PASS: absorption longitudinal trace derived; vanishes in the massless limit')
    print('F_L =',answer)
    Path(__file__).with_name('output').joinpath('longitudinal.tex').write_text(S.latex(answer)+'\n')

def main():
    expression=derive_trace()
    derive_mixing()
    derive_absorption_polarization()
    Path(__file__).with_name('output').joinpath('derived_kernel.txt').write_text(expression+'\n')

if __name__=='__main__':
    main()
