"""Exact certificates for the sharpness theorem of proof.tex.

For a grid with (N-3)(j-3) <= 4, the configurations with prescribed rectangle
constants form a family along which every edge ratio varies.  The certificate
is the rank of a Jacobian at a point over GF(p^2).  Because the defining
equations are monic in each ratio, the coordinate ring is a free module over
the polynomial ring in the positions, so a minor that is nonzero at any point
over any field is nonzero at the generic point.
"""
import random
import unittest

PRIME=1000003
assert PRIME%4==3


class GF2:
    """Elements a+bi of GF(p^2) with i^2=-1 (p = 3 mod 4)."""
    __slots__=('a','b')
    def __init__(s,a,b=0):s.a=a%PRIME;s.b=b%PRIME
    def __add__(s,o):return GF2(s.a+o.a,s.b+o.b)
    def __sub__(s,o):return GF2(s.a-o.a,s.b-o.b)
    def __neg__(s):return GF2(-s.a,-s.b)
    def __mul__(s,o):return GF2(s.a*o.a-s.b*o.b,s.a*o.b+s.b*o.a)
    def inv(s):
        n=(s.a*s.a+s.b*s.b)%PRIME;ni=pow(n,PRIME-2,PRIME)
        return GF2(s.a*ni,-s.b*ni)
    def __truediv__(s,o):return s*o.inv()
    def iszero(s):return s.a==0 and s.b==0


def sqrt_of_prime_field_element(d):
    d%=PRIME
    if d==0:return GF2(0)
    if pow(d,(PRIME-1)//2,PRIME)==1:return GF2(pow(d,(PRIME+1)//4,PRIME))
    return GF2(0,pow((-d)%PRIME,(PRIME+1)//4,PRIME))


def row_reduce(rows):
    rows=[list(r) for r in rows];m=len(rows);n=len(rows[0]);rank=0;pivots=[]
    for col in range(n):
        piv=next((i for i in range(rank,m) if not rows[i][col].iszero()),None)
        if piv is None:continue
        rows[rank],rows[piv]=rows[piv],rows[rank]
        inv=rows[rank][col].inv();rows[rank]=[x*inv for x in rows[rank]]
        for i in range(m):
            if i!=rank and not rows[i][col].iszero():
                f=rows[i][col];rows[i]=[x-f*y for x,y in zip(rows[i],rows[rank])]
        pivots.append(col);rank+=1
        if rank==m:break
    return rank,rows[:rank],pivots


def kernel(rows):
    rank,reduced,pivots=row_reduce(rows);n=len(rows[0])
    basis=[]
    for free in [c for c in range(n) if c not in pivots]:
        v=[GF2(0)]*n;v[free]=GF2(1)
        for r,pc in enumerate(pivots):v[pc]=-reduced[r][free]
        basis.append(v)
    return rank,basis


def configuration_jacobian(N,J,seed):
    """Point of the configuration variety over GF(p^2) and the Jacobian of the
    distance relations and rectangle constants there."""
    rng=random.Random(seed)
    P=[(rng.randrange(PRIME),rng.randrange(PRIME)) for _ in range(N)]
    R=[(rng.randrange(PRIME),rng.randrange(PRIME)) for _ in range(J)]
    two,one=GF2(2),GF2(1);lam={}
    for i in range(N):
        for l in range(J):
            L=((P[i][0]-R[l][0])**2+(P[i][1]-R[l][1])**2)%PRIME
            t=(2-L)%PRIME
            lam[(i,l)]=(GF2(t)-sqrt_of_prime_field_element((t*t-4)%PRIME))/two
            assert (lam[(i,l)]+lam[(i,l)].inv()-GF2(t)).iszero()
    nvar=2*(N+J)+N*J
    def ip(i,c):return 2*i+c
    def ir(l,c):return 2*N+2*l+c
    def il(i,l):return 2*(N+J)+i*J+l
    rows=[]
    for i in range(N):
        for l in range(J):
            row=[GF2(0)]*nvar;dx=GF2(P[i][0]-R[l][0]);dy=GF2(P[i][1]-R[l][1])
            row[ip(i,0)]=two*dx;row[ip(i,1)]=two*dy;row[ir(l,0)]=-(two*dx);row[ir(l,1)]=-(two*dy)
            row[il(i,l)]=one-(lam[(i,l)]*lam[(i,l)]).inv()
            rows.append(row)
    for i in range(1,N):
        for l in range(1,J):
            chi=lam[(i,l)]*lam[(0,0)]/(lam[(i,0)]*lam[(0,l)])
            row=[GF2(0)]*nvar
            row[il(i,l)]=lam[(0,0)];row[il(0,0)]=lam[(i,l)]
            row[il(i,0)]=-(chi*lam[(0,l)]);row[il(0,l)]=-(chi*lam[(i,0)])
            rows.append(row)
    return rows,nvar,il


class SharpnessCertificates(unittest.TestCase):
    MAXIMAL_SHAPES=[(7,4),(4,7),(5,5),(9,3),(3,9),(14,2),(2,14)]

    def test_every_grid_below_thirty_is_in_a_maximal_shape(self):
        for n in range(1,30):
            for j in range(1,30):
                if n*j<30 and n>1 and j>1:
                    self.assertTrue(any(n<=a and j<=b for a,b in self.MAXIMAL_SHAPES),(n,j))

    def test_dimension_count(self):
        for n in range(1,40):
            for j in range(1,40):
                internal=2*(n+j)-(n-1)*(j-1)-3
                self.assertEqual(internal>=1,(n-3)*(j-3)<=4)
        self.assertEqual(min(n*j for n in range(4,40) for j in range(4,40) if (n-3)*(j-3)>=5),30)

    def test_ratios_vary_along_the_fibres(self):
        for (N,J) in self.MAXIMAL_SHAPES:
            rows,nvar,il=configuration_jacobian(N,J,seed=100*N+J)
            rank,basis=kernel(rows)
            # Full rank: the fibre is smooth of the expected dimension at this point.
            self.assertEqual(rank,len(rows),(N,J))
            self.assertEqual(len(basis),nvar-rank)
            self.assertGreaterEqual(len(basis)-3,1)
            # Every ratio has a tangent direction along which it changes.
            for i in range(N):
                for l in range(J):
                    self.assertFalse(all(v[il(i,l)].iszero() for v in basis),(N,J,i,l))

    def test_thirty_edge_grids_are_rigid_modulo_isometry(self):
        for (N,J) in [(6,5),(5,6),(8,4)]:
            rows,nvar,il=configuration_jacobian(N,J,seed=7*N+J)
            rank,basis=kernel(rows)
            self.assertEqual(len(basis),3,(N,J))


if __name__=='__main__':unittest.main()
