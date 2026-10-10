"""Exact-rational checks for the joint 6-by-5 system of Section 4 of proof.tex.

These verify the rank statements on concrete configurations (generic, and on
a circle) and the kernel's constant scalar part. They illustrate the cited
rigidity theorem in the cases used; they are not a proof of it.
"""
from fractions import Fraction as Q
import random
import unittest
import sympy as s

N,J=6,5
EDGES=[(i,l) for i in range(N) for l in range(J)]
# Spanning tree of proof.tex: p1r4 p2r4 p3r5 p4r5 p5r1 p5r2 p6r3 p6r4 p6r5 p5r5 (1-based).
TREE={(0,3),(1,3),(2,4),(3,4),(4,0),(4,1),(5,2),(5,3),(5,4),(4,4)}
CORE=[(i,l) for (i,l) in EDGES if i<4 and l<3]


def rigidity(pts,rpts,edges):
    rows=[]
    for (i,l) in edges:
        d=(pts[i][0]-rpts[l][0],pts[i][1]-rpts[l][1])
        row=[0]*(2*(N+J))
        row[2*i],row[2*i+1]=d[0],d[1]
        row[2*N+2*l],row[2*N+2*l+1]=-d[0],-d[1]
        rows.append(row)
    return s.Matrix(rows)


def joint_system(pts,rpts,radicals):
    rows=[]
    for k,(i,l) in enumerate(EDGES):
        d=(pts[i][0]-rpts[l][0],pts[i][1]-rpts[l][1])
        row=[0]*(3*(N+J))
        row[2*i],row[2*i+1]=2*d[0],2*d[1]
        row[2*N+2*l],row[2*N+2*l+1]=-2*d[0],-2*d[1]
        row[2*(N+J)+i]=-radicals[k]
        row[2*(N+J)+N+l]=radicals[k]
        rows.append(row)
    return s.Matrix(rows)


class JointSystemTests(unittest.TestCase):
    def setUp(self):
        self.rng=random.Random(20261010)

    def q(self):return Q(self.rng.randint(-30,30),self.rng.randint(1,7))

    def _distinct(self,pt,count):
        # Grid vertices are distinct by the generic-grid lemma.
        out=[]
        while len(out)<count:
            p=pt()
            if p not in out:out.append(p)
        return out

    def generic(self):
        pt=lambda:(self.q(),self.q())
        pts=self._distinct(pt,N+J);return pts[:N],pts[N:]

    def circle(self):
        def pt():
            t=self.q();return (5*(1-t*t)/(1+t*t),10*t/(1+t*t))
        pts=self._distinct(pt,N+J);return pts[:N],pts[N:]

    def test_tree_is_spanning_and_avoids_core(self):
        self.assertEqual(len(TREE),N+J-1)
        self.assertTrue(TREE.isdisjoint(CORE))
        adj={('p',i):set() for i in range(N)};adj.update({('r',l):set() for l in range(J)})
        for (i,l) in TREE:adj[('p',i)].add(('r',l));adj[('r',l)].add(('p',i))
        seen={('p',0)};stack=[('p',0)]
        while stack:
            v=stack.pop()
            for w in adj[v]:
                if w not in seen:seen.add(w);stack.append(w)
        self.assertEqual(len(seen),N+J)
        # Each vertex outside the core keeps at least two non-tree edges.
        for v in [('p',4),('p',5),('r',3),('r',4)]:
            deg=sum(1 for (i,l) in EDGES if (i,l) not in TREE and (('p',i)==v or ('r',l)==v))
            self.assertGreaterEqual(deg,2)

    def test_rank_generic_positions(self):
        for _ in range(3):
            pts,rpts=self.generic()
            self.assertEqual(rigidity(pts,rpts,CORE).rank(),11)
            self.assertEqual(rigidity(pts,rpts,EDGES).rank(),19)
            self.assertEqual(rigidity(pts,rpts,[e for e in EDGES if e not in TREE]).rank(),19)

    def test_rank_on_a_conic(self):
        for _ in range(3):
            pts,rpts=self.circle()
            self.assertEqual(rigidity(pts,rpts,CORE).rank(),10)
            self.assertEqual(rigidity(pts,rpts,EDGES).rank(),18)
            self.assertEqual(rigidity(pts,rpts,[e for e in EDGES if e not in TREE]).rank(),18)

    def test_conic_motion_is_a_flex(self):
        # Points on x^T S x + 2 b.x + c = 0 move by H = Sp + b, U = -Sr - b.
        x1,y1,x2,y2,s11,s12,s22,b1,b2,c=s.symbols('x1 y1 x2 y2 s11 s12 s22 b1 b2 c')
        S=s.Matrix([[s11,s12],[s12,s22]]);b=s.Matrix([b1,b2])
        p=s.Matrix([x1,y1]);r=s.Matrix([x2,y2])
        G=lambda v:(v.T*S*v)[0]+2*b.dot(v)+c
        H=S*p+b;U=-S*r-b
        self.assertEqual(s.expand((p-r).dot(H-U)-(G(p)-G(r))),0)

    def test_attachment_adds_rank_two_unless_collinear(self):
        u1=(Q(0),Q(0));u2=(Q(1),Q(0))
        for v,expected in [((Q(0),Q(1)),2),((Q(2),Q(0)),1)]:
            m=s.Matrix([[v[0]-u1[0],v[1]-u1[1]],[v[0]-u2[0],v[1]-u2[1]]])
            self.assertEqual(m.rank(),expected)

    def _kernel_has_constant_scalars(self,pts,rpts,expected_rank):
        radicals=[self.q() or Q(1) for _ in EDGES]
        M=joint_system(pts,rpts,radicals)
        self.assertEqual(M.rank(),expected_rank)
        kernel=M.nullspace()
        self.assertEqual(len(kernel),3*(N+J)-expected_rank)
        for v in kernel:
            scalars=v[2*(N+J):]
            self.assertTrue(all(a==scalars[0] for a in scalars))

    def test_joint_kernel_generic(self):
        pts,rpts=self.generic()
        self._kernel_has_constant_scalars(pts,rpts,10+19)

    def test_joint_kernel_on_a_conic(self):
        pts,rpts=self.circle()
        self._kernel_has_constant_scalars(pts,rpts,10+18)

    def test_smaller_grids_fail_the_count(self):
        self.assertGreaterEqual(6*5,3*(6+5)-4)
        self.assertGreaterEqual((6-3)*(5-3),5)
        for n,j in [(7,4),(5,5),(8,3),(4,7)]:
            self.assertLess(n*j,3*(n+j)-4)
            self.assertLess((n-3)*(j-3),5)
        # Thirty is the least product with (n-3)(j-3) >= 5.
        self.assertEqual(min(n*j for n in range(1,40) for j in range(1,40) if (n-3)*(j-3)>=5 and n>3 and j>3),30)
        # With 7x4 the kernel is forced to contain a vector with nonconstant scalars.
        n,j=7,4
        pts=[(self.q(),self.q()) for _ in range(n)];rpts=[(self.q(),self.q()) for _ in range(j)]
        rows=[]
        for (i,l) in [(i,l) for i in range(n) for l in range(j)]:
            d=(pts[i][0]-rpts[l][0],pts[i][1]-rpts[l][1]);row=[0]*(3*(n+j))
            row[2*i],row[2*i+1]=2*d[0],2*d[1];row[2*n+2*l],row[2*n+2*l+1]=-2*d[0],-2*d[1]
            R=self.q() or Q(1);row[2*(n+j)+i]=-R;row[2*(n+j)+n+l]=R;rows.append(row)
        kernel=s.Matrix(rows).nullspace()
        self.assertTrue(any(len(set(v[2*(n+j):]))>1 for v in kernel))


if __name__=='__main__':unittest.main()
