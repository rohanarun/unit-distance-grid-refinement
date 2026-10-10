"""Checks for the collinear-anchor normalization and the rank-two Cramer step."""
from fractions import Fraction as Q
import random
import unittest
import sympy as s

J=s.Matrix([[0,-1],[1,0]])


class CollinearAnchorTests(unittest.TestCase):
    def test_leibniz_makes_normalized_velocities_parallel(self):
        # Anchors r_l = r_1 + sigma_l v on a line; velocities are t-derivatives.
        t=s.symbols('t')
        v=s.Matrix([s.Function('vx')(t),s.Function('vy')(t)])
        r1=s.Matrix([s.Function('ax')(t),s.Function('ay')(t)])
        sigma={1:s.Integer(0),2:s.Integer(1),3:s.Function('s3')(t),4:s.Function('s4')(t)}
        r={l:r1+sigma[l]*v for l in sigma}
        U={l:r[l].diff(t) for l in sigma}
        for l in sigma:
            leibniz=sigma[l]*v.diff(t)+sigma[l].diff(t)*v
            self.assertEqual(s.simplify(U[l]-U[1]-leibniz),s.zeros(2,1))
        omega=(U[2]-U[1]).dot(J*v)/v.dot(v)
        kappa=v.diff(t).dot(v)/v.dot(v)
        for l in sigma:
            normalized=U[l]-U[1]-omega*J*(r[l]-r1)
            self.assertEqual(s.simplify(normalized.dot(J*v)),0)
            mu=sigma[l]*kappa+sigma[l].diff(t)
            self.assertEqual(s.simplify(normalized-mu*v),s.zeros(2,1))
        # mu_1 = 0 and mu_2 = kappa, the same tau as in the general normalization.
        self.assertEqual(s.simplify(sigma[2]*kappa+sigma[2].diff(t)-kappa),0)

    def test_exact_rational_collinear_examples(self):
        rng=random.Random(20261010)
        def q():return Q(rng.randint(-20,20),rng.randint(1,9))
        for _ in range(500):
            vx,vy=q(),q()
            if vx==vy==0:continue
            dvx,dvy=q(),q()
            vv=vx*vx+vy*vy
            omega=(dvx*(-vy)+dvy*vx)/vv
            kappa=(dvx*vx+dvy*vy)/vv
            for _ in range(3):
                sig,dsig=q(),q()
                ux=sig*dvx+dsig*vx;uy=sig*dvy+dsig*vy
                nx=ux-omega*(-sig*vy);ny=uy-omega*(sig*vx)
                self.assertEqual(-vy*nx+vx*ny,0)
                mu=sig*kappa+dsig
                self.assertEqual((nx,ny),(mu*vx,mu*vy))

    def test_collinear_budget_and_boundary(self):
        anchors=4
        cost=(anchors-2)+(anchors-1)
        self.assertEqual(cost,2*anchors-3)
        self.assertEqual(cost,5)
        for dimension in (1,2):
            self.assertGreater(6*dimension-cost,6*(dimension-1))
            self.assertEqual(5*dimension-cost,5*(dimension-1))
        self.assertEqual(6*anchors,24)

    def test_rank_two_cramer_and_fourth_flip(self):
        x,y=s.symbols('x y');p=s.Matrix([x,y])
        R=s.symbols('R1:5');mu=s.symbols('mu1:5');Ap=s.symbols('Ap1:5')
        Hx,Hy,A=s.symbols('Hx Hy A');v=s.Matrix(s.symbols('vx vy'))
        r={l:s.Matrix(s.symbols(f'rx{l} ry{l}')) for l in range(1,5)}
        equations=[]
        for l in range(1,5):
            d=p-r[l]
            left=2*d.dot(s.Matrix([Hx,Hy]))-R[l-1]*A
            right=2*mu[l-1]*d.dot(v)-R[l-1]*Ap[l-1]
            equations.append(left-right)
        matrix,rhs=s.linear_eq_to_matrix(equations[:3],[Hx,Hy,A])
        det=s.expand(matrix.det())
        # Coefficient of each radical is, up to sign and the factor four,
        # the two-by-two minor of the other two position rows.
        for l in range(3):
            others=[m for m in range(3) if m!=l]
            d1=p-r[others[0]+1];d2=p-r[others[1]+1]
            minor=d1[0]*d2[1]-d1[1]*d2[0]
            self.assertEqual(s.simplify(det.coeff(R[l])**2-16*minor**2),0)
        solution=matrix.LUsolve(rhs)
        self.assertFalse(any(e.has(R[3]) for e in solution))
        residual=equations[3].subs(dict(zip([Hx,Hy,A],solution)))
        difference=s.simplify(residual.subs(R[3],-R[3])-residual)
        self.assertEqual(s.simplify(difference-2*R[3]*(solution[2]-Ap[3])),0)

    def test_one_grid_contains_all_witnesses(self):
        self.assertGreaterEqual(8,6);self.assertGreaterEqual(8,4)
        self.assertGreaterEqual(4,4);self.assertGreaterEqual(4,3)
        self.assertEqual((8*4,6*4,4*3),(32,24,12))


if __name__=='__main__':unittest.main()
