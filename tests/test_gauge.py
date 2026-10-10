"""Symbolic and exact-rational checks for the displayed gauge proof."""
from fractions import Fraction as Q
import random
import unittest
import sympy as s


class GaugeProofTests(unittest.TestCase):
    def test_symbolic_invariance(self):
        x,y,rx,ry,hx,hy,ux,uy,tx,ty,w,c=s.symbols('x y rx ry hx hy ux uy tx ty w c')
        p=s.Matrix([x,y]);r=s.Matrix([rx,ry]);H=s.Matrix([hx,hy]);U=s.Matrix([ux,uy]);T=s.Matrix([tx,ty]);J=s.Matrix([[0,-1],[1,0]])
        before=(p-r).dot(H-U)
        after=(p-r).dot((H-T-w*J*p)-(U-T-w*J*r))
        self.assertEqual(s.expand(after-before),0)
        # Scaling multiplies both sides of the system by the same factor.
        R,A,B=s.symbols('R A B')
        scaled=(p-r).dot(c*(H-T-w*J*p)-c*(U-T-w*J*r))
        self.assertEqual(s.expand(scaled-c*before),0)
        self.assertEqual(s.expand(R*(c*A-c*B)-c*R*(A-B)),0)

    def test_second_velocity_normalization(self):
        dx,dy,vx,vy=s.symbols('dx dy vx vy')
        d=s.Matrix([dx,dy]);v=s.Matrix([vx,vy]);J=s.Matrix([[0,-1],[1,0]])
        omega=v.dot(J*d)/d.dot(d)
        normalized=v-omega*J*d
        self.assertEqual(s.simplify(normalized.dot(J*d)),0)
        tau=v.dot(d)/d.dot(d)
        for component in normalized-tau*d:
            self.assertEqual(s.simplify(component),0)
        # After scaling by 1/tau the second velocity is exactly d.
        for component in s.simplify(normalized/tau-d):
            self.assertEqual(component,0)

    def test_exact_rational_examples(self):
        rng=random.Random(20261009)
        for _ in range(1000):
            dx,dy,vx,vy=(Q(rng.randint(-20,20),rng.randint(1,9)) for _ in range(4))
            if dx==dy==0:continue
            omega=(-dy*vx+dx*vy)/(dx*dx+dy*dy)
            normalized=(vx+omega*dy,vy-omega*dx)
            self.assertEqual(-dy*normalized[0]+dx*normalized[1],0)
            tau=(vx*dx+vy*dy)/(dx*dx+dy*dy)
            if tau!=0:
                self.assertEqual((normalized[0]/tau,normalized[1]/tau),(dx,dy))
            else:
                self.assertEqual(normalized,(0,0))

    def test_generator_budget_and_boundary(self):
        anchors=4
        right_velocity_data=2*(anchors-2)
        scalar_data=anchors-1
        cost=right_velocity_data+scalar_data
        self.assertEqual(cost,3*anchors-5)
        self.assertEqual(cost,7)
        for dimension in (1,2):
            self.assertGreater(8*dimension-cost,8*(dimension-1))
            # N=7 only gives equality; that count cannot prove survival.
            self.assertEqual(7*dimension-cost,7*(dimension-1))
        self.assertEqual(8*anchors,32)

    def test_wrong_rotation_sign_rejected(self):
        d=s.Matrix([1,0]);v=s.Matrix([0,1]);J=s.Matrix([[0,-1],[1,0]])
        omega=v.dot(J*d)/d.dot(d)
        self.assertNotEqual((v+omega*J*d).dot(J*d),0)

    def test_dilation_is_not_a_symmetry(self):
        # The symmetry group is five-dimensional: dilations change the
        # left side by 2*mu*L, which is not separable as R*(b_i-b'_l).
        x,y,rx,ry,hx,hy,ux,uy,mu=s.symbols('x y rx ry hx hy ux uy mu')
        p=s.Matrix([x,y]);r=s.Matrix([rx,ry]);H=s.Matrix([hx,hy]);U=s.Matrix([ux,uy])
        before=(p-r).dot(H-U)
        after=(p-r).dot((H+mu*p)-(U+mu*r))
        self.assertEqual(s.expand(after-before-mu*(p-r).dot(p-r)),0)
        self.assertNotEqual(s.expand(after-before),0)

    def test_line_case_cramer_and_sign_flip(self):
        s1,s2,s3,R1,R2,R3,u1,u2,u3,A1,A2,A3=s.symbols('s1 s2 s3 R1 R2 R3 u1 u2 u3 A1 A2 A3')
        matrix=s.Matrix([[2*s1,-R1],[2*s2,-R2]])
        det=s.expand(matrix.det())
        self.assertEqual(det.coeff(R1),2*s2)
        self.assertEqual(det.coeff(R2),-2*s1)
        rhs=s.Matrix([2*u1-R1*A1,2*u2-R2*A2])
        h,A=matrix.inv()*rhs
        self.assertEqual(s.simplify(matrix*s.Matrix([h,A])-rhs),s.zeros(2,1))
        residual=2*s3*h-R3*A-(2*u3-R3*A3)
        difference=s.expand(residual.subs(R3,-R3)-residual)
        self.assertEqual(s.simplify(difference-2*R3*(A-A3)),0)
        self.assertFalse(h.has(R3) or A.has(R3))
        # Dependent radicals can make the determinant vanish: independence
        # is essential, not an omitted hypothesis of the smaller witness.
        T=s.symbols('T')
        self.assertEqual(s.expand(det.subs({R1:s1*T,R2:s2*T})),0)

    def test_line_case_survival_budget(self):
        anchors=3
        cost=(anchors-2)+(anchors-1)
        self.assertEqual(cost,2*anchors-3)
        self.assertEqual(cost,3)
        self.assertGreater(4-cost,0)
        self.assertEqual(3-cost,0)
        self.assertEqual(4*anchors,12)


if __name__=='__main__':unittest.main()
