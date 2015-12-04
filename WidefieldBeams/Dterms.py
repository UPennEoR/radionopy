from sympy import *
init_printing(use_unicode=True)
from sympy.physics.quantum import TensorProduct

# Elements of the widefield antenna matrix
Dxa,Dxb,Dya,Dyb = symbols('D_xa D_xb D_ya D_yb')
Ax,Ay = symbols('A_x A_y')
# Coherency vector in stokes basis
I,Q,U,V = symbols('I Q U V')
Ex,Ey = symbols('E_x,E_y')

E = Matrix([Ex,Ey])

Da = Matrix([[1,Dya],[Dxa,1]])
Db = Matrix([[1,Dyb],[Dxb,1]])

Econj = Matrix([conjugate(Ex),conjugate(Ey)])
# Go from Stokes to coherency, drop the usual factor of 1/2
S = Matrix([[1,1,0,0],[0,0,1,1j],[0,0,1,-1j],[1,-1,0,0]])
# Go from coherency to Stokes
Sinv = Matrix.inv(S)

V_S = Matrix([I, Q, U, V])
V_C = TensorProduct(E,Econj)

D_C = simplify(TensorProduct(Da,Db))
D_S = simplify(Sinv*D*S)

