# Let's try out sympy
import sympy
from sympy import Matrix
from sympy import symbols

a1, a2, a3, b1, b2, b3, c1, c2, c3 = symbols('a1 a2 a3 b1 b2 b3 c1 c2 c3')
A = Matrix( [[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]])
B = Matrix( [[1,2,3],[4,5,6],[c1,c2,c3]])
Binv = B.inv()
val = Binv[0, 0].subs([(c1,0), (c2,1), (c3,5)])
val_numeric = val.evalf()
print(val_numeric)

# Turn it into a function suitable for repeat evaluation
Binv00eval = sympy.lambdify([c1, c2, c3], Binv[0, 0], "numpy")
print(Binv00eval(0, 1, 5))
