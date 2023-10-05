
def u_target(t):
     data = Expression("16*x[0]*(x[0]-1)*x[1]*(x[1]-1)*sin(pi*t)", degree=4)
     u_star = interpolate(data, V)
     return u_target

