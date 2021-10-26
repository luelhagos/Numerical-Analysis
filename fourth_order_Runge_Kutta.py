


# Write code to solve the following system of ordinary differential equations using the Python function odeint.

# dx1/dt& =1/2x1
# dx2/dt = 1/2x1 -1/4x2
# dx3/dt = 1/4x2-1/6x3
# on [0,4]

# Subject to the initial conditions x1(0) = 1, x2(0) = 1, x3(0) = 1


def RK4OdeSys(f, c, t):
    Z = np.zeros((len(t), len(c)))
    Z[0] = c
    H = t[1] - t[0]
    for k in range(len(t)-1):
        k1 = np.array(f(Z[k], t[k]))
        k2 = np.array(f(Z[k] + 0.5*H*k1, t[k] + 0.5*H))
        k3 = np.array(f(Z[k] + 0.5*H*k2, t[k] + 0.5*H))
        k4 = np.array(f(Z[k] + H*k3, t[k] + H))
        Z[k+1] = Z[k] + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    return Z
  
  # time points
t = np.linspace(0, 4, N)

# function that returns dz/dt
def model(x, t):
    x1, x2, x3 = x
    dx1dt = -0.5*x1
    dx2dt = 0.5*x1 - 0.25*x2
    dx3dt = 0.25*x2 - (1/6)*x3
    return dx1dt, dx2dt, dx3dt

# Using Runge-Kutta method
Rx = RK4OdeSys(model, [1, 1, 1], t)
