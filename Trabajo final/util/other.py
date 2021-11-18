from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import mpld3
def deriv(y, t, N, beta, gamma):
    C, I, R = y
    dSdt = -beta * C * I / N
    dIdt = beta * C * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt
def plotsir(t, S, I, R):
  f, ax = plt.subplots(1,1,figsize=(10,4))
  ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='susceptibles')
  ax.plot(t, I, 'y', alpha=0.7, linewidth=2, label='INFECTADOS')
  ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='RECUPERADOS')
  ax.set_xlabel('Time (days)')
  ax.yaxis.set_tick_params(length=0)
  ax.xaxis.set_tick_params(length=0)
  ax.grid(b=True, which='major', c='w', lw=2, ls='-')
  legend = ax.legend()
  legend.get_frame().set_alpha(0.5)
  for spine in ('top', 'right', 'bottom', 'left'):
      ax.spines[spine].set_visible(False)
  plt.show()
N = 479000
beta = 0.06  # infected person infects 1 other person per day
gamma = 0.021
S0, I0, R0 = 999, 1, 0
t = np.linspace(0, 49, 50) # Grid of time points (in days)
y0 = S0, I0, R0 # Initial conditions vector
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = ret.T
plotsir(t, S, I, R)