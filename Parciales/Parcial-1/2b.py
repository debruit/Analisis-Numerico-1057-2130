import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import approximate_taylor_polynomial


f = lambda x: (x-1)*np.log(x)


def taylor(funcion, valorX, x0, x1):
    x = np.linspace(x0, x1, num=1)
    plt.plot(x, f(valorX), label="Funcion Original")

    history = []
    for degree in np.arange(1, 15):
        n_taylor = approximate_taylor_polynomial(f, 0, degree, 1,)
        plt.plot(x, n_taylor(x), label=f"degree = {degree}")
        history.append(n_taylor)

    valoresY = []
    for xi in x:
        valoresY.append(funcion(xi))

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.0, shadow=True)
    plt.tight_layout()
    plt.axis([x0, x1, -10, 10])
    plt.show()
    
    return np.array(history)

polinomios = taylor(f, 1, -10, 10)