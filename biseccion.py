import numpy as np
from sympy.abc import x, y
from sympy import *
import time
import math
import matplotlib.pyplot as plt


class NumericalMethods:
    
    def bisectionMethod(self, a, b, tol, maxIter, func):
        fa=func.subs(x,a)
        fb=func.subs(x,b)
        if (fa*fb) < 0:
            print("Existe un cambio de signo")
        else: 
            print("No existen raices reales en el intervalo")
            exit(0) 
        error = np.inf
        iter=0

        while iter<maxIter:
            c = (a + b) /2
            fc=func.subs(x,c)
            if (abs(fc)<=tol) or (abs(a-b)<=tol):
                return c

            fa=func.subs(x,a)
            fb=func.subs(x,b)
            

            if fa*fc<0:
                b=c
            else:
                a=c
            
            iter+=1
        
        print("Se alcanzó la cantidad máxima de iteraciones, la aproximación de la raiz es x=",c)
        return c


class NaiveMethods:

    """
        Este método efectua una búsqueda incrementar sobre el intervalo [a,b]
        Entradas: Intervalo [a,b], función
        Salida: Lista con raices reales
    """
    def incrementalSearch(self, a, b, tol, func):
        X = np.linspace(a, b, 10000)
        #Y = np.zeros_like(X)
        
        for i in range(len(X)):
            if (abs(func.subs(x, X[i]))) <= tol:
                print("The real root is: ", X[i])
            
class Graph:

    def graph2DPlot(self, start, end, equation):
        X = np.arange(start, end)
        Y = np.zeros(len(X))
        for i in range(len(X)):
            Y[i] = equation.subs(x, X[i])
        plt.plot(X,Y)
     

def main():

    x0 = 0
    maxIter = 900000
    tol = 0.0000000001
    #func = 2*x**2 - 3*x - 5
    #func = x**3-67*x**2 + 790*x + 3000
    func = x**2 - 101*x - 19800
    objG = Graph()

    objG.graph2DPlot(-100, 100, func)

    # Crear objetos para el método ingenuo
    objN = NaiveMethods()

    start = time.time()
    objN.incrementalSearch(-100, 100, tol, func)
    end = time.time()
    elapsedNaive = end-start
    print("The total elapsed time of incremental serach was: ", (end-start), "secs")

    # Objeto para métodos numéticos
    objNM = NumericalMethods()

    start = time.time()
    raiz = objNM.bisectionMethod(-100, 100 , tol, maxIter, func )
    end = time.time()
    elapsedBisection = end-start
    print("Elapsed time Bisection Method: ", end - start)

    print("Naive/Bisection = ", elapsedNaive/elapsedBisection)

    if (raiz != np.inf):
        print("La raiz real es: ", float(raiz))
    else:
        print("No se encontró la raíz")


if __name__ == "__main__":
    main()