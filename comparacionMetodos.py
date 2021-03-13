import numpy as np
from sympy.abc import x, y
from sympy import *
import time
import math
import matplotlib.pyplot as plt


class NumericalMethods:
    
    def newtonMethod(self, tol, maxIter, x0, func):
        raiz = np.inf
        # Calcular derivada de la función original
        firstDer = diff(func, x)

        # Inicializar error e iteraciones
        error = np.inf
        iter =0

        while iter < maxIter:
            # Validación de división entre 0
            if  firstDer.subs(x,x0)!=0:
                X1 = x0 - (func.subs(x, x0)/ firstDer.subs(x,x0))
                x0 = X1 
            # calcular error
            if abs(func.subs(x,X1))<=tol:
                return X1
            iter += 1
        
        return raiz

    def secantMethod(self, tol, maxIter, x_1, x0, func):
        raiz = np.inf
        # Inicializar error e iteraciones
        error = np.inf
        iter =0

        while(error > tol and iter < maxIter):
            #print("Iter: ", iter )
            fx0 = func.subs(x, x0)
            x1 = x0 - ( (fx0  * (x0 - x_1)) / (fx0 - func.subs(x, x_1)))
            error = abs(x1 - x0)
            #print("Error: ", float(error))
            x_1 = x0
            x0 = x1
            iter += 1
        if error<=tol:
            raiz = x0

        return raiz
    
    def bisectionMethod(self, a, b, tol, maxIter, func):
        fa=func.subs(x,a)
        fb=func.subs(x,b)
        if (fa*fb) < 0:
            #print("Existe un cambio de signo")
            fa=fa
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
    
    def bdMethod(self, a, b, tol, maxIter, func):
        root = np.inf
        error = np.inf
        iter = 0
        fa=func.subs(x, a)
        fb=func.subs(x, b)
        if  fa*fb >= 0:
            print("No existe una raíz en el intervalo proporcionado...")
            exit(0)
        # Validar valores
        if abs(fa) < abs(fb):
            # Cambiar los valores de a por b y viceversa
            a,b = b,a
            
        c = a
        m=True
        fs=np.inf
        while iter < maxIter:
            fa=func.subs(x, a)
            fb=func.subs(x, b)
            fc= func.subs(x, c)
            if (fa !=fc) and (fb != fc):
                # Interpolación cuadrática
                s = ((a*fb*fc) / ((fa-fb)*(fa-fc))) + ((b*fa*fc) / ((fb-fa)*(fb-fc))) + ((c*fa*fb) / ((fc-fa)*(fc-fb)))
            else: 
                # Método de la secante
                s = b - (fb*((b-a) / (fb-fa)))
            w=(3*a+b)/4
            if (s>max(w,b) or s<min(w,b)) or (m and abs(s-b)>=abs(b-c)/2) or ((not m) and abs(s-b)>=abs(c-d)/2) or (m and abs(b-c)<tol) or ((not m) and abs(c-d)<tol):
                # Método de la bisección
                s=(a+b)/2
                m=True
            else:
                m=False
            
            fs=func.subs(x,s)
            d=c
            c=b
            if fa*fs<0:
                b=s
            else:
                a=s
            
            if abs(func.subs(x, a)) < (func.subs(x, b)):
                # Cambiar los valores de a por b y viceversa
                a,b = b,a
            
            if abs(fs)<tol:
                return s
            elif abs(fb)<tol:
                return b

            iter+=1
        
        return root

        
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
    maxIter = 10000000
    tol = 0.0000001
    #func = 2*x**2 - 3*x - 5
    #func = x**3-67*x**2 + 790*x + 3000
    func = x**2 - 101*x - 19800
    objG = Graph()
    linf=-100
    lsup=100

    objG.graph2DPlot(-100, 100, func)

    # Crear objetos para el método ingenuo
    objN = NaiveMethods()

    start = time.time()
    objN.incrementalSearch(linf, lsup, tol, func)
    end = time.time()
    print("Exhaustiva,", (end-start))

    # Objeto para métodos numéticos
    objNM = NumericalMethods()
    
    # Método de la bisección
    start = time.time()
    raiz1 = objNM.bisectionMethod(linf, lsup, tol, maxIter, func )
    end = time.time()
    print("Biseccion,", end - start)
    
    # Método de Newton-Raphson
    start = time.time()
    raiz2 = objNM.newtonMethod(tol, maxIter, 0, func)
    end = time.time()
    print("Newton-Raphson,", end - start)
    
    # Método de la secante
    start = time.time()
    raiz3 = objNM.secantMethod(tol, maxIter, linf, lsup, func )
    end = time.time()
    print("Secante,", end - start)
    
    # Método de Brent-Dekker
    start = time.time()
    raiz4 = objNM.bdMethod(linf, lsup, tol, maxIter, func )
    end = time.time()
    print("Brent-Dekker,", end - start)
    print("raiz encontrada por el método de la bisección: ",float(raiz1))
    print("raiz encontrada por el método de Newton-Raphson: ",float(raiz2))
    print("raiz encontrada por el método de la secante: ",float(raiz3))
    print("raiz encontrada por el método de Brent-Dekker: ",float(raiz4))


if __name__ == "__main__":
    main()