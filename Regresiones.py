import sympy
from sympy.abc  import x
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class LeastSquares():

    def linealRegression(self, X, Y):
        if(len(X) == len(Y)):
            Xi = 0
            Yi = 0
            XY = 0
            Xi2 =0
            for i in range(len(X)):
                Xi += X[i]
                Yi += Y[i]
                XY += X[i]*Y[i]
                Xi2 += X[i]**2
            # calculate a1 
            slope = ((len(X) * XY ) - (Xi *Yi)) / ((len(X) * Xi2) - Xi**2)
            # calculate a0
            intercept = (Yi/len(X)) - (slope * (Xi/len(X)))
        else:
            print("Los vectores no tienen el mismo número de elementos...")
            exit(0)

        return intercept, slope
    
    def squareMeanError(self, X, Y, func):
        meanError = 0
        # 
        for i in range (len(X)):
            meanError += (abs(Y[i] - func.subs(x, X[i])))**2
        return meanError/len(X)
    
    def determinationCoef(self, X, Y, func,sme):
        ss = 0
        yMean=np.mean(Y)
        for i in range (len(X)):
            ss += (abs(Y[i] - yMean))**2
        
        return 1-(len(X)*sme/ss)

    # TODO: Calcular coeficiente de determinación
    #       Calc. coef. de correlación

    def nOrderReg(self,X,Y,n):
        f=n+1
        c=f+1
        M = np.zeros((f,c))
        Xn=np.zeros(2*n+1)
        XYn=np.zeros(n+1)
        for i in range(2*n+1):
            for j in range(len(X)):
                Xn[i]+=X[j]**i
        for i in range(n+1):
            for j in range(len(X)):
                XYn[i]+=Y[j]*X[j]**i
        # Se llena la matriz
        for i in range(f):
            for j in range(f):
                M[i][j]=Xn[i+j]
            M[i][-1]=XYn[i]

        objG = GaussJordan() # objG es una variable de la clase "GaussJordan"
        objG.eliminacionGaussiana(f, c, M) 
        objG.GJ(f,c,M) # Se usa el algoritmo de Gauss-Jordan sobre la matriz
        
        pol=0
        for i in range(f):
            pol+=M[i][-1]*x**i

        return pol


class Graph():
    # Genera un gráfico de dispersión
    def plotScatter(self, X, Y, lab):
        plt.scatter(X, Y, c="b", label = lab)
        plt.legend(loc='best')
        plt.title("Precio de acciones Vs Ganancias por acción (S&P 500)")
        plt.xlabel("Precio de acción ($)")
        plt.ylabel("Ganancias por acción ($)")
        plt.ylim([min(Y)*0.85-3, max(Y)*1.15+3])
    
    # Genera un gráfico de línea
    def plotLine(self, func,a ,b, leg):
        X = np.linspace(a, b, 100)
        Y = np.zeros_like(X)

        for i in range(len(X)):
            Y[i] = func.subs(x, X[i])
        
        plt.plot(X,Y, ls='--', label = leg )
        plt.legend(loc='best')
    # Despliega el gráfico
    def displayPlot(self):
        plt.savefig("Graph",dpi=250)
        plt.show()

class GaussJordan:
    
    """
        Método para intercambiar dos filas de la matriz M
        Entradas: indices de la primer y segunda fila, matriz 
        Salidas: Matriz M modificada
    """
    def intercambiarFilas(self, index1, index2, M): 
        M[index1],M[index2]=M[index2],M[index1]
        return M
   
    def multiplicarFila(self, k, fila, colInicial, M):
        M[fila] = k * M[fila]
        return M

    """
        Método para restar dos filas de la matriz
        Entradas: indices de  filas 1 y 2, y matriz
        Salida: Matriz M modficada
    """
    def restarFilas(self, f1, f2, M):
        M[f1] =  M[f2] - M[f1]
        return M 

    """
        Método para buscar un elemento pivote, se implementa pivoteo parcial
        Entradas: filas, columna actual y matriz
        Salidas: Matriz M modificada
    """
    def buscarPivote(self, filas, col, M):
        indiceFila = -1
        maxNum = np.inf *-1
        for i in range(col+1, filas):
            if(M[i][col] > maxNum and M[i][col] != 0):
                indiceFila = i
                maxNum = M[i][col]
        return indiceFila

    """
        Método para efectuar la eliminación Gaussiana
        Entradas: Número de filas y columnas, matriz
        Salida: Matriz M modificada
    """
    def eliminacionGaussiana(self, f, c, M):
        # Definición de variables
        indicePiv = -1
        
        for i in range(f):
            pivote = M[i][i]
            if pivote == 0:
                indicePiv = self.buscarPivote(f, i, M) # Se implementa pivoteo parcial
                #TODO: Implementar pivoteo completo
                if indicePiv == -1:
                    print("El sistema no tiene solución")
                    exit(0)
                else:
                    M = self.intercambiarFilas(indicePiv, i, M)
                    pivote = M[i][i]
            
            for j in range(i+1, f): # Realizar la eliminación de los elementos debajo del pivote
                if M[j][i] != 0:
                    k = pivote / M[j][i]    # Multiplicador para la eliminación
                    M = self.multiplicarFila(k, j, i, M)
                    M = self.restarFilas(j, i, M)
    """
        Método para realizar el algoritmo de Gauss-Jordan
        Entradas: Número de filas y columnas, matriz
        Salida: Matriz M modificada
    """
    def GJ(self,f,c,M):
        for i in range(f):
            pivote=M[f-1-i][f-1-i]
            M[f-1-i]=M[f-1-i]/pivote
            for j in range(f-1-i):
              M[j]=M[j]-M[j][f-1-i]*M[f-1-i]

def main():
    df = pd.read_csv("https://datahub.io/core/s-and-p-500-companies-financials/r/1.csv")
    X = df["Price"]
    Y = df["Earnings/Share"]
    X, Y = zip(*sorted(zip(X, Y))) # Se ordenan ambos arreglos juntos
    print(X,Y)
    order=[1,2,4,6,8,10]
    Reg=list()
    ls = LeastSquares()
    a0, a1 = ls.linealRegression(X, Y)

    Reg.append(a0 + a1 * x)

    graph = Graph()
    graph.plotScatter(X,Y, "Puntos originales")

    graph.plotLine(Reg[0], X[0], X[-1], "Regresión lineal")
    sme=ls.squareMeanError(X, Y, Reg[0])
    print("Error cuadrático medio y coeficiente de determinación grado 1: ",  sme,", ",ls.determinationCoef(X, Y, Reg[0],sme))
    
    for i in range(1,len(order)):
        Reg.append(ls.nOrderReg(X,Y,order[i]))
        graph.plotLine(Reg[i], X[0], X[-1], str("Grado " + str(order[i])))
        sme=ls.squareMeanError(X, Y, Reg[i])
        print("Error cuadrático medio y coeficiente de determinación grado "+str(order[i])+": ", sme,", ",ls.determinationCoef(X, Y, Reg[i],sme))

    graph.displayPlot()

if __name__ == "__main__":
    main()
