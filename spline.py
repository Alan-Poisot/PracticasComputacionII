from sympy import *
from sympy.abc import x
import numpy
import matplotlib.pyplot as plt
import numpy as np
import math as m
import pandas as pd

class Graph():

    # Genera un gráfico de dispersión
    def plotScatter(self, X, Y, legend):
        plt.scatter(X, Y, label = legend)
        plt.title("Gráfica de número de elecciones nacionales por año")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.legend(loc='best')
        
    
    # Genera un gráfico de línea
    def plotLine(self, func,a ,b, legend):
        X = np.linspace(a, b, 100)
        Y = np.zeros_like(X)

        for i in range(len(X)):
            Y[i] = func.subs(x, X[i])
        
        plt.plot(X, Y)
        plt.legend(loc='best')
    
    # Despliega el gráfico (útil para ) 
    def displayPlot(self):
        plt.legend(loc='best')
        #plt.ylim([0, max(m)+15])
        plt.savefig("Graph",dpi=250)
        plt.show()

class ImportData():

    def dataCuration(self, df):
        df = df.dropna()
        return df
    
    def importDataFrame(self):
        #Leer el dataframe number-of-direct-national-elections.csv
        df = pd.read_csv("number-of-direct-national-elections.csv") #Recuperado de https://ourworldindata.org/grapher/number-of-direct-national-elections
        print("Datos:\n",df,"\n")
        X = df.Year.to_numpy()
        Y = df.Count.to_numpy()
        return X, Y

class SplineInterpolation():

    def verifyContinuity(self, fun1, fun2, fun3):
        # Calcular primer derivada de cada función
        der1 = fun1.diff()
        der2 = fun2.diff()
        der3 = fun3.diff()

        # Graficar funciones sobre un intervalo
        objG = Graph()
        objG.plotLine(fun1, -2, 0, "1st")
        objG.plotLine(fun2, 0, 1, "2nd")
        objG.plotLine(fun3, 1, 2, "3rd")
        objG.displayPlot()

        # Graficar primeras derivadas
        objG.plotLine(der1, -2, 0, "1st der")
        objG.plotLine(der2, 0, 1, "1st der")
        objG.plotLine(der3, 1, 2, "1st der")
        objG.displayPlot()

        objG.displayPlot()

    def createAugmentedMatrix(self, X, Y):
        
        #Definir vector de nodos y resultados
        knots = list()
        results = list()

        # Generación de vector de knots
        for i in range(len(X)):
            if (i == 0 or i == (len(X)-1)):
                knots.append(X[i])
                results.append(Y[i])
            else: 
                knots.append(X[i]) #TODO: verificar inserción doble
                knots.append(X[i])
                results.append(Y[i])
                results.append(Y[i])
        print("Knots vector...")
        print(knots)
        print("Results vector ...")
        print(results)

        # Definición de listas para variables simbólicas
        A = []
        B = []
        C = []

        #Producir variables simbólicas 
        #a1,...,a3, b1,...,b3, c1,...,c3 
        """for i in range(len(X)-1):
            A.append(symbols('a' +str(i+1)))
            B.append(symbols('b' +str(i+1)))
            C.append(symbols('c' +str(i+1)))
        print("Variables simbólicas")
        print(A)
        print(B)
        print(C)"""
        
        # Generación de 3*n ecuaciones
        n = len(X)-1
        
        # Matriz para almacenar sistema de ecuaciones 3*n x 3*n+1  
        eqs = np.zeros(shape=(3*n,(3*n)+1))
        
        coefIndex = 0

        #Producir las primeras 2*n ecuaciones 
        for i in range(len(knots)):
            # Knots externos
            if (i == 0 or i == len(knots)):
                # Equations value
                eqs[i][coefIndex] =  knots[i]**2
                eqs[i][coefIndex+1] =  knots[i]
                eqs[i][coefIndex+2] = 1
                #results in the last column
                eqs[i][-1] = results[i]
            # Knots internos
            else:
                eqs[i][coefIndex] =  knots[i]**2
                eqs[i][coefIndex+1] =  knots[i]
                eqs[i][coefIndex+2] = 1
                #results in the last column
                eqs[i][-1] = results[i]
                # incrementar en uno cada par de elementos
                coefIndex += i%2 * 3
            
        # Add n-1 equations
        coefIndex = 0
        for i in range(0,n-1):
            eqs[2*n+i][coefIndex] =  2*X[i+1]
            eqs[2*n+i][coefIndex+1] =  1
            eqs[2*n+i][coefIndex+3] = -2*X[i+1]
            eqs[2*n+i][coefIndex+4] =  -1
            # incrementar en 3 el indice de coeficiente
            coefIndex += 3
        
        # add final equation a1=1
        eqs[3*n-1][0] = 1

        return eqs



class GaussJordan:
    """
        Método para intercambiar dos filas de la matriz M
        Entradas: indices de la primer y segunda fila, matriz 
        Salidas: Matriz M modificada
    """
    def intercambiarFilas(self, index1, index2, M): 
        for i in range(len(M[index1])):
            temp=M[index1,i]
            M[index1,i]=M[index2,i]
            M[index2,i]=temp
        return M

    def multiplicarFila(self, k, fila, colInicial, M):
        M[fila]=M[fila]*k
        return M

    """
        Método para restar dos filas de la matriz
        Entradas: indices de  filas 1 y 2, y matriz
        Salida: Matriz M modficada
    """
    def restarFilas(self, f1, f2, M):
        for i in range(len(M)):
            M[f1,i] = M[f2,i] -M[f1,i]
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
            if M[i][col] != 0:
                return i
                maxNum = abs(M[i][col])
        return indiceFila

    """
        Método para efectuar la eliminación Gaussiana
        Entradas: Número de filas y columnas, matriz
        Salida: Matriz M modificada
    """
    def eliminacionGaussiana(self, f, c, M):
        # Definición de variables
        indicePiv = -1
        M = self.intercambiarFilas(0, -1, M) # Se sabe que la última fila ya está resuelta
        for i in range(f):
            pivote = M[i][i]
            if pivote == 0:
                indicePiv = self.buscarPivote(f, i, M) # Se implementa pivoteo parcial
                if indicePiv == -1:
                    print("El sistema no tiene solución",i)
                    exit(0)
                else:
                    M = self.intercambiarFilas(indicePiv, i, M)
                    pivote = M[i][i]
            
            for j in range(f): # Realizar la eliminación de los elementos en la columna i
                if (M[j][i] != 0) and (i!=j):
                    q = M[j][i]/pivote    # Multiplicador para la eliminación
                    for k in range(f+1):
                        M[j][k] = M[j][k] - q*M[i][k]

                    #M = self.multiplicarFila(k, j, i, M)
                    #M = self.restarFilas(j, i, M)
        for i in range(f):
            M[i][-1]/=M[i][i]
            M[i][i]/=M[i][i]



def main():
    # Object to spline class
    spline = SplineInterpolation()
    objG = Graph()
    importObj = ImportData()
    

    # Importar datos desde un CSV utilizando pandas 
    X, Y = importObj.importDataFrame()
    #X=[3,4.5,7,9]
    #Y=[2.5,1,3,0.5]
    print("Año:", X)
    print("Elecciones: ", Y)

    #Graficar Puntos originales
    objG.plotScatter(X, Y, "Puntos originales")
    #objG.displayPlot()

    eqs = spline.createAugmentedMatrix(X,Y)
    print("Ecuaciones:\n", eqs)
    f=(len(X)-1)*3
    c=f+1
    #eqs=np.linalg.inv(e)
    objGJ = GaussJordan() # objGJ es una variable de la clase "GaussJordan"
    objGJ.eliminacionGaussiana(f, c, eqs) 
    #objGJ.GJ(f,c,eqs) # Se usa el algoritmo de Gauss-Jordan sobre la matriz
    
    
    print("El resultado del método de Gauss-Jordan es:\n",eqs)
    print("Por lo tanto, los valores de cada coeficiente son:")
    for i in range(len(X)-1):
        print("a"+str(i+1)+" = "+str(eqs[i*3,-1]))
        print("b"+str(i+1)+" = "+str(eqs[i*3+1,-1]))
        print("c"+str(i+1)+" = "+str(eqs[i*3+2,-1]))

    n=float(input("Valor de x para predecir: "))
    for i in range(len(X)-1):
        fun=eqs[i*3,-1]*x**2+eqs[i*3+1,-1]*x+eqs[i*3+2,-1]
        if (X[i]<=n) and (X[i+1]>=n):
            funInt=fun
        objG.plotLine(fun, X[i], X[i+1], "Spline")

    yn=funInt.subs(x,n)
    objG.plotScatter(n, yn , "Punto interpolado")
    
    objG.displayPlot()


if __name__ == "__main__":
    main()