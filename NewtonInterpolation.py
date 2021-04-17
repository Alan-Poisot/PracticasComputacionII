import matplotlib.pyplot as plt
from sympy.abc import x
import numpy as np
import pandas as pd
import sympy


class Interpolation():

    # Generar una interpolación lineal
    def linearInterpolation(self, X, Y):
        # b0 = intecepto, b1 = pendiente
        b0 = Y[0]
        b1 = (Y[1] - Y[0]) / (abs(X[1]) - abs(X[0]))
        #print("b1", b1)
        return b0 + b1*(x - X[0])

    # 
    def tmpLinearInterpolation(self, X, Y):
        # b0 = intecepto, b1 = pendiente
        b0 = Y[0]
        b1 = (Y[1] - Y[0]) / (abs(X[1]) - abs(X[0]))

        return b0, b1

    def cuadraticInterpolation(self, X, Y):
        if (len(X) == 3 and len(Y) == 3):   
            b0, b1 = self.tmpLinearInterpolation(X, Y)
            b2 = (((Y[2] - Y[1]) / (X[2] - X[1])) - b1 ) / (X[2]-X[0])
            return b0 + b1*(x-X[0]) + b2*(x -X[0]) *(x-X[1])

    def interpolationN(self,n,X,Y):
        b = Y
        for i in range(1,n):
            for j in range(1,i):
                b[i] = (b[j] - b[i])/(X[j] - X[i])
        
        fun=b[-1]
        for i in range(len(b)-1):
            fun=fun*(x-X[-2-i])+b[-2-i]
            
        return fun

        


class Graph():
    # Genera un gráfico de dispersión
    def plotScatter(self, X, Y, lab):
        plt.scatter(X, Y, label = lab)
        plt.title("Cantidad global de elecciones nacionales anuales")
        plt.xlabel("Año")
        plt.ylabel("Elecciones")
    
    # Genera un gráfico de línea
    def plotLine(self, func,a ,b, leg,col):
        plt.title("Cantidad global de elecciones nacionales anuales")
        plt.xlabel("Año")
        plt.ylabel("Elecciones")
        X = np.linspace(a, b, 100)
        Y = np.zeros_like(X)

        for i in range(len(X)):
            Y[i] = func.subs(x, X[i])
        
        plt.plot(X,Y,  c = col, ls='--')
        
    # Despliega el gráfico (útil para ) 
    def displayPlot(self,m):
        plt.legend(loc='best')
        plt.ylim([0, max(m)+15])
        plt.savefig("Graph",dpi=250)
        plt.show()
        


def main():
    #Leer el dataframe number-of-direct-national-elections.csv
    df = pd.read_csv("number-of-direct-national-elections.csv") #Recuperado de https://ourworldindata.org/grapher/number-of-direct-national-elections
    print("Datos:\n",df,"\n")
    X = df.Year.to_numpy()
    Y = df.Count.to_numpy()

    objG = Graph()
    objG.plotScatter(X,Y, "Datos originales")

    objInt = Interpolation()


    # Calcular interpolación cuadrática

    funcList = list()  
    for i in range(len(X) - 2):
        #print(X[i:i+2])
        funcList.append(objInt.cuadraticInterpolation(X[i:i+3], Y[i:i+3]))

    # Interpolación de grado n
    grado=len(X)-1
    interN=objInt.interpolationN(grado,X,Y)
    print(sympy.simplify(interN))
    
    n=2000 #Se sabe que al año 2000 le corresponden 76 elecciones, se usa este número para probar el error de la interpolación
    if n<X[0]:
        objG.plotScatter(n,funcList[0],"Punto predicho")
    if n>X[-1]:
        objG.plotScatter(n,funcList[-1],"Punto predicho")
    for i in range(len(X)-2):
        if (n>X[i]) and (n<X[i+2]):
            yn=funcList[i].subs(x,n)
            objG.plotScatter(n,yn,"Punto predicho cuad.")
            print("Error de interpolación cuadrática= ",yn-76)
            break
    yn2=interN.subs(x,n)-76
    objG.plotScatter(n,yn2,"Punto predicho n")
    print("Error de interpolación de grado n= ",yn2-76)

    for i in range(len(X)-2):
        objG.plotLine(funcList[i], X[i], X[i+2], "Interpolación cuadrática","g")

    objG.plotLine(interN, X[0],X[-1], "Interpolación cuadrática","y")
    objG.displayPlot(Y)


if __name__ == "__main__":
    main()