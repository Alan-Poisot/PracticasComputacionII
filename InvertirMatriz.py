import numpy as np
import sympy

class GaussJordan:
    
    """
        Método para intercambiar dos filas de la matriz M
        Entradas: indices de la primer y segunda fila 
        Salidas: Matriz M modificada
    """
    def intercambiarFilas(self, index1, index2, M): 
        M[index1],M[index2]=M[index2],M[index1]
        return M
   
    def multiplicarFila(self, k, fila, colInicial, M):
        M[fila] = k * M[fila]
        return M

    """
        Método empleado para realizar la eliminación resta dos filas
        Entradas: indices de  filas 1 y 2, columna inicial y Matriz
        Salida: Matriz M modficada
    """
    def restarFilas(self, f1, f2, colInicial, M):
        M[f1] =  M[f2] - M[f1]
        return M 

    """
        Método para buscar un elemento pivote, se implementa pivoteo parcial
        Entradas: filas, columna actual y matriz
        Salidas: Matriz modificada
    """
    def buscarPivote(self, filas, col, M):
        indiceFila = -1
        maxNum = np.inf *-1
        for i in range(col+1, filas):
            if(M[i][col] > maxNum and M[i][col] != 0):
                indiceFila = i
                maxNum = M[i][col]
        return indiceFila

    #def calcularInversa(self, MAug, f, c): #Función no es necesaria
        #self.GaussJordan()

    """
        Método para efectuar la eliminación Gaussiana
        Entradas: Número de filas y columnas, Matriz
        Salida: Matriz modificada
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
                    M = self.restarFilas(j, i, i, M)
        #print("Matriz resultante EG: \n", M)

    def GJ(self,f,c,M):
        for i in range(f):
            pivote=M[f-1-i][f-1-i]
            M[f-1-i]=M[f-1-i]/pivote
            for j in range(f-1-i):
              M[j]=M[j]-M[j][f-1-i]*M[f-1-i]
        #print("Matriz resultante: \n", M)

"""
    Método para calcular la inversa de una matriz M
    Entradas: Filas, columnas y Matriz original
    Salida: Matriz inversa
"""
def calculateInverse(f, c, M,objG):
    # Definición de matriz identidad
    I = np.identity(f)

    #actualizar número de columnas
    c=c+f
    MAug = np.concatenate([M,I], axis=1)
    #print("Matriz aumentada: \n", MAug)
    # Generar matriz aumentada
    objG.eliminacionGaussiana(f,c,MAug)
    objG.GJ(f,c,MAug)
    print(MAug)
    print("Matriz inversa:")
    Minv=np.zeros_like(I)
    for i in range(f):
      Minv[i]=MAug[i][c-f:c]
      print(MAug[i][c-f:c]) #imprime las úlimas f filas, correspondientes a la matriz inversa
    return Minv
    
def main():

    # Definición de número de filas y columnas
    f = 3
    c = f+1 # +1 se debe a la columna de resultados
    
    # Inicializar una matriz de 
    M = np.random.randint(100, size = (f,c))

    print("matriz aleatoria:\n", M)
    
    # Creación de un objeto
    objG = GaussJordan() # objG es una variable de la clase "GaussJordan"
    #print("Realizando eliminación Gaussiana...\n")
    #objG.eliminacionGaussiana(f, c, M) # aplicar el método eliminacionGaussiana en objG
    #objG.GJ(f,c,M)
    print("Cálculo de soluciones y matriz inversa:\n")
    calculateInverse(f,c,M,objG)

    #aplicación para elevar matrices a potencias enteras
    print("\n--------------------------------------------------------\nAplicacion para elevar una matriz a una potencia entera\nSe usa la fórmula A^n=PD^nP^-1, donde P es la matriz formada por los eigenvectores\ny D es una matriz con los eigenvalores en su diagonal\n--------------------------------------------------------\n")
    objG = GaussJordan()
    c=f
    M = np.random.randint(100, size = (f,c))
    print("matriz aleatoria:\n", M)
    P= np.zeros_like(M)
    D= np.identity(f)
    eigv,eigvec=np.linalg.eig(M)
    print("Eigenvectores:\n",eigvec)
    

    P=np.transpose(eigvec)
    Pinv=P
    
    print("Matriz aumentada inversa:\n")
    Pinv=calculateInverse(f,c,Pinv,objG)
    #print("Matriz inversa\n",Pinv)
    n=0
    while n<1:
      n=int(input("Ingrese una potencia entera positiva para elevar la matriz: "))
    for i in range(f):
      D[i]=D[i]*(eigv[i]**n)
    R=P@D@Pinv
    print("El resultado es:\n",R)


if __name__ == "__main__":
    main()