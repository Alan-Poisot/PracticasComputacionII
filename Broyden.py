from sympy import *
import sympy
import numpy as np
import time
from numpy.linalg import inv

class Jacobian():
    def __init__(self, n):
        self.n = n;

    def generateJacobianMat(self, XSymb, F):
        # Generar una matriz jacobiana
        J = []
        for i in range (self.n):
            tmp = list()
            for j in range(self.n): 
                tmp.append(diff(F[i], XSymb[j]))
            J.append(tmp)
        return J

    def generateHessianMat(self, M, X, F):
        H = self.generateJacobianMat(M, X, F)
        H = self.generateJacobianMat(H, X, F)
        return H

    def evaluateF(self, Xsymb, F, X0):
        F1 = np.ndarray(self.n)
        for i in range(self.n):
            F1[i] = F[i].subs(dict(zip(Xsymb,X0)))
        return F1

    def evaluateJ(self, J, X0, Xsymb):
        JNum = np.zeros(shape=(self.n, self.n))
        for i in range(self.n):
            for j in range(self.n):
                JNum[i][j] = J[i][j].subs(dict(zip(Xsymb,X0)))
        return JNum


class NewtonSENL():
    def __init__(self, n):
        self.n = n

    def runNewtonMethod(self, tol, maxIter, X0, Xsymb, F, J):
        
        jac = Jacobian(self.n) 
        start = time.time()
        error = np.inf
        iter = 0

        while(error > tol and iter < maxIter):
            F1 = jac.evaluateF(Xsymb, F, X0)
            J1 = jac.evaluateJ(J, X0, Xsymb)

            # Calcular la inversa de J1 
            invJ = inv(J1)
            # Hacer el producto punto  F1 @ J1Inv para alcular y(0)
            Y0 =  -F1 @ invJ
            # X(0) + y(0) = x1
            X1=X0+Y0
            
            # calcular  error 
            error=0
            for i in range(len(X0)):
                error+=(X1[i]-X0[i])**2
            error=sqrt(error)
            # X0 = X1 
            X0=X1
            iter += 1

        # Vector de resultados
        print("Total elapsed time Newton method: ", time.time()-start, "secs") 
        return X0 

    def inv(self,n,M):
        mat=np.zeros((n,n))
        i=np.eye(n)
        for j in range(len(M)):
            mat[j]=M[j]
        mat=np.concatenate((mat, i), axis=1)

        #Gauss
        i=0
        p=True
        while i<n:
            if mat[i][i]!=0:
                for j in range(i+1,n+1):
                    mat[i][j]=mat[i][j]/mat[i][i]
                mat[i][i]=1
                for j in range(i+1,n):
                    temp=mat[j][i]
                    for k in range(i,n+1):
                        mat[j][k]-=mat[i][k]*temp
                i=i+1
            else:
                ex,ind=swap(mat,n,i)
                if ex==True:
                    mat[i],mat[ind]=mat[ind],mat[i]
                else:
                    print("El sistema no tiene solución")
                    i=n
                    p=False

        #Gauss-Jordan

        if(p):
            i=n-1
            while i>0:
                for j in range (0,i):
                    mat[j]-=mat[i]*mat[j][i]
                i=i-1

            return mat[:,n:2*n]
        exit(0)

class BroydenSENL():
    def __init__(self, n):
        self.n = n
    
    def runBroydenMethod(self, tol, maxIter, X0, Xsymb, F, J):
        
        # Objetos
        jac = Jacobian(self.n)
        
        start = time.time()

        # Paso 2 
        F0 = jac.evaluateF(Xsymb, F, X0)

        # Paso 3 
        A0Inv = jac.evaluateJ(J, X0, Xsymb) 
        A0Inv = inv(A0Inv)

        error = np.inf
        iter = 0
        # Paso 4
        X1 = X0 - (A0Inv @ F0)

        while(error > tol and iter < maxIter):
            # Paso 5
            F1 = jac.evaluateF(Xsymb, F, X1)
            
            # Paso 6
            y1 = F1 - F0

            # Paso 6.1 
            s1 = X1 - X0
            
            # Paso 7
            tmp = (s1.T@A0Inv)@y1+1
            w=s1-A0Inv@y1
            w2=s1.T@A0Inv
            # Paso 8
            A1Inv = A0Inv + ((1/tmp)*w@w2)

            # Paso 9
            X1 = X0 - (A1Inv @ F1)
            
            # Paso 10
            # Norma L2 para calcular el error
            error=0
            for i in range(len(X0)):    
                error += (X1[i] - X0[i])**2 

            A0Inv = A1Inv
            X0 = X1
            F0 = F1
            iter +=1

        print("Total elapsed time Broyden method: ", time.time()-start, "secs") 
        return X1
    
    '''
    Método para calcular la matriz inversa
    Entradas: Tamaño de la matriz y una matriz cuadrada
    Salida: Matriz inversa (si existe)
    '''
    def inv(self,n,M):
        mat=np.zeros((n,n))
        i=np.eye(n)
        for j in range(len(M)):
            mat[j]=M[j]
        mat=np.concatenate((mat, i), axis=1)

        #Gauss
        i=0
        p=True
        while i<n:
            if mat[i][i]!=0:
                for j in range(i+1,n+1):
                    mat[i][j]=mat[i][j]/mat[i][i]
                mat[i][i]=1
                for j in range(i+1,n):
                    temp=mat[j][i]
                    for k in range(i,n+1):
                        mat[j][k]-=mat[i][k]*temp
                i=i+1
            else:
                ex,ind=swap(mat,n,i)
                if ex==True:
                    mat[i],mat[ind]=mat[ind],mat[i]
                else:
                    print("El sistema no tiene solución")
                    i=n
                    p=False

        #Gauss-Jordan

        if(p):
            i=n-1
            while i>0:
                for j in range (0,i):
                    mat[j]-=mat[i]*mat[j][i]
                i=i-1

            return mat[:,n:2*n]
        exit(0)
        

def swap(m,n,x):
    y=x+1
    for i in range(y,n):
        if m[i][x]!=0:
            return True, i;
    return False, 0;

def main():
    # number of eqs and vars
    n = 3

    # Definición de objetos
    jac = Jacobian(n)
    newton = NewtonSENL(n)
    broyden = BroydenSENL(n)

    # Lista para almacenar variables simbólicas
    Xsymb = list()

    # Agrgear variables simbólicas a una lista
    for i in range(n):
        Xsymb.append(symbols('x' +str(i+1)) )
    print("Variables simbólicas", Xsymb)
    
    # Definición de sistemas de ecuaciones
    F = list()
    F.append(3*Xsymb[0] - cos(Xsymb[1]*Xsymb[2]) - 3/2)
    F.append(Xsymb[0]**2 - 81*(Xsymb[1] + 0.3)**2 + sin(Xsymb[2]) + 1.06)
    F.append(np.e**(-Xsymb[0]*Xsymb[1]) + 20* Xsymb[2] + ((12 * np.pi) -3) /3)

    # Se instancia el método para generar la matriz Jac simbólica
    J = jac.generateJacobianMat(Xsymb, F)
    print("Jacobian Mat: ", J)

    # Método de Newton
    tol = 0.00000001
    maxIter = 10000

    # Definición del vector inicial
    X0 = np.ndarray(n)
    X0[0] = 0.1
    X0[1] = 0.1
    X0[2] = -0.1

    # Instanciar el método de newton y obtener resutados
    XN = newton.runNewtonMethod(tol, maxIter, X0, Xsymb, F, J)
    XB = broyden.runBroydenMethod(tol, maxIter, X0, Xsymb, F, J)

    print("Resultados con el método de Newton:", XN)
    
    print("Resultados con el método de Broyden:", XB)



if __name__ == "__main__":
    main()