import numpy as np
import pandas as pd
from sympy.abc import x
import pygame
from pygame.locals import *
from scipy.interpolate import CubicSpline


class ImportData():

    def dataCuration(self, df):
        df = df.dropna()
        return df

    def importDataFrame(self):
        # Leer el dataframe number-of-direct-national-elections.csv
        df = pd.read_csv(
            "number-of-direct-national-elections.csv")  # Recuperado de https://ourworldindata.org/grapher/number-of-direct-national-elections
        print("Datos:\n", df, "\n")
        X = df.Year.to_numpy()
        Y = df.Count.to_numpy()
        return X, Y


class SplineInterpolation():
    def createAugmentedMatrix(self, X, Y):
        # Definir vector de nodos y resultados
        knots = list()
        results = list()

        # Generación de vector de knots
        for i in range(len(X)):
            if i == 0 or i == (len(X) - 1):
                knots.append(X[i])
                results.append(Y[i])
            else:
                knots.append(X[i])
                knots.append(X[i])
                results.append(Y[i])
                results.append(Y[i])
        print("Knots vector...")
        print(knots)
        print("Results vector ...")
        print(results)

        # Generación de 3*n ecuaciones
        n = len(X) - 1

        # Matriz para almacenar sistema de ecuaciones 3*n x 3*n+1  
        eqs = np.zeros(shape=(3 * n, (3 * n) + 1))

        coefIndex = 0
        # Producir las primeras 2*n ecuaciones
        for i in range(len(knots)):
            # Knots externos
            if (i == 0 or i == len(knots)):
                # Equations value
                eqs[i][coefIndex] = knots[i] ** 2
                eqs[i][coefIndex + 1] = knots[i]
                eqs[i][coefIndex + 2] = 1
                # results in the last column
                eqs[i][-1] = results[i]
            # Knots internos
            else:
                eqs[i][coefIndex] = knots[i] ** 2
                eqs[i][coefIndex + 1] = knots[i]
                eqs[i][coefIndex + 2] = 1
                # results in the last column
                eqs[i][-1] = results[i]
                # incrementar en uno cada par de elementos
                coefIndex += i % 2 * 3

        # Add n-1 equations
        coefIndex = 0
        for i in range(0, n - 1):
            eqs[2 * n + i][coefIndex] = 2 * X[i + 1]
            eqs[2 * n + i][coefIndex + 1] = 1
            eqs[2 * n + i][coefIndex + 3] = -2 * X[i + 1]
            eqs[2 * n + i][coefIndex + 4] = -1
            # incrementar en 3 el indice de coeficiente
            coefIndex += 3

        # add final equation a1=1
        eqs[3 * n - 1][0] = 1
        return eqs


class GaussJordan:
    """
        Método para intercambiar dos filas de la matriz M
        Entradas: indices de la primer y segunda fila, matriz 
        Salidas: Matriz M modificada
    """
    def intercambiarFilas(self, index1, index2, M):
        for i in range(len(M[index1])):
            temp = M[index1, i]
            M[index1, i] = M[index2, i]
            M[index2, i] = temp
        return M

    """
        Método para buscar un elemento pivote, se implementa pivoteo parcial
        Entradas: filas, columna actual y matriz
        Salidas: Matriz M modificada
    """
    def buscarPivote(self, filas, col, M):
        indiceFila = -1
        for i in range(col + 1, filas):
            if M[i][col] != 0:
                return i
        return indiceFila

    """
        Método para efectuar la eliminación Gaussiana
        Entradas: Número de filas y columnas, matriz
        Salida: Matriz M modificada
    """
    def eliminacionGaussiana(self, f, c, M):
        # Definición de variables
        indicePiv = -1
        M = self.intercambiarFilas(0, -1, M)  # Se sabe que la última fila ya está resuelta
        for i in range(f):
            pivote = M[i][i]
            if pivote == 0:
                indicePiv = self.buscarPivote(f, i, M)  # Se implementa pivoteo parcial
                if indicePiv == -1:
                    print("El sistema no tiene solución", i)
                    exit(0)
                else:
                    M = self.intercambiarFilas(indicePiv, i, M)
                    pivote = M[i][i]

            for j in range(f):  # Realizar la eliminación de los elementos en la columna i
                if (M[j][i] != 0) and (i != j):
                    q = M[j][i] / pivote  # Multiplicador para la eliminación
                    for k in range(f + 1):
                        M[j][k] = M[j][k] - q * M[i][k]

        for i in range(f):
            M[i][-1] /= M[i][i]
            M[i][i] /= M[i][i]


def main():
    # Object to spline class
    spline = SplineInterpolation()
    importObj = ImportData()

    # Importar datos desde un CSV utilizando pandas
    # X, Y = importObj.importDataFrame()
    X = []
    Y = []
    plot = input("Desea graficar las trayectorias? (S/N): ")
    dimX = 1000
    dimY = 1000
    pygame.init()
    screen = pygame.display.set_mode((dimX, dimY))
    pygame.display.set_caption("Presiona enter después de elegir los puntos")
    done = False
    bg = pygame.image.load("map_qro.png")
    screen.blit(bg, (0, 0))
    pygame.display.update()
    clock = pygame.time.Clock()

    while not done:
        for event in pygame.event.get():
            if event.type == KEYDOWN:
                keys = pygame.key.get_pressed()
                if keys[pygame.K_RETURN]:
                    done = True
            elif event.type == MOUSEBUTTONDOWN:
                if pygame.mouse.get_pressed(3)[0]:
                    X.append(pygame.mouse.get_pos()[0])
                    Y.append(pygame.mouse.get_pos()[1])
                    pygame.draw.circle(screen, [50, 50, 255], [X[-1], Y[-1]], 6)
                    pygame.display.update()
                    print(X[-1], Y[-1])
        clock.tick(60)
    X, Y = zip(*sorted(zip(X, Y))) # Se ordenan ambos arreglos juntos
    pygame.display.set_caption("Trayectoria del dron")
    cs = CubicSpline(X, Y, bc_type="clamped")
    eqs = spline.createAugmentedMatrix(X, Y)
    f = (len(X) - 1) * 3
    c = f + 1
    objGJ = GaussJordan()  # objGJ es una variable de la clase "GaussJordan"
    objGJ.eliminacionGaussiana(f, c, eqs)  # Se usa el algoritmo de Gauss-Jordan sobre la matriz

    done = False
    color = (255, 20, 20)
    while not done:
        for i in range(len(X) - 1):
            if done:
                break
            fun = eqs[i * 3, -1] * x ** 2 + eqs[i * 3 + 1, -1] * x + eqs[i * 3 + 2, -1]
            h = (X[i + 1] - X[i]) / 100
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    done = True
                    break

            xi = np.linspace(X[i], X[i + 1], 100, False)
            yi = cs(xi)
            for j in range(100):
                pygame.time.delay(1)
                xc = X[i] + h * j
                yc = fun.subs(x, xc)
                if plot == "N":
                    screen.blit(bg, (0, 0))
                for k in range(len(X)):
                    pygame.draw.circle(screen, [10, 10, 255], [X[k], Y[k]], 6)

                pygame.draw.rect(screen, [30, 255, 30], pygame.Rect(xi[j], yi[j], 8, 8))
                pygame.draw.rect(screen, color, pygame.Rect(xc, yc, 8, 8))
                # Mostrar coordenadas
                print(xi[j], yi[j])
                pygame.display.update()

                keys = pygame.key.get_pressed()
                if keys[pygame.K_ESCAPE]:
                    done = True
        if plot == "S":
            done = True

    done = False
    while not done:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                done = True
                pygame.quit()


if __name__ == "__main__":
    main()
