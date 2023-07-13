import openmesh
import argparse
import polyscope as ps
import numpy as np
import math


def dividir_aristas(mesh, Lmax):

    
    while True:

        edge_encontrada = False

        for arista in mesh.edges(): #iteramos sobre las aristas de la malla
            longitud_arista = mesh.calc_edge_length(arista)

            if longitud_arista > Lmax: 
                halfedge = mesh.halfedge_handle(arista, 0)  
                vertice_inicio = mesh.point(mesh.from_vertex_handle(halfedge)) 
                vertice_final = mesh.point(mesh.to_vertex_handle(halfedge))

                # calcular el punto medio
                vertice_calculado = (vertice_inicio + vertice_final) / 2.0

                # agregar el vertice
                nuevo_vertice = mesh.add_vertex(vertice_calculado)
            
                # dividir la arista en el nuevo vertice.
                mesh.split_edge(arista, nuevo_vertice)

        if not edge_encontrada:
            break

def colapsar_aristas(mesh, Lmin, Lmax):
    
    while True:
        edge_encontrada = False

        for arista in mesh.edges():
            longitud_arista = mesh.calc_edge_length(arista)

            if longitud_arista < Lmin:
                
                halfedge = mesh.halfedge_handle(arista, 0)
                vertice_a = mesh.from_vertex_handle(halfedge)
                vertice_b = mesh.point(mesh.to_vertex_handle(halfedge))

                colapsar = True

                for vecino in mesh.vv(vertice_a):
                    longitud_arista_vecina = math.sqrt((vertice_b[0] - mesh.point(vecino)[0]) ** 2 + (vertice_b[1] - mesh.point(vecino)[1]) ** 2 + (vertice_b[2] - mesh.point(vecino)[2]) ** 2)

                    if longitud_arista_vecina > Lmax:
                        colapsar = False
                        break

                if colapsar:
                    mesh.collapse(halfedge)
                    edge_encontrada = True
                    break

            
        if not edge_encontrada:
            break


    
def optimizar_valencia(mesh):
    internal_variance = 6
    external_variance = 4
    
    for arista in mesh.edges():
        if not mesh.is_boundary(arista):
            halfedge = mesh.halfedge_handle(arista, 0)
            v1 = mesh.from_vertex_handle(halfedge)
            v2 = mesh.to_vertex_handle(halfedge)

            v3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(halfedge))
            v4 = mesh.to_vertex_handle(mesh.opposite_halfedge_handle(halfedge))

            deviation_v1 = abs(mesh.valence(v1) - external_variance) if mesh.is_boundary(v1) else abs(mesh.valence(v1) - internal_variance)
            deviation_v2 = abs(mesh.valence(v2) - external_variance) if mesh.is_boundary(v2) else abs(mesh.valence(v2) - internal_variance)
            deviation_v3 = abs(mesh.valence(v3) - external_variance) if mesh.is_boundary(v3) else abs(mesh.valence(v3) - internal_variance)
            deviation_v4 = abs(mesh.valence(v4) - external_variance) if mesh.is_boundary(v4) else abs(mesh.valence(v4) - internal_variance)
            
            deviation_pre = deviation_v1 + deviation_v2 + deviation_v3 + deviation_v4

            mesh.flip(arista)

            deviation_v1 = abs(mesh.valence(v1) - external_variance) if mesh.is_boundary(v1) else abs(mesh.valence(v1) - internal_variance)
            deviation_v2 = abs(mesh.valence(v2) - external_variance) if mesh.is_boundary(v2) else abs(mesh.valence(v2) - internal_variance)
            deviation_v3 = abs(mesh.valence(v3) - external_variance) if mesh.is_boundary(v3) else abs(mesh.valence(v3) - internal_variance)
            deviation_v4 = abs(mesh.valence(v4) - external_variance) if mesh.is_boundary(v4) else abs(mesh.valence(v4) - internal_variance)

            deviation_post = deviation_v1 + deviation_v2 + deviation_v3 + deviation_v4

            if deviation_pre <= deviation_post:
                mesh.flip(arista)
            


def relajacion_tangencial(mesh):

    q = {}
    
    for vertice in mesh.vertices():
        vertice_id = vertice.idx()

        puntos_vecinos = []

        for vertice_vecino in mesh.vv(vertice):
            puntos_vecinos.append(mesh.point(vertice_vecino))
        
        baricentro = np.mean(puntos_vecinos, axis=0)

        q[vertice_id] = baricentro 

    for vertice in mesh.vertices():
        vertice_id = vertice.idx()

        p = mesh.point(vertice)
        n = mesh.normal(vertice)

        # formula del libro
        p = q[vertice_id] + np.dot(n,(p - q[vertice_id])) * n

        #actualizar posición del vertice
        mesh.set_point(vertice, p)


def proyectar_vertices(mesh):
    print('falta codigo aqui')

def remeshing(mesh, edge_length, iteraciones):

    Lmin = 4/5 * edge_length
    Lmax = 4/3 * edge_length
    
    #iterar por una cierta cantidad de iteraciones
    for i in range(iteraciones):
        print(i)
        dividir_aristas(mesh, Lmax)
        colapsar_aristas(mesh, Lmin, Lmax)
        optimizar_valencia(mesh)
        relajacion_tangencial(mesh)
        #proyectar_vertices(mesh)

    return mesh

parser = argparse.ArgumentParser()
parser.add_argument("--file", type=str, default="", help="Nombre del archivo")
opt = parser.parse_args()

mesh = openmesh.read_trimesh(opt.file)

iteraciones = 10
edge_length = 0.01
remeshing_mesh = remeshing(mesh, edge_length, iteraciones)

# Mostrar la malla original y la malla con las divisiones
ps.init()

#ps.register_surface_mesh("malla_original", mesh.points(), mesh.face_vertex_indices())
ps.register_surface_mesh("malla_remeshing", remeshing_mesh.points(), remeshing_mesh.face_vertex_indices())

ps.show()
