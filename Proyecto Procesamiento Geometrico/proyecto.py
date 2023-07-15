import openmesh
import argparse
import polyscope as ps
import numpy as np
import math
import os

from normals import calculate_vertex_normals_by_angle, calculate_face_normals


def compute_dihedral_angles(mesh: openmesh.TriMesh):
    dihedral_angles = {}
    face_normals = calculate_face_normals(mesh)
    for arista in mesh.edges():
        if not mesh.is_boundary(arista):
            halfedge = mesh.halfedge_handle(arista, 0)
            f1 = mesh.fh(halfedge)
            n1 = face_normals[f1.idx()]

            op_halfedge = mesh.opposite_halfedge_handle(halfedge)

            f2 = mesh.fh(op_halfedge)
            n2 = face_normals[f2.idx()]

            norm = np.linalg.norm(n1) * np.linalg.norm(n2)
            if norm < 1e-7:
                norm = 1e-7
            dihedral_angles[arista.idx()] = np.abs(np.dot(n1, n2)) / norm
        else:
            dihedral_angles[arista.idx()] = 0
    return dihedral_angles


def intersect_line_triangle(origin: np.ndarray,
                            direction: np.ndarray,
                            point_a: np.ndarray,
                            point_b: np.ndarray,
                            point_c: np.ndarray):
    """
    q1:
    https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
    """
    e_1 = point_b - point_a
    e_2 = point_c - point_a
    n = np.cross(e_1, e_2)
    det = -np.dot(direction, n)
    invdet = 1.0 / det
    a_o = origin - point_a
    dao = np.cross(a_o, direction)
    u = np.dot(e_2, dao) * invdet
    v = -np.dot(e_1, dao) * invdet
    t = np.dot(a_o, n) * invdet
    intersects = det >= 1e-6 and t >= 0.0 and u >= 0.0 and v >= 0.0 and (u + v) <= 1.0

    found_point = origin + t * direction
    return intersects, found_point


def intersect_line_triangle_2(q1, q2, p1, p2, p3):
    def signed_tetra_volume(a, b, c, d):
        return np.sign(np.dot(np.cross(b - a, c - a), d - a) / 6.0)

    n = np.cross(p2 - p1, p3 - p1)
    t = np.dot(p1 - q1, n) / np.dot(q2 - q1, n)
    found_point = q1 + t * (q2 - q1)

    s1 = signed_tetra_volume(q1, p1, p2, p3)
    s2 = signed_tetra_volume(q2, p1, p2, p3)
    intersects = False
    if s1 != s2:
        s3 = signed_tetra_volume(q1, q2, p1, p2)
        s4 = signed_tetra_volume(q1, q2, p2, p3)
        s5 = signed_tetra_volume(q1, q2, p3, p1)
        if s3 == s4 and s4 == s5:

            intersects = True
    return intersects, found_point


def encontrar_arista(mesh: openmesh.TriMesh, v1: openmesh.VertexHandle,
                     v2: openmesh.VertexHandle) -> openmesh.EdgeHandle:
    # https://stackoverflow.com/questions/64243444/openmesh-find-edge-connecting-two-vertices
    he = mesh.find_halfedge(v1, v2)
    if he.is_valid():
        return mesh.edge_handle(he)


def dividir_aristas(mesh: openmesh.TriMesh, Lmax: float):
    edges_left = [e for e in mesh.edges()]
    while len(edges_left) > 0:
        edges_left = list(filter(lambda arista: mesh.calc_edge_length(arista) > Lmax, edges_left))
        edges_iter = edges_left.copy()
        for arista in edges_iter:  # iteramos sobre las aristas de la malla
            halfedge = mesh.halfedge_handle(arista, 0)
            vertice_inicio = mesh.from_vertex_handle(halfedge)
            coords_inicio = mesh.point(vertice_inicio)
            vertice_final = mesh.to_vertex_handle(halfedge)
            coords_final = mesh.point(vertice_final)

            # calcular el punto medio
            vertice_calculado = (coords_inicio + coords_final) / 2.0

            # agregar el vertice
            nuevo_vertice = mesh.add_vertex(vertice_calculado)

            # dividir la arista en el nuevo vertice.
            mesh.split_edge(arista, nuevo_vertice)

            # Encontrar edge nuevo
            arista_nueva = encontrar_arista(mesh, vertice_inicio, nuevo_vertice)
            if arista_nueva == arista:
                # era la otra
                arista_nueva = encontrar_arista(mesh, vertice_final, nuevo_vertice)

            edges_left.append(arista_nueva)


def colapsar_aristas(mesh: openmesh.TriMesh, Lmin: float, Lmax: float):
    def aristas_cortas(aristas):
        return list(filter(lambda arista: mesh.calc_edge_length(arista) < Lmin, aristas))

    n_aristas_cortas = len(aristas_cortas(mesh.edges()))
    ok_count = n_aristas_cortas
    while n_aristas_cortas > 0 and ok_count > 0:
        ok_count = 0
        for arista in mesh.edges():
            longitud_arista = mesh.calc_edge_length(arista)

            if longitud_arista < Lmin:

                halfedge = mesh.halfedge_handle(arista, 0)
                vertice_a = mesh.from_vertex_handle(halfedge)
                vertice_b = mesh.point(mesh.to_vertex_handle(halfedge))

                colapsar = True

                for vecino in mesh.vv(vertice_a):
                    longitud_arista_vecina = math.sqrt(
                        ((vertice_b - mesh.point(vecino)) ** 2).sum())

                    if longitud_arista_vecina > Lmax:
                        colapsar = False
                        break
                if colapsar:
                    if mesh.is_collapse_ok(halfedge):
                        mesh.collapse(halfedge)
                        mesh.garbage_collection()
                        ok_count += 1
                    # break
            print(len(mesh.edges()))
            if len(mesh.edges()) == 3038:
                print("oal")
        n_aristas_cortas = len(aristas_cortas(mesh.edges()))
        print(n_aristas_cortas)


def optimizar_valencia(mesh: openmesh.TriMesh):
    internal_variance = 6
    external_variance = 4

    for arista in mesh.edges():
        if not mesh.is_boundary(arista) and mesh.is_flip_ok(arista):
            halfedge = mesh.halfedge_handle(arista, 0)
            v1 = mesh.from_vertex_handle(halfedge)
            v2 = mesh.to_vertex_handle(halfedge)

            v3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(halfedge))
            # es lo mismo que v1 :(
            # v4 = mesh.to_vertex_handle(mesh.opposite_halfedge_handle(halfedge))
            v4 = mesh.opposite_he_opposite_vh(halfedge)

            deviation_v1 = abs(mesh.valence(v1) - external_variance) if mesh.is_boundary(v1) else abs(
                mesh.valence(v1) - internal_variance)
            deviation_v2 = abs(mesh.valence(v2) - external_variance) if mesh.is_boundary(v2) else abs(
                mesh.valence(v2) - internal_variance)
            deviation_v3 = abs(mesh.valence(v3) - external_variance) if mesh.is_boundary(v3) else abs(
                mesh.valence(v3) - internal_variance)
            deviation_v4 = abs(mesh.valence(v4) - external_variance) if mesh.is_boundary(v4) else abs(
                mesh.valence(v4) - internal_variance)

            deviation_pre = deviation_v1 + deviation_v2 + deviation_v3 + deviation_v4

            mesh.flip(arista)

            deviation_v1 = abs(mesh.valence(v1) - external_variance) if mesh.is_boundary(v1) else abs(
                mesh.valence(v1) - internal_variance)
            deviation_v2 = abs(mesh.valence(v2) - external_variance) if mesh.is_boundary(v2) else abs(
                mesh.valence(v2) - internal_variance)
            deviation_v3 = abs(mesh.valence(v3) - external_variance) if mesh.is_boundary(v3) else abs(
                mesh.valence(v3) - internal_variance)
            deviation_v4 = abs(mesh.valence(v4) - external_variance) if mesh.is_boundary(v4) else abs(
                mesh.valence(v4) - internal_variance)

            deviation_post = deviation_v1 + deviation_v2 + deviation_v3 + deviation_v4

            if deviation_pre <= deviation_post:
                mesh.flip(arista)


def relajacion_tangencial(mesh: openmesh.TriMesh):
    vertex_normals = calculate_vertex_normals_by_angle(mesh)
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
        n = vertex_normals[vertice_id]

        # formula del libro
        p_prime = q[vertice_id] + np.dot(n, (p - q[vertice_id])) * n

        # actualizar posición del vertice
        new_p = proyectar_vertices(mesh,
                                   vertice,
                                   p_prime,
                                   n)
        if new_p is not None:
            mesh.set_point(vertice, new_p)
        # mesh.set_point(vertice, p_prime)


def proyectar_vertices(mesh: openmesh.TriMesh,
                       original_vertex: openmesh.VertexHandle,
                       tangential_pos: np.ndarray,
                       normal: np.ndarray):

    #for neighbor_face in mesh.faces():
    for neighbor_face in mesh.vf(original_vertex):
        a, b, c = [mesh.point(v) for v in mesh.fv(neighbor_face)]
        intersects, new_p = intersect_line_triangle(tangential_pos,
                                                    normal,
                                                    a,
                                                    b,
                                                    c)
        if intersects:
            return new_p
    return None


def remeshing(mesh, edge_length, iteraciones):
    Lmin = 4 / 5 * edge_length
    Lmax = 4 / 3 * edge_length

    # iterar por una cierta cantidad de iteraciones
    for i in range(iteraciones):
        dividir_aristas(mesh, Lmax)
        print("Iteración {}:\tTerminada la división de aristas".format(i))
        colapsar_aristas(mesh, Lmin, Lmax)
        print("Iteración {}:\tTerminado el colapso de aristas".format(i))
        optimizar_valencia(mesh)
        print("Iteración {}:\tTerminada la optimización de valencias".format(i))
        relajacion_tangencial(mesh)
        # proyectar_vertices(mesh)

    return mesh


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, default="", help="Nombre del archivo")
    opt = parser.parse_args()

    if not os.path.exists(opt.file):
        raise FileNotFoundError("Couldn't find file in path {}".format(opt.file))

    mesh = openmesh.read_trimesh(opt.file)

    iteraciones = 10

    lengths = np.array([mesh.calc_edge_length(arista) for arista in mesh.edges()])
    edge_length = lengths.mean()
    remeshing_mesh = remeshing(mesh, edge_length, iteraciones)

    # Mostrar la malla original y la malla con las divisiones
    ps.init()

    # ps.register_surface_mesh("malla_original", mesh.points(), mesh.face_vertex_indices())
    ps.register_surface_mesh("malla_remeshing", remeshing_mesh.points(), remeshing_mesh.face_vertex_indices())

    ps.show()
