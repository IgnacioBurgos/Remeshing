### Tarea 1: Cálculo de Curvatura

import openmesh
import numpy as np
import polyscope as ps


def calculate_face_normals(mesh: openmesh.TriMesh) -> np.ndarray:
    face_normals = np.zeros((mesh.n_faces(), 3))
    for face in mesh.faces():
        v0, v1, v2 = [mesh.points()[v.idx()] for v in mesh.fv(face)]
        normal = np.cross(v1 - v0, v2 - v0)
        normal_norm = np.linalg.norm(normal)
        if normal_norm < 1e-7:
            normal_norm = 1e-7
        normal /= normal_norm
        face_normals[face.idx(), :] = normal

    return face_normals


def angle(center: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> float:
    vec_a = p1 - center
    vec_b = p2 - center
    norm_a = max(np.linalg.norm(vec_a), 1e-8)
    norm_b = max(np.linalg.norm(vec_b), 1e-8)
    theta = np.dot(vec_b, vec_a) / (norm_a * norm_b)
    return np.arccos(theta)


def angle_sine(center_v: np.ndarray, p_r: np.ndarray, p_t: np.ndarray) -> float:
    e_vt = p_t - center_v
    e_vr = center_v - p_r
    norm_e_vt = max(np.linalg.norm(e_vt), 1e-8)
    norm_e_vr = max(np.linalg.norm(e_vr), 1e-8)
    sin_theta = np.min((1.0, np.linalg.norm(np.cross(e_vt, e_vr)) / (norm_e_vt * norm_e_vr)))
    return np.arcsin(sin_theta)


def calculate_vertex_face_angles(mesh: openmesh.TriMesh, method="cos") -> np.ndarray:
    empty_dicts = np.frompyfunc(lambda x: {}, 1, 1)
    vertex_angles = empty_dicts(np.zeros(mesh.n_vertices(), dtype=dict))
    for face in mesh.faces():
        p_vertices = [v for v in mesh.fv(face)]
        vertices = [mesh.points()[v.idx()] for v in p_vertices]
        for i in range(3):
            i_prev = (i - 1) % 3
            i_next = (i + 1) % 3
            if method == "cos":
                vert_angle = angle(vertices[i], vertices[i_prev], vertices[i_next])
            elif method == "sin":
                vert_angle = angle_sine(vertices[i], vertices[i_prev], vertices[i_next])
            vertex_angles[p_vertices[i].idx()][face.idx()] = vert_angle
    return vertex_angles


def calculate_vertex_normals_by_angle(mesh: openmesh.TriMesh, method: str = "cos") -> np.ndarray:
    face_normals = calculate_face_normals(mesh)
    vertex_face_angles = calculate_vertex_face_angles(mesh, method)

    vertex_normals = np.zeros((mesh.n_vertices(), 3))

    for vertex in mesh.vertices():
        i_v = vertex.idx()
        for face in mesh.vf(vertex):
            vertex_normals[i_v, :] += face_normals[face.idx(), :] * vertex_face_angles[i_v][face.idx()]
        normals_norm = np.max((np.linalg.norm(vertex_normals[i_v, :]), 1e-8))
        vertex_normals[i_v, :] /= normals_norm

    return vertex_normals
