#ifndef TETGEN_WRAPPER_H
#define TETGEN_WRAPPER_H

#ifdef _WIN32
#define TETGEN_API __declspec(dllexport)
#else
#define TETGEN_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Tetrahedralize a closed triangle surface mesh using TetGen.
 *
 * Input:
 *   num_points  - number of input vertices
 *   points      - flat array [num_points * 3]: x0,y0,z0, x1,y1,z1, ...
 *   num_facets  - number of input triangles
 *   facets      - flat array [num_facets * 3]: v0,v1,v2, ... (zero-based)
 *   max_volume  - maximum tet volume constraint (0 = no constraint)
 *   min_ratio   - minimum radius-edge quality ratio (0 = no quality, typical 1.2-2.0)
 *
 * Output (allocated by this function, freed by tetgen_free):
 *   out_num_points - number of output vertices (may include Steiner points)
 *   out_points     - flat array [out_num_points * 3]
 *   out_num_tets   - number of output tetrahedra
 *   out_tets       - flat array [out_num_tets * 4] (zero-based node indices)
 *   out_num_faces  - number of output boundary triangles
 *   out_faces      - flat array [out_num_faces * 3] (zero-based node indices)
 *
 * Returns 0 on success, -1 on error (call tetgen_last_error for details).
 */
TETGEN_API int tetgen_mesh(
    int num_points, const double* points,
    int num_facets, const int* facets,
    double max_volume, double min_ratio,
    int quadratic,
    int* out_num_points, double** out_points,
    int* out_num_tets, int** out_tets,
    int* out_nodes_per_tet,
    int* out_num_faces, int** out_faces
);

/* Free a buffer returned by tetgen_mesh. */
TETGEN_API void tetgen_free(void* ptr);

/* Retrieve the last error message. Returns bytes written (excluding null). */
TETGEN_API int tetgen_last_error(char* buf, int buf_len);

#ifdef __cplusplus
}
#endif

#endif /* TETGEN_WRAPPER_H */
