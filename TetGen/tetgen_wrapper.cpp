#ifndef TETLIBRARY
#define TETLIBRARY
#endif
#include "tetgen.h"
#include "tetgen_wrapper.h"

#include <cstring>
#include <cstdlib>
#include <string>

static thread_local std::string g_last_error;

static void set_error(const char* msg) {
    g_last_error = msg ? msg : "";
}

extern "C" {

TETGEN_API int tetgen_mesh(
    int num_points, const double* points,
    int num_facets, const int* facets,
    double max_volume, double min_ratio,
    int* out_num_points, double** out_points,
    int* out_num_tets, int** out_tets,
    int* out_num_faces, int** out_faces)
{
    if (!points || num_points < 4 || !facets || num_facets < 4) {
        set_error("Invalid input: need at least 4 points and 4 facets");
        return -1;
    }

    *out_num_points = 0; *out_points = nullptr;
    *out_num_tets = 0;   *out_tets = nullptr;
    *out_num_faces = 0;  *out_faces = nullptr;

    try {
        tetgenio in, out;
        in.firstnumber = 0;
        in.numberofpoints = num_points;
        in.pointlist = new REAL[num_points * 3];
        std::memcpy(in.pointlist, points, num_points * 3 * sizeof(double));

        in.numberoffacets = num_facets;
        in.facetlist = new tetgenio::facet[num_facets];
        in.facetmarkerlist = new int[num_facets];

        for (int i = 0; i < num_facets; i++) {
            tetgenio::facet* f = &in.facetlist[i];
            f->numberofpolygons = 1;
            f->polygonlist = new tetgenio::polygon[1];
            f->numberofholes = 0;
            f->holelist = nullptr;

            tetgenio::polygon* p = &f->polygonlist[0];
            p->numberofvertices = 3;
            p->vertexlist = new int[3];
            p->vertexlist[0] = facets[i * 3 + 0];
            p->vertexlist[1] = facets[i * 3 + 1];
            p->vertexlist[2] = facets[i * 3 + 2];

            in.facetmarkerlist[i] = 1;
        }

        // Build switch string: p=PLC, z=zero-based, Q=quiet
        std::string switches = "pzQ";
        if (min_ratio > 0.0) {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "q%.6g", min_ratio);
            switches += buf;
        }
        if (max_volume > 0.0) {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "a%.6g", max_volume);
            switches += buf;
        }

        tetrahedralize(const_cast<char*>(switches.c_str()), &in, &out);

        if (out.numberoftetrahedra == 0) {
            set_error("TetGen produced no tetrahedra");
            return -1;
        }

        // Copy output points
        int np = out.numberofpoints;
        *out_num_points = np;
        *out_points = (double*)std::malloc(np * 3 * sizeof(double));
        std::memcpy(*out_points, out.pointlist, np * 3 * sizeof(double));

        // Copy output tets
        int nt = out.numberoftetrahedra;
        int nc = out.numberofcorners; // should be 4
        *out_num_tets = nt;
        *out_tets = (int*)std::malloc(nt * 4 * sizeof(int));
        if (nc == 4) {
            std::memcpy(*out_tets, out.tetrahedronlist, nt * 4 * sizeof(int));
        } else {
            // Higher-order: only copy first 4 corners per tet
            for (int i = 0; i < nt; i++) {
                (*out_tets)[i * 4 + 0] = out.tetrahedronlist[i * nc + 0];
                (*out_tets)[i * 4 + 1] = out.tetrahedronlist[i * nc + 1];
                (*out_tets)[i * 4 + 2] = out.tetrahedronlist[i * nc + 2];
                (*out_tets)[i * 4 + 3] = out.tetrahedronlist[i * nc + 3];
            }
        }

        // Copy output boundary faces
        int nf = out.numberoftrifaces;
        *out_num_faces = nf;
        *out_faces = (int*)std::malloc(nf * 3 * sizeof(int));
        std::memcpy(*out_faces, out.trifacelist, nf * 3 * sizeof(int));

        set_error("");
        return 0;

    } catch (int errcode) {
        // TetGen signals errors by throwing an int
        char buf[128];
        std::snprintf(buf, sizeof(buf), "TetGen error code %d", errcode);
        set_error(buf);
        return -1;
    } catch (const std::exception& ex) {
        set_error(ex.what());
        return -1;
    } catch (...) {
        set_error("Unknown error in TetGen");
        return -1;
    }
}

TETGEN_API void tetgen_free(void* ptr) {
    std::free(ptr);
}

TETGEN_API int tetgen_last_error(char* buf, int buf_len) {
    if (!buf || buf_len <= 0) return 0;
    int n = (int)g_last_error.size();
    if (n == 0) { buf[0] = '\0'; return 0; }
    int copy = n < buf_len ? n : buf_len - 1;
    std::memcpy(buf, g_last_error.c_str(), copy);
    buf[copy] = '\0';
    return copy;
}

} // extern "C"
