#ifndef TRI2D_HPP
#define TRI2D_HPP

#include <ap_fixed.h>
#include <hls_stream.h>
#include <ap_int.h>
#include <ap_axi_sdata.h>

// ============================================
// Constant Definitions
// ============================================
#define EPSILON 0.0001

// Capacity constants
#define MAX_NO_POINTS 16003
#define MAX_NO_TRIANGLES 32000
#define MAX_FACES 32000
#define MAX_NO_HALFEDGES (3 * MAX_NO_TRIANGLES + 3)
#define MAX_POINT_PER_CELL 160
#define MAX_GRID_SIZE 32
#define HILBERT_ORDER 8

// Bit widths
#define MNPB 14  // Point index bit width

// ============================================
// Type Definitions
// ============================================
typedef ap_axis<64, 16, 5, 6> axi_i_t;
typedef ap_axis<64, 16, 5, 6> axi_o_t;
typedef ap_uint<1>  dout_t;
typedef ap_uint<14> dnode_t;
typedef ap_uint<15> didx_t;
typedef ap_uint<17> didxh_t;
typedef ap_uint<32> hilbert_t;

typedef ap_fixed<24, 4>  fixed_t;
typedef ap_fixed<33, 6>  fixed_calc_t;
typedef ap_fixed<48, 10> fixed_det_t;

const dnode_t INVALID_NODE = (1 << 13) - 1;
const didx_t INVALID_INDEX = (1 << 15) - 1;
const didxh_t INVALID_HALFEDGE = (1 << 16) - 1;

// ============================================
// Struct Definitions
// ============================================
struct FixedPoint {
    fixed_t x, y;
    FixedPoint() : x(0), y(0) {}
};

// ============================================
// HalfEdge SoA
// Edge layout:
//   face 0              -> edges 0, 1, 2
//   boundary ghost cycle -> edges 3, 4, 5
//   face f (f >= 1)     -> edges 6 + 3*(f-1) + {0,1,2}
//
// he_face[]: BRAM lookup replacing e/3 division
// ============================================
extern dnode_t  he_tail[MAX_NO_HALFEDGES];
extern didxh_t  he_twin[MAX_NO_HALFEDGES];
extern didx_t   he_face[MAX_NO_HALFEDGES];   // [OPT] face-of-edge, BRAM
extern dout_t   he_used[MAX_NO_HALFEDGES];

// ============================================
// Face SoA  (was: FaceTopo struct)
// ============================================
// face_used : 1-bit alive flag (hot, BRAM)
// face_next : freelist link (cold, LUTRAM)
// face_v1/v2/v3 : triangle vertex indices (BRAM)
// face_coord    : triangle vertex coordinates (URAM)
struct FaceCoord {
    fixed_t x1, y1;
    fixed_t x2, y2;
    fixed_t x3, y3;
};

extern ap_uint<1> face_used[MAX_NO_TRIANGLES];
extern didx_t     face_next[MAX_NO_TRIANGLES];
extern dnode_t    face_v1[MAX_NO_TRIANGLES];
extern dnode_t    face_v2[MAX_NO_TRIANGLES];
extern dnode_t    face_v3[MAX_NO_TRIANGLES];
extern FaceCoord  face_coord[MAX_NO_TRIANGLES];

// ============================================
// Vertex SoA
// Removed (zero/write-only): vtx_hilbert, vtx_edge,
//                            vtx_prev, vtx_next, vtx_processed
// ============================================
extern fixed_t     vtx_x[MAX_NO_POINTS];
extern fixed_t     vtx_y[MAX_NO_POINTS];
extern ap_uint<1>  vtx_used[MAX_NO_POINTS];

// ============================================
// Grid structure
// ============================================
struct EnhancedAdaptiveGrid {
    fixed_t x_min, y_min;
    fixed_t x_max, y_max;
    fixed_t cell_width, cell_height;
    int grid_size_x, grid_size_y;
    int total_cells;

    EnhancedAdaptiveGrid() {
#pragma HLS INLINE
        x_min = y_min = fixed_t(1e30f);
        x_max = y_max = fixed_t(-1e30f);
        cell_width = cell_height = 0;
        grid_size_x = grid_size_y = 0;
        total_cells = 0;
    }
};

// ============================================
// Triangle-centric edge helpers
// ============================================

// face -> base half-edge index
//   f == 0 : 0
//   f >= 1 : 6 + 3*(f-1) = 3*f + 3
inline didxh_t face_base_edge(didx_t f) {
#pragma HLS INLINE
    return (f == 0) ? (didxh_t)0 : (didxh_t)(6 + 3 * (f - 1));
}

// half-edge -> owning face   (1 BRAM read, replaces e/3 divider)
inline didx_t edge_to_face(didxh_t e) {
#pragma HLS INLINE
    return he_face[e];
}

// -------------------------------------------------------
// next / prev within the same triangle
// -------------------------------------------------------
inline didxh_t next_edge_of(didxh_t e) {
#pragma HLS INLINE
    if (e < 3) return (e == 2) ? (didxh_t)0 : (didxh_t)(e + 1);
    if (e < 6) return (e == 5) ? (didxh_t)3 : (didxh_t)(e + 1);
    didxh_t b = face_base_edge(he_face[e]);
    ap_uint<2> k = e - b;
    return (k == 2) ? b : (didxh_t)(e + 1);
}

inline didxh_t prev_edge_of(didxh_t e) {
#pragma HLS INLINE
    if (e < 3) return (e == 0) ? (didxh_t)2 : (didxh_t)(e - 1);
    if (e < 6) return (e == 3) ? (didxh_t)5 : (didxh_t)(e - 1);
    didxh_t b = face_base_edge(he_face[e]);
    ap_uint<2> k = e - b;
    return (k == 0) ? (didxh_t)(b + 2) : (didxh_t)(e - 1);
}

// ============================================
// Global Variable Declarations
// ============================================
extern EnhancedAdaptiveGrid enhanced_grid;
extern didx_t tri_num;
extern int num_point;

// ============================================
// Function Declarations
// ============================================
fixed_t absf(fixed_t val);
bool inCircle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c);
bool isPointInTriangle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c);

void getTriangleVertices(didx_t faceIdx, dnode_t& v1, dnode_t& v2, dnode_t& v3);
void updateFaceCache(didx_t faceIdx);

hilbert_t hilbert_xy2d(int n, int x, int y);
hilbert_t getHilbertCode(fixed_t x, fixed_t y, fixed_t x_min, fixed_t y_min,
    fixed_t x_max, fixed_t y_max, int order);

void initEnhancedGrid(EnhancedAdaptiveGrid& grid, int num_vertices,
                      int* grid_counts, int* grid_triangles);
int getGridIndex(const EnhancedAdaptiveGrid& grid, fixed_t x, fixed_t y);

didx_t allocateFace();
void freeFace(didx_t faceIdx);

int locateTriangleSimple(const FixedPoint& p, int triangleCount,
                         int* grid_counts, int* grid_triangles);
void connectTwinEdges(didxh_t e1, didxh_t e2);

void legalizeEdge(didx_t faceIdx, didxh_t edgeIdx,
                  int* grid_counts, int* grid_triangles,
                  didx_t touched_faces[], int& touched_count);

void insertSiteEnhanced(dnode_t p, didx_t containingFace,
                        int* grid_counts, int* grid_triangles);

FixedPoint calculateTriangleCentroid(didx_t triangleIdx);

void assignTriangleToGrid(didx_t triangleIdx,
                          int* grid_counts, int* grid_triangles);

void updateFaceCacheAndAssignGrid(didx_t triangleIdx,
                                  int* grid_counts, int* grid_triangles);

inline FixedPoint getFixedPoint(dnode_t pointIdx) {
#pragma HLS INLINE
    FixedPoint p;
    p.x = vtx_x[pointIdx];
    p.y = vtx_y[pointIdx];
    return p;
}

void tri2d(hls::stream<axi_i_t>& x_in_stream,
           hls::stream<axi_o_t>& cl_out_stream,
           int* grid_counts,
           int* grid_triangles,
           volatile int& debug_count);

#endif // TRI2D_HPP
