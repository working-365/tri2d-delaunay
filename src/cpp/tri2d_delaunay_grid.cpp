/*
 * KV260 Optimized DCEL-based Delaunay Triangulation
 * ====== DDR Grid SoA Version ======
 *
 * CHANGES (this revision):
 *   1. Deleted dead arrays: vtx_hilbert, vtx_edge, vtx_prev,
 *      vtx_next, vtx_processed (zero or write-only access).
 *   2. Removed redundant `v1>=0 && v2>=0 && v3>=0` check in
 *      locateTriangle fallback (face_used==1 already implies valid).
 *   3. Split FaceTopo struct into SoA arrays:
 *        face_used  (BRAM, 1-bit, hot)
 *        face_next  (LUTRAM, freelist, cold)
 *        face_v1/v2/v3 (BRAM)
 *      face_coord stays in URAM.
 */

#include "tri2d.hpp"
#include <iostream>
#include <cmath>
using namespace std;

// ============================================
// Global Variable Definitions
// ============================================

didx_t tri_num = 0;
int num_point = 16000;
int num_vertice = 16003;
static didx_t freeFaceHead = INVALID_INDEX;
static const int COMPACTION_THRESHOLD = 200;
static int processed_actual_points = 0;

// HalfEdge SoA
dnode_t  he_tail[MAX_NO_HALFEDGES];
didxh_t  he_twin[MAX_NO_HALFEDGES];
dout_t   he_used[MAX_NO_HALFEDGES];
didx_t   he_face[MAX_NO_HALFEDGES];

// Face SoA (was: FaceTopo struct)
ap_uint<1> face_used[MAX_NO_TRIANGLES];
didx_t     face_next[MAX_NO_TRIANGLES];
dnode_t    face_v1[MAX_NO_TRIANGLES];
dnode_t    face_v2[MAX_NO_TRIANGLES];
dnode_t    face_v3[MAX_NO_TRIANGLES];
FaceCoord  face_coord[MAX_NO_TRIANGLES];

// Vertex SoA  (only x/y/used kept; rest deleted)
fixed_t     vtx_x[MAX_NO_POINTS];
fixed_t     vtx_y[MAX_NO_POINTS];
ap_uint<1>  vtx_used[MAX_NO_POINTS];

// Grid metadata
EnhancedAdaptiveGrid enhanced_grid;

// --------------------------------------------------
// DDR Grid local BRAM cache (5x5 = 25 cells)
// --------------------------------------------------
static int local_counts[25];
static int local_tris[25 * MAX_POINT_PER_CELL];
static int cached_center_x = -999;
static int cached_center_y = -999;

// ============================================
// Helper Functions
// ============================================

fixed_t absf(fixed_t val) {
#pragma HLS INLINE
    return (val < 0) ? (fixed_t)(-val) : val;
}

hilbert_t getHilbertCode(fixed_t x, fixed_t y, fixed_t x_min, fixed_t y_min,
    fixed_t x_max, fixed_t y_max, int order) {
    float x_f = x.to_float();
    float y_f = y.to_float();
    float xmin_f = x_min.to_float();
    float ymin_f = y_min.to_float();
    float xmax_f = x_max.to_float();
    float ymax_f = y_max.to_float();

    float norm_x = (x_f - xmin_f) / (xmax_f - xmin_f);
    float norm_y = (y_f - ymin_f) / (ymax_f - ymin_f);
    norm_x = std::max(0.0f, std::min(1.0f, norm_x));
    norm_y = std::max(0.0f, std::min(1.0f, norm_y));

    int grid_size = 1 << order;
    int hx = (int)(norm_x * (grid_size - 1));
    int hy = (int)(norm_y * (grid_size - 1));

    return hilbert_xy2d(order, hx, hy);
}

hilbert_t hilbert_xy2d(int n, int x, int y) {
#pragma HLS INLINE
    hilbert_t d = 0;
    for (int s = 1; s < (1 << n); s *= 2) {
#pragma HLS PIPELINE
        int rx = (x & s) > 0;
        int ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        if (ry == 0) {
            if (rx == 1) { x = s - 1 - x; y = s - 1 - y; }
            int t = x; x = y; y = t;
        }
    }
    return d;
}

static inline void getTrianglePoints(didx_t tri,
    FixedPoint& c, FixedPoint& m12, FixedPoint& m23, FixedPoint& m31) {
#pragma HLS INLINE
    FaceCoord& fc = face_coord[tri];
    c.x = (fc.x1 + fc.x2 + fc.x3) / fixed_t(3.0);
    c.y = (fc.y1 + fc.y2 + fc.y3) / fixed_t(3.0);
    m12.x = (fc.x1 + fc.x2) / fixed_t(2.0);
    m12.y = (fc.y1 + fc.y2) / fixed_t(2.0);
    m23.x = (fc.x2 + fc.x3) / fixed_t(2.0);
    m23.y = (fc.y2 + fc.y3) / fixed_t(2.0);
    m31.x = (fc.x3 + fc.x1) / fixed_t(2.0);
    m31.y = (fc.y3 + fc.y1) / fixed_t(2.0);
}

// ============================================
// updateFaceCache: base+1/+2 direct
// ============================================
inline void updateFaceCache(didx_t faceIdx) {
#pragma HLS INLINE
    didxh_t e1 = face_base_edge(faceIdx);
    didxh_t e2 = (didxh_t)(e1 + 1);
    didxh_t e3 = (didxh_t)(e1 + 2);

    dnode_t v1 = he_tail[e1];
    dnode_t v2 = he_tail[e2];
    dnode_t v3 = he_tail[e3];

    face_v1[faceIdx] = v1;
    face_v2[faceIdx] = v2;
    face_v3[faceIdx] = v3;

    face_coord[faceIdx].x1 = vtx_x[v1]; face_coord[faceIdx].y1 = vtx_y[v1];
    face_coord[faceIdx].x2 = vtx_x[v2]; face_coord[faceIdx].y2 = vtx_y[v2];
    face_coord[faceIdx].x3 = vtx_x[v3]; face_coord[faceIdx].y3 = vtx_y[v3];
}

// ============================================
// getTriangleVertices: base+1/+2 direct
// ============================================
void getTriangleVertices(didx_t faceIdx, dnode_t& v1, dnode_t& v2, dnode_t& v3) {
#pragma HLS INLINE
    didxh_t e1 = face_base_edge(faceIdx);
    didxh_t e2 = (didxh_t)(e1 + 1);
    didxh_t e3 = (didxh_t)(e1 + 2);

    v1 = he_tail[e1];
    v2 = he_tail[e2];
    v3 = he_tail[e3];
}

bool isPointInTriangle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c) {
    fixed_calc_t cross1 = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
    fixed_calc_t cross2 = (c.x - b.x) * (p.y - b.y) - (c.y - b.y) * (p.x - b.x);
    fixed_calc_t cross3 = (a.x - c.x) * (p.y - c.y) - (a.y - c.y) * (p.x - c.x);
    bool pos = (cross1 >= 0) && (cross2 >= 0) && (cross3 >= 0);
    bool neg = (cross1 <= 0) && (cross2 <= 0) && (cross3 <= 0);
    return pos || neg;
}

inline int getGridIndex(const EnhancedAdaptiveGrid& grid, fixed_t x, fixed_t y) {
    int grid_x = max(0, min(grid.grid_size_x - 1,
        int((x - grid.x_min) / grid.cell_width)));
    int grid_y = max(0, min(grid.grid_size_y - 1,
        int((y - grid.y_min) / grid.cell_height)));
    return grid_y * grid.grid_size_x + grid_x;
}

// ============================================
// initEnhancedGrid
// ============================================
void initEnhancedGrid(EnhancedAdaptiveGrid& grid, int num_vertice,
                      int* grid_counts, int* grid_triangles) {
#pragma HLS INLINE off

    grid.x_min = fixed_t(1e30f);
    grid.y_min = fixed_t(1e30f);
    grid.x_max = fixed_t(-1e30f);
    grid.y_max = fixed_t(-1e30f);

    for (int i = 3; i < num_vertice; i++) {
        if (vtx_used[i]) {
            fixed_t x = vtx_x[i];
            fixed_t y = vtx_y[i];
            grid.x_min = (x < grid.x_min) ? x : grid.x_min;
            grid.y_min = (y < grid.y_min) ? y : grid.y_min;
            grid.x_max = (x > grid.x_max) ? x : grid.x_max;
            grid.y_max = (y > grid.y_max) ? y : grid.y_max;
        }
    }

    fixed_calc_t margin_x = (grid.x_max - grid.x_min) * fixed_calc_t(0.2);
    fixed_calc_t margin_y = (grid.y_max - grid.y_min) * fixed_calc_t(0.2);
    grid.x_min -= fixed_t(margin_x);
    grid.y_min -= fixed_t(margin_y);
    grid.x_max += fixed_t(margin_x);
    grid.y_max += fixed_t(margin_y);

    grid.grid_size_x = MAX_GRID_SIZE;
    grid.grid_size_y = MAX_GRID_SIZE;
    grid.total_cells = grid.grid_size_x * grid.grid_size_y;

    float x_range = (grid.x_max - grid.x_min).to_float();
    float y_range = (grid.y_max - grid.y_min).to_float();
    grid.cell_width = fixed_t(x_range / grid.grid_size_x);
    grid.cell_height = fixed_t(y_range / grid.grid_size_y);

INIT_COUNTS: for (int i = 0; i < grid.total_cells; i++) {
#pragma HLS PIPELINE II=1
        grid_counts[i] = 0;
    }

    int total_slots = grid.total_cells * MAX_POINT_PER_CELL;
INIT_TRIS: for (int i = 0; i < total_slots; i++) {
#pragma HLS PIPELINE II=1
        grid_triangles[i] = (int)INVALID_INDEX;
    }
}

FixedPoint calculateTriangleCentroid(didx_t triangleIdx) {
#pragma HLS INLINE
    dnode_t v1, v2, v3;
    getTriangleVertices(triangleIdx, v1, v2, v3);
    FixedPoint p1 = getFixedPoint(v1);
    FixedPoint p2 = getFixedPoint(v2);
    FixedPoint p3 = getFixedPoint(v3);
    FixedPoint centroid;
    centroid.x = (p1.x + p2.x + p3.x) / fixed_t(3.0);
    centroid.y = (p1.y + p2.y + p3.y) / fixed_t(3.0);
    return centroid;
}

// ============================================
// tryInsertTriToCell_DDR
// ============================================
static inline void tryInsertTriToCell_DDR(int* grid_counts, int* grid_triangles,
                                          int grid_idx, didx_t triangleIdx,
                                          int total_cells) {
#pragma HLS INLINE
    if (grid_idx < 0 || grid_idx >= total_cells) return;

    int count = grid_counts[grid_idx];

    int base = grid_idx * MAX_POINT_PER_CELL;
    if (count < 0) return;

    int write_idx = 0;
    for (int i = 0; i < MAX_POINT_PER_CELL; i++) {
#pragma HLS PIPELINE II=1
        if (i >= count) break;
        int existing = grid_triangles[base + i];
        if (existing == (int)triangleIdx) return;
        if (existing != (int)INVALID_INDEX &&
            existing < MAX_NO_TRIANGLES &&
            face_used[(didx_t)existing]) {
            grid_triangles[base + write_idx] = existing;
            write_idx++;
        }
    }

    grid_counts[grid_idx] = write_idx;

    if (write_idx < MAX_POINT_PER_CELL) {
        grid_triangles[base + write_idx] = (int)triangleIdx;
        grid_counts[grid_idx] = write_idx + 1;
    }
}

// ============================================
// assignTriangleToGrid
// ============================================
void assignTriangleToGrid(didx_t triangleIdx,
                          int* grid_counts, int* grid_triangles) {
#pragma HLS INLINE off

    if (triangleIdx >= MAX_NO_TRIANGLES) return;

    if (!face_used[triangleIdx]) return;

    // assumes face_coord[triangleIdx] already updated by updateFaceCache.
    FixedPoint c, m12, m23, m31;
    getTrianglePoints(triangleIdx, c, m12, m23, m31);

    int ids[4];
#pragma HLS ARRAY_PARTITION variable=ids complete
    ids[0] = getGridIndex(enhanced_grid, c.x,   c.y);
    ids[1] = getGridIndex(enhanced_grid, m12.x, m12.y);
    ids[2] = getGridIndex(enhanced_grid, m23.x, m23.y);
    ids[3] = getGridIndex(enhanced_grid, m31.x, m31.y);

    int uniq[4];
#pragma HLS ARRAY_PARTITION variable=uniq complete
    int n = 0;

UNIQ_IDS:
    for (int k = 0; k < 4; k++) {
#pragma HLS UNROLL
        int cid = ids[k];
        if (cid < 0 || cid >= enhanced_grid.total_cells) continue;

        bool dup = false;
    CHECK_DUP:
        for (int j = 0; j < 4; j++) {
#pragma HLS UNROLL
            if (j >= n) break;
            if (uniq[j] == cid) {
                dup = true;
                break;
            }
        }

        if (!dup) {
            uniq[n] = cid;
            n++;
        }
    }

INSERT_UNIQ:
    for (int i = 0; i < 4; i++) {
#pragma HLS UNROLL
        if (i >= n) break;
        tryInsertTriToCell_DDR(grid_counts, grid_triangles,
                               uniq[i], triangleIdx,
                               enhanced_grid.total_cells);
    }
}

void updateFaceCacheAndAssignGrid(didx_t triangleIdx,
                                  int* grid_counts, int* grid_triangles) {
#pragma HLS INLINE off

    if (triangleIdx >= MAX_NO_TRIANGLES) return;

    if (!face_used[triangleIdx]) return;

    // 1) take 3 vertices straight from half-edge
    didxh_t e0 = face_base_edge(triangleIdx);
    didxh_t e1 = (didxh_t)(e0 + 1);
    didxh_t e2 = (didxh_t)(e0 + 2);

    dnode_t v1 = he_tail[e0];
    dnode_t v2 = he_tail[e1];
    dnode_t v3 = he_tail[e2];
    fixed_t x1 = vtx_x[v1], y1 = vtx_y[v1];
    fixed_t x2 = vtx_x[v2], y2 = vtx_y[v2];
    fixed_t x3 = vtx_x[v3], y3 = vtx_y[v3];
    face_v1[triangleIdx] = v1;
    face_v2[triangleIdx] = v2;
    face_v3[triangleIdx] = v3;

    face_coord[triangleIdx].x1 = x1; face_coord[triangleIdx].y1 = y1;
    face_coord[triangleIdx].x2 = x2; face_coord[triangleIdx].y2 = y2;
    face_coord[triangleIdx].x3 = x3; face_coord[triangleIdx].y3 = y3;
    FixedPoint c, m12, m23, m31;
    c.x   = (x1 + x2 + x3) / fixed_t(3.0);
    c.y   = (y1 + y2 + y3) / fixed_t(3.0);
    m12.x = (x1 + x2) / fixed_t(2.0);
    m12.y = (y1 + y2) / fixed_t(2.0);
    m23.x = (x2 + x3) / fixed_t(2.0);
    m23.y = (y2 + y3) / fixed_t(2.0);
    m31.x = (x3 + x1) / fixed_t(2.0);
    m31.y = (y3 + y1) / fixed_t(2.0);

    int ids[4];
#pragma HLS ARRAY_PARTITION variable=ids complete
    ids[0] = getGridIndex(enhanced_grid, c.x,   c.y);
    ids[1] = getGridIndex(enhanced_grid, m12.x, m12.y);
    ids[2] = getGridIndex(enhanced_grid, m23.x, m23.y);
    ids[3] = getGridIndex(enhanced_grid, m31.x, m31.y);

    int uniq[4];
#pragma HLS ARRAY_PARTITION variable=uniq complete
    int n = 0;

UNIQ_IDS_MERGED:
    for (int k = 0; k < 4; k++) {
#pragma HLS UNROLL
        int cid = ids[k];
        if (cid < 0 || cid >= enhanced_grid.total_cells) continue;

        bool dup = false;
    CHECK_DUP_MERGED:
        for (int j = 0; j < 4; j++) {
#pragma HLS UNROLL
            if (j >= n) break;
            if (uniq[j] == cid) {
                dup = true;
                break;
            }
        }

        if (!dup) {
            uniq[n] = cid;
            n++;
        }
    }

INSERT_UNIQ_MERGED:
    for (int i = 0; i < 4; i++) {
#pragma HLS UNROLL
        if (i >= n) break;
        tryInsertTriToCell_DDR(grid_counts, grid_triangles,
                               uniq[i], triangleIdx,
                               enhanced_grid.total_cells);
    }
}

// ============================================
// loadGridNeighborhood
// ============================================
void loadGridNeighborhood(int* grid_counts, int* grid_triangles,
                          int center_x, int center_y,
                          int grid_size_x, int grid_size_y) {
#pragma HLS INLINE off

    if (center_x == cached_center_x && center_y == cached_center_y)
        return;

    cached_center_x = center_x;
    cached_center_y = center_y;

    const int R = 2;
    int local_idx = 0;

LOAD_DY:
    for (int dy = -R; dy <= R; dy++) {
        int ty = center_y + dy;

    LOAD_DX:
        for (int dx = -R; dx <= R; dx++) {
            int tx = center_x + dx;
            int lbase = local_idx * MAX_POINT_PER_CELL;

            if (ty < 0 || ty >= grid_size_y || tx < 0 || tx >= grid_size_x) {
                local_counts[local_idx] = 0;
            } else {
                int cell_idx = ty * grid_size_x + tx;
                int count = grid_counts[cell_idx];

                if (count < 0) count = 0;
                if (count > MAX_POINT_PER_CELL) count = MAX_POINT_PER_CELL;

                local_counts[local_idx] = count;

                int ddr_base = cell_idx * MAX_POINT_PER_CELL;

            LOAD_CELL_USED:
                for (int k = 0; k < MAX_POINT_PER_CELL; k++) {
#pragma HLS PIPELINE II=1
                    if (k >= count) break;

                    local_tris[lbase + k] = grid_triangles[ddr_base + k];
                }
            }

            local_idx++;
        }
    }
}

// ============================================
// locateTriangleSimple
// ============================================
#ifndef __SYNTHESIS__
static int total_locate_calls = 0;
static int fallback_searches = 0;
#endif

int locateTriangleSimple(const FixedPoint& p, int triangleCount,
                         int* grid_counts, int* grid_triangles) {
#pragma HLS INLINE off
#ifndef __SYNTHESIS__
    total_locate_calls++;
#endif

    int center_x = max(0, min(enhanced_grid.grid_size_x - 1,
        (int)((p.x - enhanced_grid.x_min) / enhanced_grid.cell_width).to_float()));
    int center_y = max(0, min(enhanced_grid.grid_size_y - 1,
        (int)((p.y - enhanced_grid.y_min) / enhanced_grid.cell_height).to_float()));

    loadGridNeighborhood(grid_counts, grid_triangles,
                         center_x, center_y,
                         enhanced_grid.grid_size_x,
                         enhanced_grid.grid_size_y);

SEARCH_LOCAL: for (int idx = 0; idx < 25; idx++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=25 avg=12
        int count = local_counts[idx];
        int lbase = idx * MAX_POINT_PER_CELL;
        didx_t found_tri = INVALID_INDEX;

    SEARCH_CELL: for (int ti = 0; ti < MAX_POINT_PER_CELL; ti++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=0 max=60 avg=12
            if (ti >= count) continue;
            didx_t tri = (didx_t)local_tris[lbase + ti];

            if (tri == INVALID_INDEX || tri >= triangleCount) continue;

            if (!face_used[tri]) continue;

            if (found_tri == INVALID_INDEX) {
                FaceCoord& fc = face_coord[tri];
                FixedPoint p1; p1.x = fc.x1; p1.y = fc.y1;
                FixedPoint p2; p2.x = fc.x2; p2.y = fc.y2;
                FixedPoint p3; p3.x = fc.x3; p3.y = fc.y3;
                if (isPointInTriangle(p, p1, p2, p3)) {
                    found_tri = tri;
                }
            }
        }

        if (found_tri != INVALID_INDEX) {
            return found_tri;
        }
    }

#ifndef __SYNTHESIS__
    fallback_searches++;
    if (total_locate_calls % 10 == 0) {
        printf("[FALLBACK] calls=%d fallbacks=%d rate=%.1f%%\n",
            total_locate_calls, fallback_searches,
            (float)fallback_searches / total_locate_calls * 100.0f);
    }
#endif

FALLBACK_SEARCH_ALL: for (didx_t tri = 0; tri < MAX_NO_TRIANGLES; tri++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=100 max=14000 avg=5000
    if (tri >= triangleCount) break;

    // face_used==1 already implies v1/v2/v3 are valid.
    if (face_used[tri]) {
        FaceCoord& fc = face_coord[tri];
        FixedPoint p1; p1.x = fc.x1; p1.y = fc.y1;
        FixedPoint p2; p2.x = fc.x2; p2.y = fc.y2;
        FixedPoint p3; p3.x = fc.x3; p3.y = fc.y3;
        if (isPointInTriangle(p, p1, p2, p3)) {
            return tri;
        }
    }
}
    return -1;
}

// ============================================
// connectTwinEdges, inCircle, allocators
// ============================================
void connectTwinEdges(didxh_t e1, didxh_t e2) {
#pragma HLS INLINE
    if (e1 < MAX_NO_HALFEDGES && e2 < MAX_NO_HALFEDGES) {
        he_twin[e1] = e2;
        he_twin[e2] = e1;
    }
}

bool inCircle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1
#pragma HLS LATENCY min=4

    fixed_calc_t a1 = a.x - p.x;
    fixed_calc_t a2 = a.y - p.y;
    fixed_calc_t b1 = b.x - p.x;
    fixed_calc_t b2 = b.y - p.y;
    fixed_calc_t c1 = c.x - p.x;
    fixed_calc_t c2 = c.y - p.y;

    fixed_det_t a1_sq = a1 * a1;
    fixed_det_t a2_sq = a2 * a2;
    fixed_det_t b1_sq = b1 * b1;
    fixed_det_t b2_sq = b2 * b2;
    fixed_det_t c1_sq = c1 * c1;
    fixed_det_t c2_sq = c2 * c2;

    fixed_det_t sum_a = a1_sq + a2_sq;
    fixed_det_t sum_b = b1_sq + b2_sq;
    fixed_det_t sum_c = c1_sq + c2_sq;

    fixed_det_t cross_a = b1 * c2 - b2 * c1;
    fixed_det_t cross_b = a1 * c2 - a2 * c1;
    fixed_det_t cross_c = a1 * b2 - a2 * b1;

    fixed_det_t term1 = sum_a * cross_a;
    fixed_det_t term2 = sum_b * cross_b;
    fixed_det_t term3 = sum_c * cross_c;

    fixed_det_t det = term1 - term2 + term3;
    return det > 0;
}

didx_t allocateFace() {
#pragma HLS INLINE
    if (freeFaceHead != INVALID_INDEX) {
        didx_t result = freeFaceHead;
        freeFaceHead = face_next[result];
        face_next[result] = INVALID_INDEX;
        face_used[result] = 1;
        return result;
    }
    return tri_num++;
}

void freeFace(didx_t faceIdx) {
#pragma HLS INLINE
    if (faceIdx >= MAX_NO_TRIANGLES) return;
    face_used[faceIdx] = 0;
    face_next[faceIdx] = freeFaceHead;
    freeFaceHead = faceIdx;
}

// ============================================
// legalizeEdge
// ============================================
void legalizeEdge(didx_t faceIdx, didxh_t edgeIdx,
                  int* grid_counts, int* grid_triangles,
                  didx_t touched_faces[], int& touched_count) {
#pragma HLS INLINE off

    (void)grid_counts;
    (void)grid_triangles;

    const int MAX_STACK = 30;
    didxh_t edgeStack[MAX_STACK];
    didx_t  faceStack[MAX_STACK];
    const int MAX_FLIPS = 20;
    int flip_count = 0;

    int stackSize = 0;
    if (edgeIdx != INVALID_HALFEDGE && faceIdx != INVALID_INDEX) {
        edgeStack[0] = edgeIdx;
        faceStack[0] = faceIdx;
        stackSize = 1;
    }

LEGALIZE_EDGE_LOOP:
    for (int iter = 0; iter < MAX_FLIPS; iter++) {
#pragma HLS LOOP_TRIPCOUNT min=0 max=20 avg=5
        if (stackSize <= 0) break;

        stackSize--;
        didxh_t currentEdge = edgeStack[stackSize];
        didx_t  currentFace = faceStack[stackSize];

        if (!he_used[currentEdge] || !face_used[currentFace]) continue;

        didxh_t baseCur = face_base_edge(currentFace);
        ap_uint<2> localCur = (ap_uint<2>)(currentEdge - baseCur);

        didxh_t nextEdge = (localCur == 2) ? baseCur : (didxh_t)(currentEdge + 1);
        didxh_t prevEdge = (localCur == 0) ? (didxh_t)(baseCur + 2) : (didxh_t)(currentEdge - 1);

        dnode_t v1 = he_tail[currentEdge];
        dnode_t v2 = he_tail[nextEdge];
        dnode_t v3 = he_tail[prevEdge];

        didxh_t twinEdge = he_twin[currentEdge];
        if (twinEdge == INVALID_HALFEDGE) continue;

        didx_t neighborFace = he_face[twinEdge];
        if (neighborFace == INVALID_INDEX) continue;

        if (!face_used[neighborFace]) continue;

        didxh_t baseNb = face_base_edge(neighborFace);
        ap_uint<2> localNb = (ap_uint<2>)(twinEdge - baseNb);

        didxh_t twinNextEdge = (localNb == 2) ? baseNb : (didxh_t)(twinEdge + 1);
        didxh_t twinPrevEdge = (localNb == 0) ? (didxh_t)(baseNb + 2) : (didxh_t)(twinEdge - 1);

        dnode_t v4 = he_tail[twinPrevEdge];

        FixedPoint p1 = getFixedPoint(v1);
        FixedPoint p2 = getFixedPoint(v2);
        FixedPoint p3 = getFixedPoint(v3);
        FixedPoint p4 = getFixedPoint(v4);

        if (!inCircle(p4, p1, p2, p3)) continue;

        flip_count++;

        didxh_t outerEdgeV1V4 = he_twin[twinNextEdge];
        didxh_t outerEdgeV3V1 = he_twin[prevEdge];
        didxh_t outerEdgeV4V2 = he_twin[twinPrevEdge];
        didxh_t outerEdgeV2V3 = he_twin[nextEdge];

        freeFace(currentFace);
        freeFace(neighborFace);

        // -----------------------------
        // New triangle 1: v1-v4-v3
        // -----------------------------
        didx_t  f1_new = allocateFace();
        didxh_t e1 = face_base_edge(f1_new);
        didxh_t e2 = (didxh_t)(e1 + 1);
        didxh_t e3 = (didxh_t)(e1 + 2);

        face_used[f1_new] = 1;

        he_face[e1] = f1_new;
        he_face[e2] = f1_new;
        he_face[e3] = f1_new;

        he_tail[e1] = v1;
        he_tail[e2] = v4;
        he_tail[e3] = v3;

        he_used[e1] = 1;
        he_used[e2] = 1;
        he_used[e3] = 1;

        he_twin[e1] = INVALID_HALFEDGE;
        he_twin[e2] = INVALID_HALFEDGE;
        he_twin[e3] = INVALID_HALFEDGE;

        if (touched_count < 128) touched_faces[touched_count++] = f1_new;

        // New triangle 2: v4-v2-v3
        // -----------------------------
        didx_t  f2_new = allocateFace();
        didxh_t e4 = face_base_edge(f2_new);
        didxh_t e5 = (didxh_t)(e4 + 1);
        didxh_t e6 = (didxh_t)(e4 + 2);

        face_used[f2_new] = 1;

        he_face[e4] = f2_new;
        he_face[e5] = f2_new;
        he_face[e6] = f2_new;

        he_tail[e4] = v4;
        he_tail[e5] = v2;
        he_tail[e6] = v3;

        he_used[e4] = 1;
        he_used[e5] = 1;
        he_used[e6] = 1;

        he_twin[e4] = INVALID_HALFEDGE;
        he_twin[e5] = INVALID_HALFEDGE;
        he_twin[e6] = INVALID_HALFEDGE;

        if (touched_count < 128) touched_faces[touched_count++] = f2_new;

        connectTwinEdges(e2, e6);
        connectTwinEdges(e1, outerEdgeV1V4);
        connectTwinEdges(e3, outerEdgeV3V1);
        connectTwinEdges(e4, outerEdgeV4V2);
        connectTwinEdges(e5, outerEdgeV2V3);

        // (Removed: vtx_edge[v1..v4] writes — vtx_edge deleted as dead code)

        if (stackSize + 4 <= MAX_STACK) {
            edgeStack[stackSize] = e1; faceStack[stackSize] = f1_new; stackSize++;
            edgeStack[stackSize] = e3; faceStack[stackSize] = f1_new; stackSize++;
            edgeStack[stackSize] = e4; faceStack[stackSize] = f2_new; stackSize++;
            edgeStack[stackSize] = e5; faceStack[stackSize] = f2_new; stackSize++;
        }
    }

}

// ============================================
// insertSiteEnhanced
// ============================================
void insertSiteEnhanced(dnode_t p, didx_t containingFace,
                        int* grid_counts, int* grid_triangles) {
    if (containingFace >= MAX_NO_TRIANGLES) return;

    bool isFirstPoint = (processed_actual_points == 0);
    processed_actual_points++;

    didx_t touched_faces[128];
#pragma HLS ARRAY_PARTITION variable=touched_faces cyclic factor=4
    int touched_count = 0;

    didxh_t orig_e1 = face_base_edge(containingFace);
    didxh_t orig_e2 = (didxh_t)(orig_e1 + 1);
    didxh_t orig_e3 = (didxh_t)(orig_e1 + 2);

    dnode_t v1 = he_tail[orig_e1];
    dnode_t v2 = he_tail[orig_e2];
    dnode_t v3 = he_tail[orig_e3];

    didxh_t outerEdge12 = he_twin[orig_e1];
    didxh_t outerEdge23 = he_twin[orig_e2];
    didxh_t outerEdge31 = he_twin[orig_e3];

    freeFace(containingFace);

    // Triangle 1: p-v1-v2
    didx_t f1 = allocateFace();
    didxh_t e1 = face_base_edge(f1);
    didxh_t e2 = (didxh_t)(e1 + 1);
    didxh_t e3 = (didxh_t)(e1 + 2);

    face_used[f1] = 1;

    he_face[e1] = f1;
    he_face[e2] = f1;
    he_face[e3] = f1;

    he_tail[e1] = p;
    he_tail[e2] = v1;
    he_tail[e3] = v2;
    he_used[e1] = he_used[e2] = he_used[e3] = 1;
    he_twin[e1] = he_twin[e2] = he_twin[e3] = INVALID_HALFEDGE;

    didx_t f2 = allocateFace();
    didxh_t e4 = face_base_edge(f2);
    didxh_t e5 = (didxh_t)(e4 + 1);
    didxh_t e6 = (didxh_t)(e4 + 2);

    face_used[f2] = 1;

    he_face[e4] = f2;
    he_face[e5] = f2;
    he_face[e6] = f2;

    he_tail[e4] = p;
    he_tail[e5] = v2;
    he_tail[e6] = v3;
    he_used[e4] = he_used[e5] = he_used[e6] = 1;
    he_twin[e4] = he_twin[e5] = he_twin[e6] = INVALID_HALFEDGE;

    didx_t f3 = allocateFace();
    didxh_t e7 = face_base_edge(f3);
    didxh_t e8 = (didxh_t)(e7 + 1);
    didxh_t e9 = (didxh_t)(e7 + 2);

    face_used[f3] = 1;

    he_face[e7] = f3;
    he_face[e8] = f3;
    he_face[e9] = f3;

    he_tail[e7] = p;
    he_tail[e8] = v3;
    he_tail[e9] = v1;
    he_used[e7] = he_used[e8] = he_used[e9] = 1;
    he_twin[e7] = he_twin[e8] = he_twin[e9] = INVALID_HALFEDGE;

    connectTwinEdges(e1, e9);
    connectTwinEdges(e4, e3);
    connectTwinEdges(e7, e6);
    connectTwinEdges(e2, outerEdge12);
    connectTwinEdges(e5, outerEdge23);
    connectTwinEdges(e8, outerEdge31);

    touched_faces[touched_count++] = f1;
    touched_faces[touched_count++] = f2;
    touched_faces[touched_count++] = f3;

    if (!isFirstPoint) {
        legalizeEdge(f1, e2, grid_counts, grid_triangles, touched_faces, touched_count);
        legalizeEdge(f2, e5, grid_counts, grid_triangles, touched_faces, touched_count);
        legalizeEdge(f3, e8, grid_counts, grid_triangles, touched_faces, touched_count);
    }

    didx_t unique_faces[128];
    int unique_count = 0;

    DEDUP_TOUCHED:
    for (int i = 0; i < 128; i++) {
    #pragma HLS LOOP_TRIPCOUNT min=3 max=128 avg=24
        if (i >= touched_count) break;

        didx_t f = touched_faces[i];
        if (f >= MAX_NO_TRIANGLES) continue;
        if (!face_used[f]) continue;

        bool dup = false;
    DEDUP_CHECK:
        for (int j = 0; j < 128; j++) {
    #pragma HLS UNROLL
            if (j >= unique_count) break;
            if (unique_faces[j] == f) {
                dup = true;
                break;
            }
        }

        if (!dup) {
            unique_faces[unique_count++] = f;
        }
    }

    FINAL_UPDATE:
    for (int i = 0; i < 128; i++) {
    #pragma HLS LOOP_TRIPCOUNT min=3 max=128 avg=24
        if (i >= unique_count) break;
        didx_t f = unique_faces[i];
        updateFaceCacheAndAssignGrid(f, grid_counts, grid_triangles);
    }

    cached_center_x = -999;
    cached_center_y = -999;
}

// ============================================
// Top-level function: tri2d
// ============================================
void tri2d(hls::stream<axi_i_t>& x_in_stream,
           hls::stream<axi_o_t>& cl_out_stream,
           int* grid_counts,
           int* grid_triangles,
           volatile int& debug_count) {

#pragma HLS INTERFACE axis port=x_in_stream
#pragma HLS INTERFACE axis port=cl_out_stream
#pragma HLS INTERFACE m_axi port=grid_counts    offset=slave bundle=gmem0 \
    depth=1024 max_read_burst_length=64 max_write_burst_length=64
#pragma HLS INTERFACE m_axi port=grid_triangles offset=slave bundle=gmem0 \
    depth=81920 max_read_burst_length=64 max_write_burst_length=64
#pragma HLS INTERFACE s_axilite port=grid_counts    bundle=control
#pragma HLS INTERFACE s_axilite port=grid_triangles bundle=control
#pragma HLS INTERFACE s_axilite port=debug_count    bundle=control
#pragma HLS INTERFACE s_axilite port=return         bundle=control

    // ============================================
    // Storage Mapping Pragmas
    // ============================================
    // HalfEdge
#pragma HLS BIND_STORAGE variable=he_tail  type=ram_1p  impl=lutram
#pragma HLS BIND_STORAGE variable=he_twin  type=ram_s2p impl=bram
#pragma HLS BIND_STORAGE variable=he_used  type=ram_1p  impl=lutram
#pragma HLS BIND_STORAGE variable=he_face  type=ram_s2p impl=bram

    // Face SoA  --- KEY CHANGE: v1/v2/v3 + used in BRAM, next in LUTRAM
#pragma HLS BIND_STORAGE variable=face_used  type=ram_s2p impl=bram
#pragma HLS BIND_STORAGE variable=face_next  type=ram_1p  impl=lutram
#pragma HLS BIND_STORAGE variable=face_v1 type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=face_v2    type=ram_s2p impl=bram
#pragma HLS BIND_STORAGE variable=face_v3    type=ram_s2p impl=bram
#pragma HLS BIND_STORAGE variable=face_coord type=ram_1p  impl=uram

    // Vertex (vtx_x/y kept in URAM to leave BRAM headroom for face_v*)
#pragma HLS BIND_STORAGE variable=vtx_x    type=ram_1p impl=uram
#pragma HLS BIND_STORAGE variable=vtx_y    type=ram_1p impl=uram
#pragma HLS BIND_STORAGE variable=vtx_used type=ram_1p impl=lutram

    // Local cache
#pragma HLS BIND_STORAGE variable=local_counts type=ram_1p impl=bram
#pragma HLS BIND_STORAGE variable=local_tris   type=ram_1p impl=bram

    // ============================================
    // Initialize
    // ============================================

    debug_count = -100;

    tri_num = 0;
    freeFaceHead = INVALID_INDEX;
    processed_actual_points = 0;
    cached_center_x = -999;
    cached_center_y = -999;
#ifndef __SYNTHESIS__
    total_locate_calls = 0;
    fallback_searches = 0;
#endif

    debug_count = -99;

INIT_FACE_TOPO: for (int i = 0; i < MAX_NO_TRIANGLES; i++) {
#pragma HLS PIPELINE II=1
    face_used[i] = 0;
    face_next[i] = INVALID_INDEX;
    face_v1[i]   = 0;
    face_v2[i]   = 0;
    face_v3[i]   = 0;

    face_coord[i].x1 = 0; face_coord[i].y1 = 0;
    face_coord[i].x2 = 0; face_coord[i].y2 = 0;
    face_coord[i].x3 = 0; face_coord[i].y3 = 0;
}

    debug_count = -98;

INIT_HALFEDGE: for (int i = 0; i < MAX_NO_HALFEDGES; i++) {
#pragma HLS PIPELINE II=1
    he_tail[i] = 0;
    he_twin[i] = INVALID_HALFEDGE;
    he_used[i] = 0;
    he_face[i] = INVALID_INDEX;
}

    debug_count = -97;

INIT_VERTEX: for (int i = 0; i < MAX_NO_POINTS; i++) {
#pragma HLS PIPELINE II=1
    vtx_x[i]    = 0;
    vtx_y[i]    = 0;
    vtx_used[i] = 0;
}

    debug_count = -96;

    // ============================================
    // Read input data
    // ============================================
    axi_i_t in_val;

READ_SUPER_TRIANGLE: for (int i = 0; i < 3; i++) {
#pragma HLS PIPELINE II=1
    in_val = x_in_stream.read();
    uint32_t x_bits = in_val.data.range(31, 0);
    uint32_t y_bits = in_val.data.range(63, 32);
    float x_val, y_val;
    memcpy(&x_val, &x_bits, sizeof(float));
    memcpy(&y_val, &y_bits, sizeof(float));
    vtx_x[i] = x_val;
    vtx_y[i] = y_val;
    vtx_used[i] = 1;
}

    debug_count = -93;

READ_NORMAL_POINTS: for (int i = 0; i < num_point; i++) {
#pragma HLS PIPELINE II=1
    in_val = x_in_stream.read();
    uint32_t x_bits = in_val.data.range(31, 0);
    uint32_t y_bits = in_val.data.range(63, 32);
    float x_val, y_val;
    memcpy(&x_val, &x_bits, sizeof(float));
    memcpy(&y_val, &y_bits, sizeof(float));
    vtx_x[i + 3] = x_val;
    vtx_y[i + 3] = y_val;
    vtx_used[i + 3] = 1;
}

    debug_count = -92;

    initEnhancedGrid(enhanced_grid, num_point + 3, grid_counts, grid_triangles);

    debug_count = -91;

    // ============================================
    // Build initial super triangle
    // ============================================
    if (num_point >= 3) {
        didxh_t e0 = (didxh_t)0;
        didxh_t e1 = (didxh_t)1;
        didxh_t e2 = (didxh_t)2;

        face_used[0] = 1;

        he_face[e0] = 0;
        he_face[e1] = 0;
        he_face[e2] = 0;

        he_tail[e0] = 0; he_used[e0] = 1;
        he_tail[e1] = 1; he_used[e1] = 1;
        he_tail[e2] = 2; he_used[e2] = 1;

        // boundary ghost cycle: 3,4,5
        he_tail[3] = 1; he_used[3] = 1;
        he_tail[4] = 2; he_used[4] = 1;
        he_tail[5] = 0; he_used[5] = 1;

        he_face[3] = INVALID_INDEX;
        he_face[4] = INVALID_INDEX;
        he_face[5] = INVALID_INDEX;

        he_twin[e0] = he_twin[e1] = he_twin[e2] = INVALID_HALFEDGE;
        he_twin[3]  = he_twin[4]  = he_twin[5]  = INVALID_HALFEDGE;

        connectTwinEdges(e0, 3);
        connectTwinEdges(e1, 4);
        connectTwinEdges(e2, 5);

        tri_num = 1;
        updateFaceCache(0);
        assignTriangleToGrid(0, grid_counts, grid_triangles);
    }

    debug_count = -90;

    // ============================================
    // Main insertion loop
    // ============================================
INSERT_POINTS: for (int i = 3; i < num_vertice; i++) {
    debug_count = i;
    if (vtx_used[i]) {
        FixedPoint p = getFixedPoint(i);
        int found = locateTriangleSimple(p, tri_num, grid_counts, grid_triangles);
        if (found < 0) continue;
        insertSiteEnhanced(i, (didx_t)found, grid_counts, grid_triangles);
    }
}

    didx_t valid_faces[MAX_FACES];
    int output_count = 0;

COLLECT_INDICES: for (int i = 0; i < MAX_NO_TRIANGLES; i++) {
#pragma HLS PIPELINE II=1
    if (face_used[i]) {
        if (face_v1[i] > 2 && face_v2[i] > 2 && face_v3[i] > 2) {
            valid_faces[output_count++] = i;
        }
    }
}

    // ============================================
    // Output
    // ============================================
OUTPUT_STREAM: for (int i = 0; i < output_count; i++) {
#pragma HLS PIPELINE II=1
    didx_t idx = valid_faces[i];
    dnode_t v1_o = face_v1[idx];
    dnode_t v2_o = face_v2[idx];
    dnode_t v3_o = face_v3[idx];

    axi_o_t val_out;
    val_out.data = 0;
    val_out.data.range(MNPB - 1, 0) = v1_o;
    val_out.data.range(2 * MNPB - 1, MNPB) = v2_o;
    val_out.data.range(3 * MNPB - 1, 2 * MNPB) = v3_o;
    val_out.keep = 0xFF;
    val_out.strb = 0xFF;
    val_out.user = 0;
    val_out.id = 0;
    val_out.dest = 0;
    val_out.last = (i == output_count - 1) ? 1 : 0;

    cl_out_stream.write(val_out);
}

    debug_count = -60;
}
