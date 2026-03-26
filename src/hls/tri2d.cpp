/*
 * KV260 Optimized DCEL-based Delaunay Triangulation
 *
 * Core changes:
 * 1. HalfEdge AoS -> SoA, each field in an independent array -> LUTRAM
 * 2. Face split into FaceTopo (LUTRAM) + FaceCoord (URAM)
 * 3. Vertex split into coordinates (URAM) + topology (LUTRAM)
 * 4. Grid remains in BRAM
 */

#include "tri2d.hpp"
#include <iostream>
#include <cmath>
using namespace std;

// ============================================
// Global Variable Definitions + Storage Mapping Pragmas
// ============================================

didx_t tri_num = 0;
#ifndef __SYNTHESIS__
unsigned int is_cnt = 0;
#endif
int num_point = 1024;
int num_vertice = 1027;
didxh_t global_edge_idx = 6;
static didx_t freeFaceHead = INVALID_INDEX;
#ifndef __SYNTHESIS__
static int inserted_point_count = 0;
#endif
static const int COMPACTION_THRESHOLD = 200;
static int processed_actual_points = 0;

// --------------------------------------------------
// HalfEdge SoA arrays -> LUTRAM (distributed RAM)
//
// Rationale:
//   - Depth 11000, width 12~14 bits, ~154 Kbits per array
//   - Walking algorithm needs to read he_next[e] and he_twin[e] in the same cycle
//   - SoA allows HLS to give each field an independent read port, avoiding arbitration
//   - Depth 11000 may cause HLS to auto-select BRAM; pragma forces LUTRAM
//     If timing closure fails, some arrays can fall back to BRAM
//   - Priority: he_next > he_twin > he_prev > he_tail > he_face > he_used
// --------------------------------------------------
dnode_t  he_tail[MAX_NO_HALFEDGES];
didxh_t  he_twin[MAX_NO_HALFEDGES];
didxh_t  he_prev[MAX_NO_HALFEDGES];
didxh_t  he_next[MAX_NO_HALFEDGES];
didx_t   he_face[MAX_NO_HALFEDGES];
dout_t   he_used[MAX_NO_HALFEDGES];

// --------------------------------------------------
// Face split: topology (LUTRAM) + coordinate cache (URAM)
// --------------------------------------------------
FaceTopo  face_topo[MAX_NO_TRIANGLES];
FaceCoord face_coord[MAX_NO_TRIANGLES];

// --------------------------------------------------
// Vertex SoA: coordinates (URAM) + topology/metadata (LUTRAM)
// --------------------------------------------------
fixed_t     vtx_x[MAX_NO_POINTS];
fixed_t     vtx_y[MAX_NO_POINTS];
hilbert_t   vtx_hilbert[MAX_NO_POINTS];
didxh_t     vtx_edge[MAX_NO_POINTS];
didx_t      vtx_prev[MAX_NO_POINTS];
didx_t      vtx_next[MAX_NO_POINTS];
ap_uint<1>  vtx_used[MAX_NO_POINTS];
ap_uint<1>  vtx_processed[MAX_NO_POINTS];

// Grid (BRAM)
EnhancedAdaptiveGrid enhanced_grid;

// ============================================
// Helper Functions
// ============================================

fixed_t absf(fixed_t val) {
#pragma HLS INLINE
    return (val < 0) ? (fixed_t)(-val) : val;
}

// Hilbert curve encoding
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

// ============================================
// [Updated] getTrianglePoints - uses split arrays
// ============================================
static inline void getTrianglePoints(didx_t tri,
    FixedPoint& c, FixedPoint& m12, FixedPoint& m23, FixedPoint& m31) {
#pragma HLS INLINE
    // Read all coordinates from URAM face_coord in one access
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
// [Updated] updateFaceCache - uses SoA
// ============================================
inline void updateFaceCache(didx_t faceIdx) {
#pragma HLS INLINE
    // Read topology from LUTRAM
    didxh_t e1 = face_topo[faceIdx].edge;
    didxh_t e2 = he_next[e1];
    didxh_t e3 = he_prev[e1];

    dnode_t v1 = he_tail[e1];
    dnode_t v2 = he_tail[e2];
    dnode_t v3 = he_tail[e3];

    // Write vertex indices into face_topo
    face_topo[faceIdx].v1 = v1;
    face_topo[faceIdx].v2 = v2;
    face_topo[faceIdx].v3 = v3;

    // Write coordinate cache into face_coord (URAM)
    face_coord[faceIdx].x1 = vtx_x[v1];
    face_coord[faceIdx].y1 = vtx_y[v1];
    face_coord[faceIdx].x2 = vtx_x[v2];
    face_coord[faceIdx].y2 = vtx_y[v2];
    face_coord[faceIdx].x3 = vtx_x[v3];
    face_coord[faceIdx].y3 = vtx_y[v3];
}

void getTriangleVertices(didx_t faceIdx, dnode_t& v1, dnode_t& v2, dnode_t& v3) {
#pragma HLS INLINE
    didxh_t edgeIdx = face_topo[faceIdx].edge;
    v1 = he_tail[edgeIdx];
    didxh_t nextEdge = he_next[edgeIdx];
    v2 = he_tail[nextEdge];
    didxh_t prevEdge = he_prev[edgeIdx];
    v3 = he_tail[prevEdge];
}

// Point-in-triangle test (unchanged)
bool isPointInTriangle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c) {
    fixed_calc_t cross1 = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
    fixed_calc_t cross2 = (c.x - b.x) * (p.y - b.y) - (c.y - b.y) * (p.x - b.x);
    fixed_calc_t cross3 = (a.x - c.x) * (p.y - c.y) - (a.y - c.y) * (p.x - c.x);
    bool pos = (cross1 >= 0) && (cross2 >= 0) && (cross3 >= 0);
    bool neg = (cross1 <= 0) && (cross2 <= 0) && (cross3 <= 0);
    return pos || neg;
}

// Grid functions
inline int getGridIndex(const EnhancedAdaptiveGrid& grid, fixed_t x, fixed_t y) {
    int grid_x = max(0, min(grid.grid_size_x - 1,
        int((x - grid.x_min) / grid.cell_width)));
    int grid_y = max(0, min(grid.grid_size_y - 1,
        int((y - grid.y_min) / grid.cell_height)));
    return grid_y * grid.grid_size_x + grid_x;
}

// [Updated] initEnhancedGrid - uses vtx_x / vtx_y
void initEnhancedGrid(EnhancedAdaptiveGrid& grid, int num_vertice) {
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

    for (int flat_idx = 0; flat_idx < grid.total_cells; flat_idx++) {
#pragma HLS PIPELINE II=2
        grid.flat_cells[flat_idx].triangle_count = 0;
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

// tryInsertTriToCell - uses face_topo
static inline void tryInsertTriToCell(int grid_idx, didx_t triangleIdx) {
#pragma HLS INLINE
    if (grid_idx < 0 || grid_idx >= enhanced_grid.total_cells) return;
    EnhancedGridCell& cell = enhanced_grid.flat_cells[grid_idx];
    int local_tri_count = cell.triangle_count;

CHECK_EXIST: for (int i = 0; i < MAX_POINT_PER_CELL; i++) {
#pragma HLS UNROLL factor=4
    if (i >= local_tri_count) break;
    if (cell.triangles[i] == triangleIdx) return;
}

if (local_tri_count < MAX_POINT_PER_CELL) {
    cell.triangles[local_tri_count] = triangleIdx;
    cell.triangle_count = local_tri_count + 1;
    return;
}

FIND_INVALID: for (int i = 0; i < MAX_POINT_PER_CELL; i++) {
#pragma HLS PIPELINE II=1
didx_t existing = cell.triangles[i];
if (existing == INVALID_INDEX || existing >= MAX_NO_TRIANGLES || !face_topo[existing].used) {
    cell.triangles[i] = triangleIdx;
    return;
}
}
}

void assignTriangleToGrid(didx_t triangleIdx) {
#pragma HLS INLINE off
    if (triangleIdx >= MAX_NO_TRIANGLES) return;
    if (!face_topo[triangleIdx].used) return;

    FixedPoint c, m12, m23, m31;
    getTrianglePoints(triangleIdx, c, m12, m23, m31);

    int ids[4];
#pragma HLS ARRAY_PARTITION variable=ids complete
    ids[0] = getGridIndex(enhanced_grid, c.x, c.y);
    ids[1] = getGridIndex(enhanced_grid, m12.x, m12.y);
    ids[2] = getGridIndex(enhanced_grid, m23.x, m23.y);
    ids[3] = getGridIndex(enhanced_grid, m31.x, m31.y);

    int uniq[4];
#pragma HLS ARRAY_PARTITION variable=uniq complete
    int n = 0;

    for (int k = 0; k < 4; k++) {
#pragma HLS UNROLL
        int cid = ids[k];
        if (cid < 0 || cid >= enhanced_grid.total_cells) continue;
        bool dup = false;
        for (int j = 0; j < n; j++) {
#pragma HLS UNROLL
            if (uniq[j] == cid) { dup = true; break; }
        }
        if (!dup) uniq[n++] = cid;
    }

    for (int i = 0; i < n; i++) {
#pragma HLS UNROLL
        tryInsertTriToCell(uniq[i], triangleIdx);
    }
}

// ============================================
// Point Location - uses face_topo + face_coord
// ============================================
#ifndef __SYNTHESIS__
static int global_hits = 0;
static int total_hits = 0;
static int total_locate_calls = 0;
static int fallback_searches = 0;
#endif

int locateTriangleSimple(const FixedPoint& p, int triangleCount) {
#pragma HLS INLINE off
#ifndef __SYNTHESIS__
    total_locate_calls++;
#endif

    int center_x = max(0, min(enhanced_grid.grid_size_x - 1,
        (int)((p.x - enhanced_grid.x_min) / enhanced_grid.cell_width).to_float()));
    int center_y = max(0, min(enhanced_grid.grid_size_y - 1,
        (int)((p.y - enhanced_grid.y_min) / enhanced_grid.cell_height).to_float()));

    const int search_range = 5;

SEARCH_DX_LOOP: for (int dx = -search_range; dx <= search_range; dx++) {
SEARCH_DY_LOOP: for (int dy = -search_range; dy <= search_range; dy++) {
#pragma HLS LOOP_FLATTEN

    int target_x = center_x + dx;
    int target_y = center_y + dy;

    if (target_x >= 0 && target_x < enhanced_grid.grid_size_x &&
        target_y >= 0 && target_y < enhanced_grid.grid_size_y) {

        int cell_idx = target_y * enhanced_grid.grid_size_x + target_x;
        EnhancedGridCell& cell = enhanced_grid.flat_cells[cell_idx];
        int old_count = cell.triangle_count;
        int write_idx = 0;
        didx_t found_tri = INVALID_INDEX;

    SEARCH_AND_COMPACT: for (int ti = 0; ti < MAX_POINT_PER_CELL; ti++) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=0 max=60 avg=12
        if (ti >= old_count) continue;
        didx_t tri = cell.triangles[ti];

        if (tri == INVALID_INDEX || tri >= triangleCount || !face_topo[tri].used) {
            continue;
        }

        cell.triangles[write_idx] = tri;
        write_idx++;

        if (found_tri == INVALID_INDEX) {
            // Read coordinate cache from URAM
            FaceCoord& fc = face_coord[tri];
            FixedPoint p1; p1.x = fc.x1; p1.y = fc.y1;
            FixedPoint p2; p2.x = fc.x2; p2.y = fc.y2;
            FixedPoint p3; p3.x = fc.x3; p3.y = fc.y3;

            if (isPointInTriangle(p, p1, p2, p3)) {
                found_tri = tri;
            }
        }
    }

    cell.triangle_count = write_idx;

    if (found_tri != INVALID_INDEX) {
#ifndef __SYNTHESIS__
        total_hits++;
#endif
        return found_tri;
    }
    }
}
}

// Fallback: linear search over all triangles
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
if (face_topo[tri].used) {
    FaceTopo& ft = face_topo[tri];
    if (ft.v1 >= 0 && ft.v2 >= 0 && ft.v3 >= 0) {
        FaceCoord& fc = face_coord[tri];
        FixedPoint p1; p1.x = fc.x1; p1.y = fc.y1;
        FixedPoint p2; p2.x = fc.x2; p2.y = fc.y2;
        FixedPoint p3; p3.x = fc.x3; p3.y = fc.y3;
        if (isPointInTriangle(p, p1, p2, p3)) {
#ifndef __SYNTHESIS__
            global_hits++;
            total_hits++;
#endif
            return tri;
        }
    }
}
}
return -1;
}

// ============================================
// Twin edge connection - uses he_twin SoA
// ============================================
void connectTwinEdges(didxh_t e1, didxh_t e2) {
#pragma HLS INLINE
    if (e1 < MAX_NO_HALFEDGES && e2 < MAX_NO_HALFEDGES) {
        he_twin[e1] = e2;
        he_twin[e2] = e1;
    }
}

bool inCircle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c) {
#pragma HLS INLINE off      // Never inline: prevents flattening into oversized combinational logic
#pragma HLS PIPELINE II=1   // Enable pipelining: HLS may use multiple cycles but accepts new data every cycle
#pragma HLS LATENCY min=4

#ifndef __SYNTHESIS__
    static int call_count = 0;
    call_count++;
#endif

    // Stage 1: compute relative coordinates (subtraction)
    fixed_calc_t a1 = a.x - p.x;
    fixed_calc_t a2 = a.y - p.y;
    fixed_calc_t b1 = b.x - p.x;
    fixed_calc_t b2 = b.y - p.y;
    fixed_calc_t c1 = c.x - p.x;
    fixed_calc_t c2 = c.y - p.y;

    // Stage 2: compute squared terms (multiplication)
    fixed_det_t a1_sq = a1 * a1;
    fixed_det_t a2_sq = a2 * a2;
    fixed_det_t b1_sq = b1 * b1;
    fixed_det_t b2_sq = b2 * b2;
    fixed_det_t c1_sq = c1 * c1;
    fixed_det_t c2_sq = c2 * c2;

    // Stage 3: compute sum of squares and cross products
    fixed_det_t sum_a = a1_sq + a2_sq;
    fixed_det_t sum_b = b1_sq + b2_sq;
    fixed_det_t sum_c = c1_sq + c2_sq;

    fixed_det_t cross_a = b1 * c2 - b2 * c1;
    fixed_det_t cross_b = a1 * c2 - a2 * c1;
    fixed_det_t cross_c = a1 * b2 - a2 * b1;

    // Stage 4: core multiplications (decomposed)
    fixed_det_t term1 = sum_a * cross_a;
    fixed_det_t term2 = sum_b * cross_b;
    fixed_det_t term3 = sum_c * cross_c;

    // Stage 5: final determinant
    fixed_det_t det = term1 - term2 + term3;

    return det > 0;
}

// ============================================
// Edge/Face allocation - uses SoA
// ============================================
const int EDGE_POOL_SIZE = 480;
static didxh_t freeEdges[EDGE_POOL_SIZE];
static int freeEdgeCount = 0;
static int replaceIdx = 0;

didxh_t allocateHalfEdge(didxh_t& newEdgeIdx) {
    if (freeEdgeCount > 0) {
        return freeEdges[--freeEdgeCount];
    }
    return newEdgeIdx++;
}

void freeHalfEdge(didxh_t edgeId) {
    he_used[edgeId] = 0;
    if (freeEdgeCount < EDGE_POOL_SIZE) {
        freeEdges[freeEdgeCount++] = edgeId;
    }
    else {
        freeEdges[replaceIdx] = edgeId;
        replaceIdx = (replaceIdx + 1) % EDGE_POOL_SIZE;
    }
}

didx_t allocateFace() {
#pragma HLS INLINE
    if (freeFaceHead != INVALID_INDEX) {
        didx_t result = freeFaceHead;
        freeFaceHead = face_topo[result].next;
        face_topo[result].next = INVALID_INDEX;
        face_topo[result].used = 1;
        return result;
    }
    return tri_num++;
}

void freeFace(didx_t faceIdx) {
#pragma HLS INLINE
    if (faceIdx >= MAX_NO_TRIANGLES) return;
    face_topo[faceIdx].used = 0;
    face_topo[faceIdx].edge = INVALID_HALFEDGE;
    face_topo[faceIdx].prev = INVALID_INDEX;
    face_topo[faceIdx].next = freeFaceHead;
    freeFaceHead = faceIdx;
}

// ============================================
// [Updated] legalizeEdge - uses SoA
// ============================================
void legalizeEdge(didx_t faceIdx, didxh_t edgeIdx, didxh_t& newEdgeIdx) {
#pragma HLS INLINE off

    const int MAX_STACK = 30;
    didxh_t edgeStack[MAX_STACK];
    didx_t faceStack[MAX_STACK];
    const int MAX_FLIPS = 20;
    int flip_count = 0;
#ifndef __SYNTHESIS__
    static int total_flips = 0;
    static int max_flips_per_call = 0;
    static int collected_triangles = 0;
#endif

    int stackSize = 0;
    if (edgeIdx != INVALID_HALFEDGE && faceIdx != INVALID_INDEX) {
        edgeStack[0] = edgeIdx;
        faceStack[0] = faceIdx;
        stackSize = 1;
    }

LEGALIZE_EDGE_LOOP: for (int iter = 0; iter < MAX_FLIPS; iter++) {
#pragma HLS LOOP_TRIPCOUNT min=0 max=20 avg=5
    if (stackSize <= 0) break;

    stackSize--;
    didxh_t currentEdge = edgeStack[stackSize];
    didx_t currentFace = faceStack[stackSize];

    if (!he_used[currentEdge] || !face_topo[currentFace].used) continue;

    // Read half-edge topology from SoA
    dnode_t v1 = he_tail[currentEdge];
    didxh_t nextEdge = he_next[currentEdge];
    dnode_t v2 = he_tail[nextEdge];
    didxh_t prevEdge = he_prev[currentEdge];
    dnode_t v3 = he_tail[prevEdge];

    didxh_t twinEdge = he_twin[currentEdge];
    if (twinEdge == INVALID_HALFEDGE) continue;

    didx_t neighborFace = he_face[twinEdge];
    if (neighborFace == INVALID_INDEX || !face_topo[neighborFace].used) continue;

    didxh_t twinNextEdge = he_next[twinEdge];
    didxh_t twinPrevEdge = he_prev[twinEdge];
    dnode_t v4 = he_tail[twinPrevEdge];

    FixedPoint p1 = getFixedPoint(v1);
    FixedPoint p2 = getFixedPoint(v2);
    FixedPoint p3 = getFixedPoint(v3);
    FixedPoint p4 = getFixedPoint(v4);

    if (!inCircle(p4, p1, p2, p3)) continue;

    flip_count++;
#ifndef __SYNTHESIS__
    total_flips++;
#endif

    // Save outer twin edges before freeing
    didxh_t outerEdgeV1V4 = he_twin[twinNextEdge];
    didxh_t outerEdgeV3V1 = he_twin[prevEdge];
    didxh_t outerEdgeV4V2 = he_twin[twinPrevEdge];
    didxh_t outerEdgeV2V3 = he_twin[nextEdge];

    freeFace(currentFace);
    freeFace(neighborFace);
    freeHalfEdge(currentEdge);
    freeHalfEdge(nextEdge);
    freeHalfEdge(prevEdge);
    freeHalfEdge(twinEdge);
    freeHalfEdge(twinNextEdge);
    freeHalfEdge(twinPrevEdge);

    // New triangle 1: v1-v4-v3
    didx_t f1_new = allocateFace();
    didxh_t e1 = allocateHalfEdge(newEdgeIdx);
    didxh_t e2 = allocateHalfEdge(newEdgeIdx);
    didxh_t e3 = allocateHalfEdge(newEdgeIdx);

    face_topo[f1_new].edge = e1;
    face_topo[f1_new].used = 1;

    // Write SoA half-edge fields
    he_tail[e1] = v1; he_face[e1] = f1_new; he_next[e1] = e2; he_prev[e1] = e3; he_used[e1] = 1;
    he_tail[e2] = v4; he_face[e2] = f1_new; he_next[e2] = e3; he_prev[e2] = e1; he_used[e2] = 1;
    he_tail[e3] = v3; he_face[e3] = f1_new; he_next[e3] = e1; he_prev[e3] = e2; he_used[e3] = 1;

    // New triangle 2: v4-v2-v3
    didx_t f2_new = allocateFace();
    didxh_t e4 = allocateHalfEdge(newEdgeIdx);
    didxh_t e5 = allocateHalfEdge(newEdgeIdx);
    didxh_t e6 = allocateHalfEdge(newEdgeIdx);

    face_topo[f2_new].edge = e4;
    face_topo[f2_new].used = 1;

    he_tail[e4] = v4; he_face[e4] = f2_new; he_next[e4] = e5; he_prev[e4] = e6; he_used[e4] = 1;
    he_tail[e5] = v2; he_face[e5] = f2_new; he_next[e5] = e6; he_prev[e5] = e4; he_used[e5] = 1;
    he_tail[e6] = v3; he_face[e6] = f2_new; he_next[e6] = e4; he_prev[e6] = e5; he_used[e6] = 1;

    connectTwinEdges(e2, e6);
    connectTwinEdges(e1, outerEdgeV1V4);
    connectTwinEdges(e3, outerEdgeV3V1);
    connectTwinEdges(e4, outerEdgeV4V2);
    connectTwinEdges(e5, outerEdgeV2V3);

    vtx_edge[v1] = e1;
    vtx_edge[v2] = e5;
    vtx_edge[v3] = e6;
    vtx_edge[v4] = e4;

    assignTriangleToGrid(f1_new);
    assignTriangleToGrid(f2_new);
    updateFaceCache(f1_new);
    updateFaceCache(f2_new);
#ifndef __SYNTHESIS__
    collected_triangles += 2;
#endif

    if (stackSize + 4 <= MAX_STACK) {
        edgeStack[stackSize] = e1; faceStack[stackSize] = f1_new; stackSize++;
        edgeStack[stackSize] = e3; faceStack[stackSize] = f1_new; stackSize++;
        edgeStack[stackSize] = e4; faceStack[stackSize] = f2_new; stackSize++;
        edgeStack[stackSize] = e5; faceStack[stackSize] = f2_new; stackSize++;
    }
}

#ifndef __SYNTHESIS__
if (flip_count > max_flips_per_call) max_flips_per_call = flip_count;
#endif
}

// ============================================
// [Updated] insertSiteEnhanced - uses SoA
// ============================================
void insertSiteEnhanced(dnode_t p, didx_t containingFace) {
    if (containingFace >= MAX_NO_TRIANGLES) return;
    FixedPoint newPoint = getFixedPoint(p);
#ifndef __SYNTHESIS__
    static int insert_count = 0;
    insert_count++;
#endif

    bool isFirstPoint = (processed_actual_points == 0);
    processed_actual_points++;
    didxh_t newEdgeIdx = global_edge_idx;

    // Read original triangle topology from SoA
    didxh_t orig_e1 = face_topo[containingFace].edge;
    didxh_t orig_e2 = he_next[orig_e1];
    didxh_t orig_e3 = he_prev[orig_e1];

    dnode_t v1 = he_tail[orig_e1];
    dnode_t v2 = he_tail[orig_e2];
    dnode_t v3 = he_tail[orig_e3];

    freeFace(containingFace);

    didxh_t outerEdge12 = he_twin[orig_e1];
    didxh_t outerEdge23 = he_twin[orig_e2];
    didxh_t outerEdge31 = he_twin[orig_e3];

    freeHalfEdge(orig_e1);
    freeHalfEdge(orig_e2);
    freeHalfEdge(orig_e3);

    // === Triangle 1: p-v1-v2 ===
    didx_t f1 = allocateFace();
    didxh_t e1 = allocateHalfEdge(newEdgeIdx);
    didxh_t e2 = allocateHalfEdge(newEdgeIdx);
    didxh_t e3 = allocateHalfEdge(newEdgeIdx);

    face_topo[f1].edge = e1;
    face_topo[f1].used = 1;

    he_tail[e1] = p;  he_face[e1] = f1; he_next[e1] = e2; he_prev[e1] = e3; he_used[e1] = 1;
    he_tail[e2] = v1; he_face[e2] = f1; he_next[e2] = e3; he_prev[e2] = e1; he_used[e2] = 1;
    he_tail[e3] = v2; he_face[e3] = f1; he_next[e3] = e1; he_prev[e3] = e2; he_used[e3] = 1;

    // === Triangle 2: p-v2-v3 ===
    didx_t f2 = allocateFace();
    didxh_t e4 = allocateHalfEdge(newEdgeIdx);
    didxh_t e5 = allocateHalfEdge(newEdgeIdx);
    didxh_t e6 = allocateHalfEdge(newEdgeIdx);

    face_topo[f2].edge = e4;
    face_topo[f2].used = 1;

    he_tail[e4] = p;  he_face[e4] = f2; he_next[e4] = e5; he_prev[e4] = e6; he_used[e4] = 1;
    he_tail[e5] = v2; he_face[e5] = f2; he_next[e5] = e6; he_prev[e5] = e4; he_used[e5] = 1;
    he_tail[e6] = v3; he_face[e6] = f2; he_next[e6] = e4; he_prev[e6] = e5; he_used[e6] = 1;

    // === Triangle 3: p-v3-v1 ===
    didx_t f3 = allocateFace();
    didxh_t e7 = allocateHalfEdge(newEdgeIdx);
    didxh_t e8 = allocateHalfEdge(newEdgeIdx);
    didxh_t e9 = allocateHalfEdge(newEdgeIdx);

    face_topo[f3].edge = e7;
    face_topo[f3].used = 1;

    he_tail[e7] = p;  he_face[e7] = f3; he_next[e7] = e8; he_prev[e7] = e9; he_used[e7] = 1;
    he_tail[e8] = v3; he_face[e8] = f3; he_next[e8] = e9; he_prev[e8] = e7; he_used[e8] = 1;
    he_tail[e9] = v1; he_face[e9] = f3; he_next[e9] = e7; he_prev[e9] = e8; he_used[e9] = 1;

    // Connect internal twin pairs
    connectTwinEdges(e1, e9);
    connectTwinEdges(e4, e3);
    connectTwinEdges(e7, e6);

    // Connect external twin pairs
    connectTwinEdges(e2, outerEdge12);
    connectTwinEdges(e5, outerEdge23);
    connectTwinEdges(e8, outerEdge31);

    assignTriangleToGrid(f1);
    assignTriangleToGrid(f2);
    assignTriangleToGrid(f3);
    updateFaceCache(f1);
    updateFaceCache(f2);
    updateFaceCache(f3);

    if (!isFirstPoint) {
        legalizeEdge(f1, e2, newEdgeIdx);
        legalizeEdge(f2, e5, newEdgeIdx);
        legalizeEdge(f3, e8, newEdgeIdx);
    }
    else {
#ifndef __SYNTHESIS__
        cout << "First point insertion completed, skip edge flipping" << endl;
#endif
    }
    global_edge_idx = newEdgeIdx;
}

// ============================================
// Top-level function: tri2d
// ============================================
void tri2d(hls::stream<axi_i_t>& x_in_stream,
    hls::stream<axi_o_t>& cl_out_stream,
    volatile int& debug_count) {
#pragma HLS INTERFACE axis port=x_in_stream
#pragma HLS INTERFACE axis port=cl_out_stream
#pragma HLS INTERFACE s_axilite port=debug_count bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    // ============================================
    // Storage Mapping Pragmas
    // ============================================

    // HalfEdge SoA -> LUTRAM (distributed RAM)
    // Note: depth 11000 is large for LUTRAM and will consume significant LUTs.
    // If post-synthesis LUT utilization exceeds 80%, move he_face back to BRAM.
    // Priority: he_next > he_twin > he_prev > he_tail > he_face > he_used
#pragma HLS BIND_STORAGE variable=he_tail  type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=he_twin  type=ram_s2p impl=bram
#pragma HLS BIND_STORAGE variable=he_prev  type=ram_s2p impl=bram
#pragma HLS BIND_STORAGE variable=he_next  type=ram_s2p impl=bram
#pragma HLS BIND_STORAGE variable=he_face  type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=he_used  type=ram_1p impl=lutram

    // Face topology -> LUTRAM
#pragma HLS BIND_STORAGE variable=face_topo type=ram_1p impl=lutram

    // Face coordinate cache -> URAM (large capacity, wide data, read-heavy)
#pragma HLS BIND_STORAGE variable=face_coord type=ram_1p impl=uram

    // Vertex coordinates -> URAM
#pragma HLS BIND_STORAGE variable=vtx_x type=ram_1p impl=uram
#pragma HLS BIND_STORAGE variable=vtx_y type=ram_1p impl=uram

    // Vertex topology/metadata -> LUTRAM (depth 2003 is ideal)
#pragma HLS BIND_STORAGE variable=vtx_edge      type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=vtx_used       type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=vtx_processed  type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=vtx_prev       type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=vtx_next       type=ram_1p impl=lutram
#pragma HLS BIND_STORAGE variable=vtx_hilbert    type=ram_1p impl=lutram

    // Grid -> BRAM (default; explicit pragma for documentation clarity)
#pragma HLS BIND_STORAGE variable=enhanced_grid.flat_cells type=ram_1p impl=bram

    // freeEdges small array -> LUTRAM
#pragma HLS BIND_STORAGE variable=freeEdges type=ram_1p impl=lutram

    // ============================================
    // Initialize global variables
    // ============================================
    debug_count = -100;  // entered function

    tri_num = 0;
#ifndef __SYNTHESIS__
    is_cnt = 0;
    inserted_point_count = 0;
#endif
    freeFaceHead = INVALID_INDEX;
    replaceIdx = 0;
    global_edge_idx = 6;
    freeEdgeCount = 0;
    processed_actual_points = 0;

    debug_count = -99;  // scalar init done

    // Initialize face topology and coordinate cache
INIT_FACE_TOPO: for (int i = 0; i < MAX_NO_TRIANGLES; i++) {
#pragma HLS PIPELINE II=1
    face_topo[i].edge = INVALID_HALFEDGE;
    face_topo[i].used = 0;
    face_topo[i].prev = INVALID_INDEX;
    face_topo[i].next = INVALID_INDEX;
    face_coord[i].x1 = 0; face_coord[i].y1 = 0;
    face_coord[i].x2 = 0; face_coord[i].y2 = 0;
    face_coord[i].x3 = 0; face_coord[i].y3 = 0;
}

debug_count = -98;  // face init done

// Initialize HalfEdge SoA arrays
INIT_HALFEDGE: for (int i = 0; i < MAX_NO_HALFEDGES; i++) {
#pragma HLS PIPELINE II=1
he_tail[i] = 0;
he_twin[i] = INVALID_HALFEDGE;
he_prev[i] = INVALID_HALFEDGE;
he_next[i] = INVALID_HALFEDGE;
he_face[i] = INVALID_INDEX;
he_used[i] = 0;
}

debug_count = -97;  // halfedge init done

// Initialize Vertex SoA arrays
INIT_VERTEX: for (int i = 0; i < MAX_NO_POINTS; i++) {
#pragma HLS PIPELINE II=1
vtx_x[i] = 0;
vtx_y[i] = 0;
vtx_used[i] = 0;
vtx_edge[i] = INVALID_HALFEDGE;
vtx_prev[i] = INVALID_INDEX;
vtx_next[i] = INVALID_INDEX;
vtx_processed[i] = 0;
vtx_hilbert[i] = 0;
}

debug_count = -96;  // vertex init done

// Initialize grid cells
INIT_GRID: for (int i = 0; i < MAX_GRID_SIZE * MAX_GRID_SIZE; i++) {
#pragma HLS PIPELINE II=1
enhanced_grid.flat_cells[i].triangle_count = 0;
}

debug_count = -95;  // grid init done

// Initialize free edge pool
INIT_FREE_EDGES: for (int i = 0; i < EDGE_POOL_SIZE; i++) {
#pragma HLS PIPELINE II=1
freeEdges[i] = INVALID_HALFEDGE;
}

debug_count = -94;  // freeEdges init done

// ============================================
// Read input data
// ============================================
axi_i_t in_val;
axi_o_t val_out;

// Read super-triangle vertices
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

debug_count = -93;  // super triangle read done

// Read regular input points
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

debug_count = -92;  // all points read done

// ============================================
// Initialize grid and super-triangle
// ============================================
initEnhancedGrid(enhanced_grid, num_point + 3);

debug_count = -91;  // grid setup done

if (num_point >= 3) {
    face_topo[0].edge = 0;
    face_topo[0].used = 1;

    // Super-triangle half-edges (SoA style)
    he_tail[0] = 0; he_face[0] = 0; he_next[0] = 1; he_prev[0] = 2; he_used[0] = 1;
    he_tail[1] = 1; he_face[1] = 0; he_next[1] = 2; he_prev[1] = 0; he_used[1] = 1;
    he_tail[2] = 2; he_face[2] = 0; he_next[2] = 0; he_prev[2] = 1; he_used[2] = 1;

    he_tail[3] = 1; he_face[3] = INVALID_INDEX; he_next[3] = 4; he_prev[3] = 5; he_used[3] = 1;
    he_tail[4] = 2; he_face[4] = INVALID_INDEX; he_next[4] = 5; he_prev[4] = 3; he_used[4] = 1;
    he_tail[5] = 0; he_face[5] = INVALID_INDEX; he_next[5] = 3; he_prev[5] = 4; he_used[5] = 1;

    connectTwinEdges(0, 3);
    connectTwinEdges(1, 4);
    connectTwinEdges(2, 5);

    tri_num = 1;
    global_edge_idx = 6;
    updateFaceCache(0);
}

debug_count = -90;  // super triangle setup done, entering INSERT_POINTS

// ============================================
// Main point insertion loop
// ============================================
INSERT_POINTS: for (int i = 3; i < num_vertice; i++) {
    debug_count = i;  // update before processing each point
    if (vtx_used[i]) {
        FixedPoint p = getFixedPoint(i);
        int found = locateTriangleSimple(p, tri_num);
        if (found < 0) continue;
        insertSiteEnhanced(i, (didx_t)found);
    }
}

didx_t valid_faces[MAX_FACES];
int output_count = 0;

COLLECT_INDICES: for (int i = 0; i < MAX_NO_TRIANGLES; i++) {
#pragma HLS PIPELINE II=1
if (face_topo[i].used && face_topo[i].edge != INVALID_HALFEDGE) {
    if (face_topo[i].v1 > 2 && face_topo[i].v2 > 2 && face_topo[i].v3 > 2) {
        valid_faces[output_count++] = i;
    }
}
}

// ============================================
// Output results
// ============================================
OUTPUT_STREAM: for (int i = 0; i < output_count; i++) {
#pragma HLS PIPELINE II=1
didx_t idx = valid_faces[i];
FaceTopo& ft = face_topo[idx];

axi_o_t val_out;
val_out.data = 0;
val_out.data.range(MNPB - 1, 0) = ft.v1;
val_out.data.range(2 * MNPB - 1, MNPB) = ft.v2;
val_out.data.range(3 * MNPB - 1, 2 * MNPB) = ft.v3;
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