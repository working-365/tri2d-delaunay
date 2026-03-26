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
#define MAX_NO_POINTS 1027
#define MAX_NO_TRIANGLES 8000
#define MAX_FACES 8000
#define MAX_NO_HALFEDGES 20000
#define MAX_POINT_PER_CELL 40
#define MAX_GRID_SIZE 32
#define HILBERT_ORDER 8

// Bit widths
#define MNPB 13  // Point index bit width

// ============================================
// Type Definitions
// ============================================
typedef ap_axis<64, 16, 5, 6> axi_i_t;
typedef ap_axis<64, 16, 5, 6> axi_o_t;
typedef ap_uint<1> dout_t;
typedef ap_uint<13> dnode_t;
typedef ap_uint<14> didx_t;
typedef ap_uint<16> didxh_t;
typedef ap_uint<32> hilbert_t;

typedef ap_fixed<24, 4> fixed_t;
typedef ap_fixed<33, 6> fixed_calc_t;
typedef ap_fixed<48, 10> fixed_det_t;

const dnode_t INVALID_NODE = (1 << 14) - 1;
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
// [Core Change 1] HalfEdge: AoS -> SoA (Structure of Arrays)
//
// Before: AoS (Array of Structures)
//   HalfEdge he[11000];  // 68 bits each, HLS maps as a monolithic BRAM block
//
// After: Split into independent arrays, each mappable separately
// Benefits:
//   - Each field array: depth 11000, width 12~14 bits
//   - HLS can map each array independently, allowing multi-field reads in one cycle
//   - With ram_style=distributed, consumes LUTs instead of BRAMs
//     If timing closure fails, some arrays can fall back to BRAM
// ============================================

// HalfEdge SoA field declarations (defined in .cpp)
extern dnode_t  he_tail[MAX_NO_HALFEDGES];     // 12-bit x 11000
extern didxh_t  he_twin[MAX_NO_HALFEDGES];     // 14-bit x 11000
extern didxh_t  he_prev[MAX_NO_HALFEDGES];     // 14-bit x 11000
extern didxh_t  he_next[MAX_NO_HALFEDGES];     // 14-bit x 11000
extern didx_t   he_face[MAX_NO_HALFEDGES];     // 13-bit x 11000
extern dout_t   he_used[MAX_NO_HALFEDGES];     // 1-bit  x 11000

// ============================================
// [Core Change 2] Face: Split into topology + coordinate cache
//
// Topology (FaceTopo): frequent read/write -> LUTRAM
//   - edge, used, next, prev, v1, v2, v3
//   - ~14+1+13+13+12x3 = 77 bits per entry
//   - 5000 x 77 = 385,000 bits ~= 47 KB -> fits in LUTRAM
//
// Coordinate cache (FaceCoord): large but read-heavy -> URAM
//   - x1,y1,x2,y2,x3,y3 (6 x 24 = 144 bits per entry)
//   - 5000 x 144 = 720,000 bits ~= 88 KB -> fits comfortably in one URAM
// ============================================

struct FaceTopo {
    didxh_t edge;
    ap_uint<1> used;
    didx_t next, prev;
    dnode_t v1, v2, v3;
};

struct FaceCoord {
    fixed_t x1, y1;
    fixed_t x2, y2;
    fixed_t x3, y3;
};

extern FaceTopo face_topo[MAX_NO_TRIANGLES];
extern FaceCoord face_coord[MAX_NO_TRIANGLES];

// ============================================
// [Core Change 3] Vertex: Separate coordinates from topology
//
// Coordinate arrays: read-heavy -> URAM (packed as 48-bit: x[24]+y[24])
// Topology arrays: small, suitable for LUTRAM or BRAM
// ============================================

// Vertex coordinate SoA (-> URAM)
extern fixed_t vtx_x[MAX_NO_POINTS];          // 24-bit x 2003
extern fixed_t vtx_y[MAX_NO_POINTS];          // 24-bit x 2003

// Vertex topology/metadata SoA (-> LUTRAM, depth 2003 is ideal)
extern hilbert_t    vtx_hilbert[MAX_NO_POINTS];  // 32-bit x 2003
extern didxh_t      vtx_edge[MAX_NO_POINTS];     // 14-bit x 2003
extern didx_t       vtx_prev[MAX_NO_POINTS];     // 13-bit x 2003
extern didx_t       vtx_next[MAX_NO_POINTS];     // 13-bit x 2003
extern ap_uint<1>   vtx_used[MAX_NO_POINTS];     // 1-bit  x 2003
extern ap_uint<1>   vtx_processed[MAX_NO_POINTS]; // 1-bit x 2003

// ============================================
// Grid structure unchanged (-> BRAM)
// Medium depth (1024), medium width: BRAM is the best fit
// ============================================

struct EnhancedGridCell {
    didx_t triangles[MAX_POINT_PER_CELL];
    int triangle_count;

    EnhancedGridCell() {
#pragma HLS INLINE
        triangle_count = 0;
    }
};

struct EnhancedAdaptiveGrid {
    fixed_t x_min, y_min;
    fixed_t x_max, y_max;
    fixed_t cell_width, cell_height;
    int grid_size_x, grid_size_y;
    int total_cells;

    EnhancedGridCell flat_cells[MAX_GRID_SIZE * MAX_GRID_SIZE];

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

void initEnhancedGrid(EnhancedAdaptiveGrid& grid, int num_vertices);
int getGridIndex(const EnhancedAdaptiveGrid& grid, fixed_t x, fixed_t y);

didxh_t allocateHalfEdge(didxh_t& newEdgeIdx);
void freeHalfEdge(didxh_t edgeId);
didx_t allocateFace();
void freeFace(didx_t faceIdx);

int locateTriangleSimple(const FixedPoint& p, int triangleCount);
void connectTwinEdges(didxh_t e1, didxh_t e2);
void legalizeEdge(didx_t faceIdx, didxh_t edgeIdx, didxh_t& newEdgeIdx);
void insertSiteEnhanced(dnode_t p, didx_t containingFace);

FixedPoint calculateTriangleCentroid(didx_t triangleIdx);
void assignTriangleToGrid(didx_t triangleIdx);

// Helper inline function: read FixedPoint from SoA
inline FixedPoint getFixedPoint(dnode_t pointIdx) {
#pragma HLS INLINE
    FixedPoint p;
    p.x = vtx_x[pointIdx];
    p.y = vtx_y[pointIdx];
    return p;
}

void tri2d(hls::stream<axi_i_t>& x_in_stream,
    hls::stream<axi_o_t>& cl_out_stream,
    volatile int& debug_count);
#endif // TRI2D_HPP
