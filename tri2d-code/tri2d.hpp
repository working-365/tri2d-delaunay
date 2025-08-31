/**
 * @file tri2d.hpp
 * @brief 2D Delaunay Triangulation (Tri2D) IP core with Hilbert curve-based spatial indexing
 *
 * Header file for Delaunay triangulation implementation with:
 * - Dynamic grid sizing based on point count square root
 * - Hilbert curve ordering for spatial coherence
 * - Grid-based point location with Hilbert locality search
 *
 * This code was originally based on repository
 * https://github.com/PranjalSahu/Point-Cloud-Triangulation
 * Enhanced by adding DCEL structure and Hilbert curve-based spatial indexing
 */

#ifndef TRI2D_HPP
#define TRI2D_HPP

#include <ap_fixed.h>
#include <hls_stream.h>
#include <ap_int.h>
#include <ap_axi_sdata.h>


 // Definitions and constants
#ifndef M_PI
#define M_PI	3.14159265358979323846  /* pi */
#endif

#define DEG_TO_RAD(angle_deg) ((angle_deg) * M_PI / 180.0) // 角度转弧度宏
#define EPSILON 0.0001  // 浮点数比较容��?

// 网格容量常量 - 基于5000个点和约4万个三角��?
#define MAX_NO_POINTS 5003      // ��?大点数（顶点数）
#define MAX_VERTICES 5000       // 同MAX_POINTS，保持一致�??
#define MAX_FACES 12000        // ��?大面数（大约为顶点数��?8-9倍）
#define MAX_NO_HALFEDGES 33000 // ��?大半边数（约��?3*MAX_FACES��?
#define MAX_POINT_PER_CELL 25 // 每个网格单元��?大点��?/三角形数
#define MAX_NO_TRIANGLES 12000  // 支持接近4万个三角��?


// 算法特定常量
#define MAX_GRID_SIZE 64       // 网格��?大尺寸（每个维度��?
#define MAX_HILBERT_ORDER 8    // Hilbert曲线的最大阶��?

// 位宽计算常量 - 用于硬件实现
// 计算clog2 - 向上取整的log2，用于确定位��?
#define MNPB 12  // MAX_NO_POINTS的位��?

// 数组访问安全边界
#define SAFE_MARGIN 18        // 为防止越界提供的额外安全边界

#define VISUAL
#define FIXPT

#ifdef FIXPT
typedef ap_axis<64, 16, 5, 6> axi_i_t;   // 保持不变，数据位宽足��?
typedef ap_axis<64, 16, 5, 6> axi_o_t;   // 保持不变
typedef ap_uint<1> dout_t;               // 保持不变，布尔�??
typedef ap_uint<1> dchg_t;               // 保持不变，布尔�??
typedef ap_uint<13> dnode_t;             // 支持��?��?8,191个点（�?�合5000个点��?
typedef ap_uint<14> didx_t;              // 三角形索引，支持��?��?32,768个三角形
typedef ap_uint<16> didxh_t;              // 半边索引，支持最��?131,072个半��?
typedef ap_uint<32> hilbert_t;           // 保持不变，足够表示Hilbert编码

// 定点数计算类��? - 增加精度以处理更大规模的计算
typedef ap_fixed<24, 4> fixed_t;        // 增加整数部分��?16位，以容纳左��?12位后的�??
typedef ap_fixed<33, 6> fixed_calc_t;   // 中间计算精度相应增加
typedef ap_fixed<48, 10> fixed_det_t;    // 行列式计算精度相应增��?

#else
typedef unsigned int dout_t;
typedef unsigned int dchg_t;
typedef unsigned int dnode_t;
typedef unsigned int dnum_t;
typedef unsigned int didx_t;
typedef unsigned int didxh_t;
typedef unsigned long long hilbert_t;
#endif

// 修正无效值的计算方式
const dnode_t INVALID_NODE = (1 << 11) - 1;   // 8191
const didx_t INVALID_INDEX = (1 << 12) - 1;   // 32767
const didxh_t INVALID_HALFEDGE = (1 << 14) - 1; // 131071
const int INVALID_CELL = -1;  // 保持使用-1作为普�?�int的无效�??

// 用于几何计算的浮点点结构
struct FloatPoint {
    float x, y;
    dnode_t index;
    hilbert_t hilbert_code;  // Hilbert曲线编码

    // 默认构�?�函��?
    FloatPoint() {
#pragma HLS INLINE
        x = 0.0f;
        y = 0.0f;
        index = INVALID_NODE;  // 使用常量
        hilbert_code = 0;
    }

    // 添加列表初始化构造函��?
    FloatPoint(float _x, float _y, dnode_t _index) {
#pragma HLS INLINE
        x = _x;
        y = _y;
        index = _index;
        hilbert_code = 0;
    }
};

// 定点数点结构��?
struct FixedPoint {
    fixed_t x, y;
    dnode_t index;


    FixedPoint() {
        x = 0;
        y = 0;
        index = INVALID_NODE;

    }

    FixedPoint(fixed_t _x, fixed_t _y, dnode_t _index) {
        x = _x;
        y = _y;
        index = _index;

    }
};



// 半边结构定义
struct HalfEdge {
    dnode_t tail;     // 半边尾部顶点索引
    didxh_t twin;     // 对偶半边索引 - 修改为didxh_t
    didxh_t prev;     // 边界中前��?条半边的索引 - 修改为didxh_t
    didxh_t next;     // 边界中下��?条半边的索引 - 修改为didxh_t
    didx_t face;      // 半边左侧的面索引 - 保持didx_t
    dout_t used;      // 使用标志

    HalfEdge() {
#pragma HLS INLINE
        tail = (dnode_t)0;
        twin = INVALID_HALFEDGE;  // 使用新常��?
        prev = INVALID_HALFEDGE;  // 使用新常��?
        next = INVALID_HALFEDGE;  // 使用新常��?
        face = INVALID_INDEX;     // 使用三角形索引常��?
        used = (dout_t)0;
    }
};

struct Face {
    didxh_t edge;     // 面边界上的一条半边索��? - 修改为didxh_t
    didx_t prev;      // 链表中前��?个面的索��? - 保持didx_t
    didx_t next;      // 链表中下��?个面的索��? - 保持didx_t
    dout_t used;      // 使用标志

    Face() {
#pragma HLS INLINE
        edge = INVALID_HALFEDGE;  // 使用新常��?
        prev = INVALID_INDEX;     // 使用三角形索引常��?
        next = INVALID_INDEX;     // 使用三角形索引常��?
        used = (dout_t)0;
    }
};

struct Vertex {
    fixed_t x, y;
    hilbert_t hilbert_code;  // Hilbert曲线编码
    didxh_t edge;    // 修改为didxh_t
    didx_t prev, next;
    ap_uint<1> used;
    ap_uint<1> processed;

    Vertex() {
        x = (dnode_t)0;
        y = (dnode_t)0;
        hilbert_code = 0;
        edge = INVALID_HALFEDGE;  // 使用新常��?
        prev = INVALID_INDEX;     // 使用三角形索引常��?
        next = INVALID_INDEX;     // 使用三角形索引常��?
        used = 0;
        processed = 0;
    }
};

// 网格单元的点信息 - 用于排序
struct PointWithZ {
    dnode_t index;         // 点的索引
    hilbert_t z_value;     // Hilbert曲线值（空间填充曲线编码��?
};



// 在EnhancedGridCell结构中添加三角形数组和计数器
struct EnhancedGridCell {
    // 顶点数组
    dnode_t vertices[MAX_POINT_PER_CELL];
    int vertex_count;
    
    // 新增：三角形数组
    didx_t triangles[MAX_POINT_PER_CELL]; 
    int triangle_count;

    EnhancedGridCell() {
#pragma HLS INLINE
        vertex_count = 0;
        triangle_count = 0;
    }
};

// 转换后的拍平网格结构
struct EnhancedAdaptiveGrid {
    fixed_t x_min, y_min;
    fixed_t x_max, y_max;
    fixed_t cell_width, cell_height;
    int grid_size_x, grid_size_y;
    int total_cells; // 总单元格数量

    // 拍平的一维单元格数组
    EnhancedGridCell flat_cells[MAX_GRID_SIZE * MAX_GRID_SIZE];

    EnhancedAdaptiveGrid() {
#pragma HLS INLINE
        x_min = y_min = fixed_t(1e30f);  // 使用定点数最大�??
        x_max = y_max = fixed_t(-1e30f); // 使用定点数最小�??
        cell_width = cell_height = 0;
        grid_size_x = grid_size_y = 0;
        total_cells = 0;
    }

    // 二维索引转一维索引的辅助函数
    inline int getFlatIndex(int x, int y) const {
        return y * grid_size_x + x;
    }

    // ��?维索引转二维坐标的辅助函��?
    inline void getGridCoords(int flat_index, int& x, int& y) const {
        x = flat_index % grid_size_x;
        y = flat_index / grid_size_x;
    }

    // 获取指定坐标的单元格（兼容旧代码��?
    inline EnhancedGridCell& getCell(int x, int y) {
        return flat_cells[getFlatIndex(x, y)];
    }

    // 通过线�?�索引获取单元格
    inline EnhancedGridCell& getCellByIndex(int idx) {
        return flat_cells[idx];
    }
};

// 全局变量声明
extern Face face[MAX_NO_TRIANGLES];
extern HalfEdge he[MAX_NO_HALFEDGES];
extern Vertex vertex[MAX_NO_POINTS];
extern EnhancedAdaptiveGrid enhanced_grid;
extern didx_t tri_num;
extern unsigned int is_cnt;
extern int num_point;

// 绝对值函数声��?
fixed_t absf(fixed_t val);


// 几何工具函数声明
bool inCircle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c);

// 三角形顶点获取函数声��?
void getTriangleVertices(didx_t faceIdx, dnode_t& v1, dnode_t& v2, dnode_t& v3);
fixed_t orientTest(const FixedPoint& a, const FixedPoint& b, const FixedPoint& p);

// Hilbert曲线函数声明
hilbert_t hilbert_xy2d(int n, int x, int y);  // 新增：Hilbert XY到D映射函数
hilbert_t getHilbertCode(fixed_t x, fixed_t y, fixed_t x_min, fixed_t y_min, 
                         fixed_t x_max, fixed_t y_max, int order);  // 改：从Z曲线到Hilbert曲线

// 保留morton_encode用于兼容性或其他用�??
hilbert_t morton_encode(int x, int y);
// 保留getZCurveCode用于兼容性，但最好�?�步替换为Hilbert版本
hilbert_t getZCurveCode(fixed_t x, fixed_t y, fixed_t x_min, fixed_t y_min, 
                        fixed_t x_max, fixed_t y_max, int order);
// 网格函数声明
void initEnhancedGrid(EnhancedAdaptiveGrid& grid, const Vertex* vertices, int num_vertices);
int getTrianglesFromPoint(dnode_t pointIdx, didx_t* triangleList, int maxTriangles);
int getGridIndex(const EnhancedAdaptiveGrid& grid, fixed_t x, fixed_t y);
didxh_t allocateHalfEdge(didxh_t& newEdgeIdx) ;
void freeHalfEdge(didxh_t edgeId);
void compactTriangleIndices();



// 三角形定位函数声��?
bool isPointInTriangle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c);
int locateTriangleSimple(const FixedPoint& p, int triangleCount);  // 添加此函数声��?

// DCEL操作函数声明
void connectTwinEdges(didxh_t e1, didxh_t e2);  // 修改为didxh_t
void legalizeEdge(didx_t faceIdx, didxh_t edgeIdx, didx_t& newFaceIdx, didxh_t& newEdgeIdx);  // 混合使用两种类型
void insertSiteEnhanced(dnode_t p);

// 新增：三角形网格管理函数
FixedPoint calculateTriangleCentroid(didx_t triangleIdx);
void assignTriangleToGrid(didx_t triangleIdx);

// Delaunay属�?�检查函��?
bool checkDelaunayProperty();  // 添加此函数声��?

// 主函数声��? - 移除了不再使用的edge_mat参数
void tri2d(hls::stream<axi_i_t>& x_in_stream,
    hls::stream<axi_o_t>& cl_out_stream);

#endif // TRI2D_HPP
