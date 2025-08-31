/*
 * This code implements a DCEL-based Delaunay triangulation algorithm with:
 * - Dynamic grid sizing based on point count square root
 * - Grid-based point location with simplified search
 *
 * This code was originally based on repository
 * https://github.com/PranjalSahu/Point-Cloud-Triangulation
 * Enhanced by adding DCEL structure and spatial indexing
 */

#include "tri2d.hpp"
#include <iostream>
#include <cmath>
using namespace std;

// 全局变量定义
didx_t tri_num = 0; // 或其他初始�??; // 三角形数�????
unsigned int is_cnt = 0;
Face face[MAX_NO_TRIANGLES];
HalfEdge he[MAX_NO_HALFEDGES];
Vertex vertex[MAX_NO_POINTS];
EnhancedAdaptiveGrid enhanced_grid;
int num_point = 5000;  // 实际点数（假设与输入数据�????致）
int num_vertice = 5003;
didxh_t global_edge_idx = 6; 
static int inserted_point_count = 0;
static const int COMPACTION_THRESHOLD = 200;

// 定点版本
fixed_t absf(fixed_t val) {
#pragma HLS INLINE
    if (val < 0)
        return -val;
    else
        return val;
}

// 优化getFixedPoint函数 - 添加管道化和强制内联
FixedPoint getFixedPoint(dnode_t pointIdx) {
    FixedPoint p;
    if (pointIdx >= 0 && pointIdx < MAX_NO_POINTS) {
        // 直接使用坐标，不再进行缩�????
        p.x = vertex[pointIdx].x;
        p.y = vertex[pointIdx].y;
    }
    else {
        p.x = fixed_t(0);
        p.y = fixed_t(0);
    }
    return p;
}

// Hilbert曲线编码计算函数 - 为点计算Hilbert�????
hilbert_t getHilbertCode(fixed_t x, fixed_t y, fixed_t x_min, fixed_t y_min,
                     fixed_t x_max, fixed_t y_max, int order) {
    
    
    // 转换为浮点数进行计算
    float x_f = x.to_float();
    float y_f = y.to_float();
    float xmin_f = x_min.to_float();
    float ymin_f = y_min.to_float();
    float xmax_f = x_max.to_float();
    float ymax_f = y_max.to_float();

    // 归一化坐标到[0,1]范围
    float norm_x = (x_f - xmin_f) / (xmax_f - xmin_f);
    float norm_y = (y_f - ymin_f) / (ymax_f - ymin_f);

    // 限制范围
    norm_x = std::max(0.0f, std::min(1.0f, norm_x));
    norm_y = std::max(0.0f, std::min(1.0f, norm_y));

    // 映射到网�????
    int grid_size = 1 << order;
    int hx = (int)(norm_x * (grid_size - 1));
    int hy = (int)(norm_y * (grid_size - 1));

    // 计算Hilbert曲线编码
    hilbert_t code = hilbert_xy2d(order, hx, hy);

    return code;
}

// Hilbert曲线的XY到D（一维）映射
hilbert_t hilbert_xy2d(int n, int x, int y) {
#pragma HLS INLINE

    hilbert_t d = 0;
    for (int s = 1; s < (1 << n); s *= 2) {
#pragma HLS PIPELINE
        int rx = (x & s) > 0;
        int ry = (y & s) > 0;
        
        d += s * s * ((3 * rx) ^ ry);
        
        // 旋转
        if (ry == 0) {
            if (rx == 1) {
                x = s - 1 - x;
                y = s - 1 - y;
            }
            // 交换x和y
            int t = x;
            x = y;
            y = t;
        }
    }
    return d;
}

// 便捷函数 - �????次�?�获取三个顶�????
void getTriangleVertices(didx_t faceIdx, dnode_t& v1, dnode_t& v2, dnode_t& v3) {
#pragma HLS INLINE
    // 获取半边索引
    didxh_t edgeIdx = face[faceIdx].edge;

    // 获取三个顶点
    v1 = he[edgeIdx].tail;
    didxh_t nextEdge = he[edgeIdx].next;
    v2 = he[nextEdge].tail;
    didxh_t nextEdge2 = he[nextEdge].next;
    v3 = he[nextEdge2].tail;
}

// 点在三角形内判断
bool isPointInTriangle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c) {
    // 计算三个叉积
    fixed_calc_t cross1 = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
    fixed_calc_t cross2 = (c.x - b.x) * (p.y - b.y) - (c.y - b.y) * (p.x - b.x);
    fixed_calc_t cross3 = (a.x - c.x) * (p.y - c.y) - (a.y - c.y) * (p.x - c.x);

    // 判断三个叉积的符号是否相�????
    bool pos = (cross1 >= 0) && (cross2 >= 0) && (cross3 >= 0);
    bool neg = (cross1 <= 0) && (cross2 <= 0) && (cross3 <= 0);

    return pos || neg;
}

// 计算点所在网格的线�?�索�???? (而不是坐�????)
inline int getGridIndex(const EnhancedAdaptiveGrid& grid, fixed_t x, fixed_t y) {
    // 计算网格坐标，确保在有效范围�????
    int grid_x = max(0, min(grid.grid_size_x - 1,
        int((x - grid.x_min) / grid.cell_width)));
    int grid_y = max(0, min(grid.grid_size_y - 1,
        int((y - grid.y_min) / grid.cell_height)));

    // 直接返回线�?�索�????
    return grid_y * grid.grid_size_x + grid_x;
}

void initEnhancedGrid(EnhancedAdaptiveGrid& grid, const Vertex* vertices, int num_vertice) {
    #pragma HLS INLINE off

    // Find bounds
    grid.x_min = fixed_t(1e30f);
    grid.y_min = fixed_t(1e30f);
    grid.x_max = fixed_t(-1e30f);
    grid.y_max = fixed_t(-1e30f);

    // Only use non-initial triangle vertices (index >= 3)
    for (int i = 3; i < num_vertice; i++) {
        if (vertices[i].used) {
            fixed_t x = vertices[i].x;
            fixed_t y = vertices[i].y;
            
            grid.x_min = (x < grid.x_min) ? x : grid.x_min;
            grid.y_min = (y < grid.y_min) ? y : grid.y_min;
            grid.x_max = (x > grid.x_max) ? x : grid.x_max;
            grid.y_max = (y > grid.y_max) ? y : grid.y_max;
        }
    }

    // Add margin
    fixed_calc_t margin_x = (grid.x_max - grid.x_min) * fixed_calc_t(0.2);
    fixed_calc_t margin_y = (grid.y_max - grid.y_min) * fixed_calc_t(0.2);
    grid.x_min -= fixed_t(margin_x);
    grid.y_min -= fixed_t(margin_y);
    grid.x_max += fixed_t(margin_x);
    grid.y_max += fixed_t(margin_y);

    // Set grid size
    grid.grid_size_x = MAX_GRID_SIZE;
    grid.grid_size_y = MAX_GRID_SIZE;
    grid.total_cells = grid.grid_size_x * grid.grid_size_y;

    // Calculate cell size
    float x_range = (grid.x_max - grid.x_min).to_float();
    float y_range = (grid.y_max - grid.y_min).to_float();
    grid.cell_width = fixed_t(x_range / grid.grid_size_x);
    grid.cell_height = fixed_t(y_range / grid.grid_size_y);

    // Initialize cells
    for (int flat_idx = 0; flat_idx < grid.total_cells; flat_idx++) {
        #pragma HLS PIPELINE II=2
        grid.flat_cells[flat_idx].vertex_count = 0;
        grid.flat_cells[flat_idx].triangle_count = 0;
    }
}

// 计算三角形重�????
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

void assignTriangleToGrid(didx_t triangleIdx) {
    #pragma HLS INLINE off
    
    // 先检查三角形自身是否有效
    if (!face[triangleIdx].used) {
        return;  // 跳过未使用的三角�??
    }
    
    // 计算重心
    FixedPoint centroid = calculateTriangleCentroid(triangleIdx);
    
    // 计算网格索引
    int grid_idx = getGridIndex(enhanced_grid, centroid.x, centroid.y);
    
    // �??查索引是否有�??
    if (grid_idx < 0 || grid_idx >= enhanced_grid.total_cells) {
        return;
    }
    
    // 获取网格单元引用
    EnhancedGridCell& cell = enhanced_grid.flat_cells[grid_idx];
    
    // 1. 首先�??查单元格是否未满，如果未满直接添加到末尾
    if (cell.triangle_count < MAX_POINT_PER_CELL) {
        cell.triangles[cell.triangle_count] = triangleIdx;
        cell.triangle_count++;
        return;
    }
    
    // 2. 只有当单元格已满时，才尝试替换无效三角形
    for (int i = 0; i < cell.triangle_count; i++) {
        didx_t existing_tri = cell.triangles[i];
        
        // �??查三角形是否无效
        if (!face[existing_tri].used) {
            // 替换无效三角�??
            cell.triangles[i] = triangleIdx;
            return;
        }
    }
    
    // 3. 如果未能添加（单元格已满且无无效三角形），直接返�??
}


// Global counters to track statistics

// ȫ�ֻ�̬������
static int global_hits = 0;
static int total_hits  = 0;
static int total_locate_calls = 0;
static int fallback_searches = 0;

int locateTriangleSimple(const FixedPoint& p, int triangleCount) {
    #pragma HLS INLINE off
    
    total_locate_calls++; // ͳ���ܵ��ô���

    int center_x = max(0, min(enhanced_grid.grid_size_x - 1,
                   (int)((p.x - enhanced_grid.x_min) / enhanced_grid.cell_width).to_float()));
    int center_y = max(0, min(enhanced_grid.grid_size_y - 1,
                   (int)((p.y - enhanced_grid.y_min) / enhanced_grid.cell_height).to_float()));
    
    const int search_range = 2;
    
    // Search nearby grid cells
    for (int dx = -search_range; dx <= search_range; dx++) {
        for (int dy = -search_range; dy <= search_range; dy++) {
            int target_x = center_x + dx;
            int target_y = center_y + dy;
            
            if (target_x >= 0 && target_x < enhanced_grid.grid_size_x &&
                target_y >= 0 && target_y < enhanced_grid.grid_size_y) {
                
                int cell_idx = target_y * enhanced_grid.grid_size_x + target_x;
                EnhancedGridCell* current_cell = &enhanced_grid.flat_cells[cell_idx];
                
                for (int ti = 0; ti < current_cell->triangle_count; ti++) {
                    didx_t tri = current_cell->triangles[ti];
                    
                    if (tri != INVALID_INDEX && face[tri].used) {
                        dnode_t v1, v2, v3;
                        getTriangleVertices(tri, v1, v2, v3);
                        
                        FixedPoint p1 = getFixedPoint(v1);
                        FixedPoint p2 = getFixedPoint(v2);
                        FixedPoint p3 = getFixedPoint(v3);
                        
                        if (isPointInTriangle(p, p1, p2, p3)) {
                            total_hits++;
                            return tri;
                        }
                    }
                }
            }
        }
    }
    
    // ����fallbackǰͳ��
    fallback_searches++;

    // Fallback: search all triangles
    for (didx_t tri = 0; tri < triangleCount; tri++) {
        if (face[tri].used) {
            dnode_t v1, v2, v3;
            getTriangleVertices(tri, v1, v2, v3);
            if (v1 >= 0 && v2 >= 0 && v3 >= 0) {
                FixedPoint p1 = getFixedPoint(v1);
                FixedPoint p2 = getFixedPoint(v2);
                FixedPoint p3 = getFixedPoint(v3);
                
                if (isPointInTriangle(p, p1, p2, p3)) {
                    global_hits++;
                    total_hits++;
                    // �����������ӡͳ�����?
                    if (total_locate_calls % 100 == 0) {
                        printf("Fallback ratio: %.2f%% (%d/%d)\n",
                               (double)fallback_searches / total_locate_calls * 100.0,
                               fallback_searches, total_locate_calls);
                    }
                    return tri;
                }
            }
        }
    }
    
    return -1;
}


// 连接两条半边作为对偶
void connectTwinEdges(didxh_t e1, didxh_t e2) {
#pragma HLS INLINE 
    // �????查边索引有效�????
    if (e1 < MAX_NO_HALFEDGES && e2 < MAX_NO_HALFEDGES) {
        he[e1].twin = e2;
        he[e2].twin = e1;
    }
}

// 点在外接圆内测试
bool inCircle(const FixedPoint& p, const FixedPoint& a, const FixedPoint& b, const FixedPoint& c) {
    // 添加静�?�计数器
    static int call_count = 0;
    call_count++;  // 增加计数

    // 使用行列式方法计�????
    fixed_calc_t a1 = a.x - p.x;
    fixed_calc_t a2 = a.y - p.y;
    fixed_calc_t b1 = b.x - p.x;
    fixed_calc_t b2 = b.y - p.y;
    fixed_calc_t c1 = c.x - p.x;
    fixed_calc_t c2 = c.y - p.y;

    // 计算平方�????(只计算一�????)
    fixed_det_t a1_sq = a1 * a1;
    fixed_det_t a2_sq = a2 * a2;
    fixed_det_t b1_sq = b1 * b1;
    fixed_det_t b2_sq = b2 * b2;
    fixed_det_t c1_sq = c1 * c1;
    fixed_det_t c2_sq = c2 * c2;

    // 使用预计算的平方�????
    fixed_det_t det =
        (a1_sq + a2_sq) * (b1 * c2 - b2 * c1) -
        (b1_sq + b2_sq) * (a1 * c2 - a2 * c1) +
        (c1_sq + c2_sq) * (a1 * b2 - a2 * b1);

    // 结果比较
    bool result = det > 0;

    return result;
}



// �????查整个三角剖分是否满足Delaunay性质（DCEL版本�????
 /*bool checkDelaunayProperty() {
    int invalid_triangles = 0;
    int valid_triangles = 0;

    // 遍历�????有有效的三角�????
    for (didx_t i = 0; i < tri_num; i++) {
        if (!face[i].used) continue;

        valid_triangles++;

        // 获取三角形的三条�????
        didxh_t e1 = face[i].edge;
        didxh_t e2 = he[e1].next;
        didxh_t e3 = he[e2].next;

        // 获取三个顶点
        dnode_t v1 = he[e1].tail;
        dnode_t v2 = he[e2].tail;
        dnode_t v3 = he[e3].tail;

        // 获取顶点坐标
        FixedPoint p1 = getFixedPoint(v1);
        FixedPoint p2 = getFixedPoint(v2);
        FixedPoint p3 = getFixedPoint(v3);

        // �????查其他所有点是否在这个三角形的外接圆�????
        bool violates_delaunay = false;

        // 直接遍历�????有顶�????
        for (dnode_t test_vertex = 0; test_vertex < MAX_NO_POINTS; test_vertex++) {
            // 跳过未使用的顶点
            if (!vertex[test_vertex].used) continue;

            // 跳过三角形自身的顶点
            if (test_vertex == v1 || test_vertex == v2 || test_vertex == v3) {
                continue;
            }

            FixedPoint test_point = getFixedPoint(test_vertex);

            // 使用inCircle函数�????查点是否在外接圆�????
            if (inCircle(test_point, p1, p2, p3)) {
                violates_delaunay = true;
                break;  // 找到�????个违反的点就足够�????
            }
        }

        if (violates_delaunay) {
            invalid_triangles++;
        }
    }

    // 只打印统计结�????
    cout << "aDelaunay=" << valid_triangles
         << ", bDelaunay=" << invalid_triangles
         << ", 结果=" << (invalid_triangles == 0 ? "通过" : "失败") << endl;

    return (invalid_triangles == 0);
}*/

// === 全局复用池定�??? ===
const int EDGE_POOL_SIZE = 480;
static didxh_t freeEdges[EDGE_POOL_SIZE];
static int freeEdgeCount = 0;
static int replaceIdx = 0;


// === 全局半边管理函数 ===
didxh_t allocateHalfEdge(didxh_t& newEdgeIdx) {
    if (freeEdgeCount > 0) {
        return freeEdges[--freeEdgeCount];
    }
    return newEdgeIdx++;
}

void freeHalfEdge(didxh_t edgeId) {
    he[edgeId].used = 0;
    if (freeEdgeCount < EDGE_POOL_SIZE) {
        freeEdges[freeEdgeCount++] = edgeId;
    } else {
        freeEdges[replaceIdx] = edgeId;
        replaceIdx = (replaceIdx + 1) % EDGE_POOL_SIZE;
    }
}

void compactTriangleIndices() {
    #pragma HLS INLINE off
    
    didx_t index_map[MAX_NO_TRIANGLES];
    didx_t new_tri_count = 0;
    
    // Initialize mapping
    for (didx_t i = 0; i < MAX_NO_TRIANGLES; i++) {
        index_map[i] = INVALID_INDEX;
    }
    
    // Build mapping
    for (didx_t old_idx = 0; old_idx < tri_num; old_idx++) {
        if (face[old_idx].used) {
            index_map[old_idx] = new_tri_count++;
        }
    }
    
    // Compact face array
    didx_t write_ptr = 0;
    for (didx_t read_ptr = 0; read_ptr < tri_num; read_ptr++) {
        if (face[read_ptr].used) {
            if (read_ptr != write_ptr) {
                face[write_ptr] = face[read_ptr];
                face[read_ptr].used = 0;
            }
            write_ptr++;
        }
    }
    
    // Update halfedge face references
    for (didxh_t he_idx = 0; he_idx < global_edge_idx; he_idx++) {
        if (he[he_idx].used && he[he_idx].face != INVALID_INDEX) {
            didx_t old_face_idx = he[he_idx].face;
            if (old_face_idx < tri_num && index_map[old_face_idx] != INVALID_INDEX) {
                he[he_idx].face = index_map[old_face_idx];
            } else {
                he[he_idx].face = INVALID_INDEX;
            }
        }
    }
    
    // Update grid triangle references
    for (int cell_idx = 0; cell_idx < enhanced_grid.total_cells; cell_idx++) {
        EnhancedGridCell& cell = enhanced_grid.flat_cells[cell_idx];
        int write_pos = 0;
        
        for (int read_pos = 0; read_pos < cell.triangle_count; read_pos++) {
            didx_t old_tri_idx = cell.triangles[read_pos];
            
            if (old_tri_idx < tri_num && index_map[old_tri_idx] != INVALID_INDEX) {
                cell.triangles[write_pos++] = index_map[old_tri_idx];
            }
        }
        
        cell.triangle_count = write_pos;
    }
    
    tri_num = new_tri_count;
}


void legalizeEdge(didx_t faceIdx, didxh_t edgeIdx, didxh_t& newEdgeIdx) {
    #pragma HLS INLINE off
    
    // 用于存储待检查的�????
    const int MAX_STACK = 30;
    didxh_t edgeStack[MAX_STACK];
    didx_t faceStack[MAX_STACK];

    // 翻转次数限制
    const int MAX_FLIPS = 20;   
    int flip_count = 0;         
    
    // 统计计数�????
    static int total_flips = 0;
    static int max_flips_per_call = 0;
    static int collected_triangles = 0;
    
    // 初始化栈
    int stackSize = 0;
    if (edgeIdx != INVALID_HALFEDGE && faceIdx != INVALID_INDEX) {
        edgeStack[0] = edgeIdx;
        faceStack[0] = faceIdx;
        stackSize = 1;
    }

    // 处理栈中�????有边
    while (stackSize > 0 && flip_count < MAX_FLIPS) {
        // 弹出栈顶�????
        stackSize--;
        didxh_t currentEdge = edgeStack[stackSize];
        didx_t currentFace = faceStack[stackSize];

        // �????查有效�??
        if (!he[currentEdge].used || !face[currentFace].used)
            continue;

        // 获取三角形顶�????
        dnode_t v1 = he[currentEdge].tail;
        didxh_t nextEdge = he[currentEdge].next;
        dnode_t v2 = he[nextEdge].tail;
        didxh_t prevEdge = he[currentEdge].prev;
        dnode_t v3 = he[prevEdge].tail;

        // 获取对偶边和邻接三角�????
        didxh_t twinEdge = he[currentEdge].twin;
        if (twinEdge == INVALID_HALFEDGE)
            continue;

        didx_t neighborFace = he[twinEdge].face;
        if (neighborFace == INVALID_INDEX || !face[neighborFace].used)
            continue;

        // 获取第四个顶�????
        didxh_t twinNextEdge = he[twinEdge].next;
        didxh_t twinPrevEdge = he[twinEdge].prev;
        dnode_t v4 = he[twinPrevEdge].tail;

        // 获取各顶点坐�????
        FixedPoint p1 = getFixedPoint(v1);
        FixedPoint p2 = getFixedPoint(v2);
        FixedPoint p3 = getFixedPoint(v3);
        FixedPoint p4 = getFixedPoint(v4);

        // �????查是否需要翻�???? - Delaunay条件
        if (!inCircle(p4, p1, p2, p3))
            continue;

        // 增加翻转计数
        flip_count++;
        total_flips++;

        // 保存外部边的对偶关系
        didxh_t outerEdgeV1V4 = he[twinNextEdge].twin;
        didxh_t outerEdgeV3V1 = he[prevEdge].twin;
        didxh_t outerEdgeV4V2 = he[twinPrevEdge].twin;
        didxh_t outerEdgeV2V3 = he[nextEdge].twin;

        // 标记原三角形为未使用
       face[currentFace].used = 0;
       face[neighborFace].used = 0;

        freeHalfEdge(currentEdge);
        freeHalfEdge(nextEdge);
        freeHalfEdge(prevEdge);
        freeHalfEdge(twinEdge);
        freeHalfEdge(twinNextEdge);
        freeHalfEdge(twinPrevEdge);

        // 创建两个新三角形和对应的半边
        // 第一个三角形: v1-v4-v3
        didx_t f1_new = tri_num++ ;
        
        // 分配半边
        didxh_t e1 = allocateHalfEdge(newEdgeIdx);
        didxh_t e2 = allocateHalfEdge(newEdgeIdx);
        didxh_t e3 = allocateHalfEdge(newEdgeIdx);

        face[f1_new].edge = e1;
        face[f1_new].used = 1;

        he[e1].tail = v1;  he[e1].face = f1_new;  he[e1].next = e2;  he[e1].prev = e3;  he[e1].used = 1;
        he[e2].tail = v4;  he[e2].face = f1_new;  he[e2].next = e3;  he[e2].prev = e1;  he[e2].used = 1;
        he[e3].tail = v3;  he[e3].face = f1_new;  he[e3].next = e1;  he[e3].prev = e2;  he[e3].used = 1;

        // 第二个三角形: v4-v2-v3
        didx_t f2_new = tri_num++ ; 
        
        // 分配半边
        didxh_t e4 = allocateHalfEdge(newEdgeIdx);
        didxh_t e5 = allocateHalfEdge(newEdgeIdx);
        didxh_t e6 = allocateHalfEdge(newEdgeIdx);

        face[f2_new].edge = e4;
        face[f2_new].used = 1;

        he[e4].tail = v4;  he[e4].face = f2_new;  he[e4].next = e5;  he[e4].prev = e6;  he[e4].used = 1;
        he[e5].tail = v2;  he[e5].face = f2_new;  he[e5].next = e6;  he[e5].prev = e4;  he[e5].used = 1;
        he[e6].tail = v3;  he[e6].face = f2_new;  he[e6].next = e4;  he[e6].prev = e5;  he[e6].used = 1;

        // 设置twin关系
        connectTwinEdges(e2, e6);
        connectTwinEdges(e1, outerEdgeV1V4);
        connectTwinEdges(e3, outerEdgeV3V1);
        connectTwinEdges(e4, outerEdgeV4V2);
        connectTwinEdges(e5, outerEdgeV2V3);

        // 更新顶点的出边引�????
        vertex[v1].edge = e1;
        vertex[v2].edge = e5;
        vertex[v3].edge = e6;
        vertex[v4].edge = e4;

        // 立即将新三角形添加到网格
        assignTriangleToGrid(f1_new);
        assignTriangleToGrid(f2_new);
        collected_triangles += 2; // 维护统计信息
        
        // 将新边加入栈中，以便�????�????
        if (stackSize + 4 <= MAX_STACK && flip_count < MAX_FLIPS) {
            edgeStack[stackSize] = e1;
            faceStack[stackSize] = f1_new;
            stackSize++;

            edgeStack[stackSize] = e3;
            faceStack[stackSize] = f1_new;
            stackSize++;

            edgeStack[stackSize] = e4;
            faceStack[stackSize] = f2_new;
            stackSize++;

            edgeStack[stackSize] = e5;
            faceStack[stackSize] = f2_new;
            stackSize++;
        }
    }
    
    // 更新�????大翻转次数统�????
    if (flip_count > max_flips_per_call) {
        max_flips_per_call = flip_count;
    }
    
}

// === 点插入函�???? ===
void insertSiteEnhanced(dnode_t p, didx_t containingFace = -1) {
    // 1. 转换点坐�????
    FixedPoint newPoint = getFixedPoint(p);
    static int insert_count = 0;
    insert_count++;

    // 2. 跟踪已处理点�????
    static int processed_actual_points = 0;
    
    // 3. 如果没有传入包含三角形，�????要自己查�????
    if (containingFace < 0) {
        // 第一个点：直接用超级三角�????
        if (processed_actual_points == 0) {
            containingFace = 0;
        } else {
            // 正常查找
            containingFace = locateTriangleSimple(newPoint, tri_num);
            // 如果找不到，直接返回
            if (containingFace < 0) return;
        }
    }
    
    // �????查是否是第一个点
    bool isFirstPoint = (processed_actual_points == 0);
    processed_actual_points++;
    didxh_t newEdgeIdx = global_edge_idx;

    // —�?? 以下为�?�用�????"拆分 containingFace"逻辑 —�??

    // 获取包含三角形的三个顶点 - �????化版�????
    didxh_t orig_e1 = face[containingFace].edge;  // 第一条边
    didxh_t orig_e2 = he[orig_e1].next;           // 第二条边
    didxh_t orig_e3 = he[orig_e2].next;           // 第三条边

    dnode_t v1 = he[orig_e1].tail;                // 第一个顶�????
    dnode_t v2 = he[orig_e2].tail;                // 第二个顶�????
    dnode_t v3 = he[orig_e3].tail;                // 第三个顶�????

    // 获取三角形顶点的坐标 - 修改为使用定点数版本
    FixedPoint p1 = getFixedPoint(v1);
    FixedPoint p2 = getFixedPoint(v2);
    FixedPoint p3 = getFixedPoint(v3);

    // 首先标记原三角形为未使用
    face[containingFace].used = 0;

    // 保存原外边对偶关系，供后续合法化使用
    didxh_t outerEdge12 = he[orig_e1].twin;  // v1-v2的对偶边
    didxh_t outerEdge23 = he[orig_e2].twin;  // v2-v3的对偶边
    didxh_t outerEdge31 = he[orig_e3].twin;  // v3-v1的对偶边

    // 使用全局函数释放半边
    freeHalfEdge(orig_e1);
    freeHalfEdge(orig_e2);
    freeHalfEdge(orig_e3);

    // 第一个三角形: p-v1-v2
    didx_t f1 = tri_num++;
    
    // 使用全局函数分配半tri_num++;    
    didxh_t e1 = allocateHalfEdge(newEdgeIdx);
    didxh_t e2 = allocateHalfEdge(newEdgeIdx);
    didxh_t e3 = allocateHalfEdge(newEdgeIdx);

    face[f1].edge = e1;
    face[f1].used = 1;

    he[e1].tail = p;
    he[e1].face = f1;
    he[e1].next = e2;
    he[e1].prev = e3;
    he[e1].used = 1;

    he[e2].tail = v1;
    he[e2].face = f1;
    he[e2].next = e3;
    he[e2].prev = e1;
    he[e2].used = 1;

    he[e3].tail = v2;
    he[e3].face = f1;
    he[e3].next = e1;
    he[e3].prev = e2;
    he[e3].used = 1;

    // 第二个三角形: p-v2-v3
    didx_t f2 = tri_num++; 
    
    // 使用全局函数分配半边
    didxh_t e4 = allocateHalfEdge(newEdgeIdx);
    didxh_t e5 = allocateHalfEdge(newEdgeIdx);
    didxh_t e6 = allocateHalfEdge(newEdgeIdx);

    face[f2].edge = e4;
    face[f2].used = 1;

    he[e4].tail = p;
    he[e4].face = f2;
    he[e4].next = e5;
    he[e4].prev = e6;
    he[e4].used = 1;

    he[e5].tail = v2;
    he[e5].face = f2;
    he[e5].next = e6;
    he[e5].prev = e4;
    he[e5].used = 1;

    he[e6].tail = v3;
    he[e6].face = f2;
    he[e6].next = e4;
    he[e6].prev = e5;
    he[e6].used = 1;

    // 第三个三角形: p-v3-v1
    didx_t f3 = tri_num++;
    
    // 使用全局函数分配半边
    didxh_t e7 = allocateHalfEdge(newEdgeIdx);
    didxh_t e8 = allocateHalfEdge(newEdgeIdx);
    didxh_t e9 = allocateHalfEdge(newEdgeIdx);

    face[f3].edge = e7;
    face[f3].used = 1;

    he[e7].tail = p;
    he[e7].face = f3;
    he[e7].next = e8;
    he[e7].prev = e9;
    he[e7].used = 1;

    he[e8].tail = v3;
    he[e8].face = f3;
    he[e8].next = e9;
    he[e8].prev = e7;
    he[e8].used = 1;

    he[e9].tail = v1;
    he[e9].face = f3;
    he[e9].next = e7;
    he[e9].prev = e8;
    he[e9].used = 1;

    // 连接内部半边对偶
    connectTwinEdges(e1, e9); // p-v1的对偶关�????
    connectTwinEdges(e4, e3); // p-v2的对偶关�????
    connectTwinEdges(e7, e6); // p-v3的对偶关�????

    // 连接外部半边对偶
    connectTwinEdges(e2, outerEdge12); // v1-v2边与外部连接
    connectTwinEdges(e5, outerEdge23); // v2-v3边与外部连接
    connectTwinEdges(e8, outerEdge31); // v3-v1边与外部连接

    // 将新三角形添加到网格
    assignTriangleToGrid(f1);
    assignTriangleToGrid(f2);
    assignTriangleToGrid(f3);

    // 仅对非第�????个点执行合法�????
    if (!isFirstPoint) {
        // 第一个三角形外边合法�????
        legalizeEdge(f1, e2,  newEdgeIdx);
        // 第二个三角形外边合法�????
        legalizeEdge(f2, e5,  newEdgeIdx);
        // 第三个三角形外边合法�????
        legalizeEdge(f3, e8,  newEdgeIdx);
    } else {
        // 第一个点不需要边翻转，只�????要记录一�????
        cout << "第一个点插入完成，跳过边翻转" << endl;
    }
    global_edge_idx = newEdgeIdx;
    
    // 插入完成后增加计数并�??查是否需要压�??
    inserted_point_count++;
    if (inserted_point_count >= COMPACTION_THRESHOLD) {
        compactTriangleIndices();
        inserted_point_count = 0;  // 重置计数�??
    }
}
       
     // 主函�???? - 移除了不�????要的edge_mat参数
void tri2d(hls::stream<axi_i_t>& x_in_stream,
    hls::stream<axi_o_t>& cl_out_stream) {
#pragma HLS INTERFACE axis port=x_in_stream
#pragma HLS INTERFACE axis port=cl_out_stream
#pragma HLS INTERFACE ap_ctrl_hs port=return

    tri_num = 0;
    is_cnt = 0;
    // 添加边索引初始化
    global_edge_idx = 6; // 重置为初始�??
    freeEdgeCount = 0;   // 清空复用�??
    inserted_point_count = 0; // 重置插入点计�??

    // 读取输入数据
    axi_i_t in_val; 
    axi_o_t val_out;
    unsigned int count = 0;

    // 读取�????有点（包括超大三角形顶点和普通点�????
    int point_index = 0;

    // 保证至少读取3个点（超大三角形顶点�????
    // 读取超大三角形顶�????
    for (int i = 0; i < 3; i++) {
        in_val = x_in_stream.read();
        
        // 提取位模�????
        uint32_t x_bits = in_val.data.range(31, 0);
        uint32_t y_bits = in_val.data.range(63, 32);
        
        // 转换回浮点�??
        float x_val, y_val;
        memcpy(&x_val, &x_bits, sizeof(float));
        memcpy(&y_val, &y_bits, sizeof(float));
        
        // 赋�?�给顶点
        vertex[i].x = x_val;
        vertex[i].y = y_val;
        vertex[i].used = 1;
    }

    // 读取普�?�点
    for (int i = 0; i < num_point; i++) {
        in_val = x_in_stream.read();
        
        // 提取位模�????
        uint32_t x_bits = in_val.data.range(31, 0);
        uint32_t y_bits = in_val.data.range(63, 32);
        
        // 转换回浮点�??
        float x_val, y_val;
        memcpy(&x_val, &x_bits, sizeof(float));
        memcpy(&y_val, &y_bits, sizeof(float));
        
        // 赋�?�给顶点
        vertex[i+3].x = x_val;
        vertex[i+3].y = y_val;
        vertex[i+3].used = 1;
    }

    // 初始化自适应网格 - 现在只计算边界和网格参数
    initEnhancedGrid(enhanced_grid, vertex, num_point+3);

    // ������ʼ������
    if (num_point >= 3) {
        // ��ʼ����һ�������� (0,1,2)
        face[0].edge = 0;
        face[0].used = 1;

        // �����������?
        he[0].tail = 0;  he[0].face = 0;  he[0].next = 1;  he[0].prev = 2;  he[0].used = 1;
        he[1].tail = 1;  he[1].face = 0;  he[1].next = 2;  he[1].prev = 0;  he[1].used = 1;
        he[2].tail = 2;  he[2].face = 0;  he[2].next = 0;  he[2].prev = 1;  he[2].used = 1;

        // �����ⲿ������
        he[3].tail = 1;  he[3].face = INVALID_INDEX;  he[3].next = 5;  he[3].prev = 4;  he[3].used = 1;
        he[4].tail = 2;  he[4].face = INVALID_INDEX;  he[4].next = 3;  he[4].prev = 5;  he[4].used = 1;
        he[5].tail = 0;  he[5].face = INVALID_INDEX;  he[5].next = 4;  he[5].prev = 3;  he[5].used = 1;

        // ����twin��ϵ
        he[0].twin = 3;  he[3].twin = 0;
        he[1].twin = 4;  he[4].twin = 1;
        he[2].twin = 5;  he[5].twin = 2;

        tri_num = 1;
    }

    // ����Hilbertֵ������
    PointWithZ points_with_z[MAX_VERTICES];
    int valid_point_count = 0;

    // �������е��Hilbertֵ
    for (int i = 3; i < num_vertice; i++) {
        #pragma HLS PIPELINE
        if (vertex[i].used) {
            FixedPoint p = getFixedPoint(i);

            hilbert_t hilbert_value = getHilbertCode(
                p.x, p.y,
                enhanced_grid.x_min, enhanced_grid.y_min,
                enhanced_grid.x_max, enhanced_grid.y_max,
                MAX_HILBERT_ORDER
            );
            
            points_with_z[valid_point_count].index = i;
            points_with_z[valid_point_count].z_value = hilbert_value;
            valid_point_count++;
        }
    }

    // ��������
    for (int i = 0; i < valid_point_count; i++) {
        PointWithZ temp = points_with_z[i];
        int j = i;
        
        while (j > 0 && points_with_z[j - 1].z_value > temp.z_value) {
            points_with_z[j] = points_with_z[j - 1];
            j--;
        }
        
        points_with_z[j] = temp;
    }

    // ��Hilbert˳������
    for (int i = 0; i < valid_point_count; i++) {
        int point_index = points_with_z[i].index;
        FixedPoint p = getFixedPoint(point_index);

        didx_t containing_face = locateTriangleSimple(p, tri_num);

        if (containing_face < 0) {
            continue;  // �����޷���λ�ĵ�
        }

        insertSiteEnhanced(point_index, containing_face);
    }

    // ����ѹ��
    if (inserted_point_count > 0) {
        compactTriangleIndices();
        inserted_point_count = 0;
    }

    // �ռ��������?
    ap_uint<MNPB * 3> con_list[MAX_FACES];
    int output_count = 0;

    for (int i = 0; i < MAX_NO_TRIANGLES; i++) {
        #pragma HLS PIPELINE

        if (face[i].used && face[i].edge != INVALID_HALFEDGE) {
            dnode_t v1, v2, v3;
            getTriangleVertices(i, v1, v2, v3);

            // ֻ����ǳ��������εĶ���?
            if (v1 > 2 && v2 > 2 && v3 > 2) {
                con_list[output_count] = 0;
                con_list[output_count].range(MNPB - 1, 0) = v1;
                con_list[output_count].range(2 * MNPB - 1, MNPB) = v2;
                con_list[output_count].range(3 * MNPB - 1, 2 * MNPB) = v3;
                output_count++;
            }
        }
    }

    // �����AXI��
    for (int i = 0; i < output_count; i++) {
        #pragma HLS PIPELINE II=1
        val_out.data = con_list[i];
        val_out.keep = 1;
        val_out.last = (i == output_count - 1) ? 1 : 0;
        val_out.user = 0;
        val_out.id = 0;
        val_out.dest = 0;
        val_out.strb = 1;

        cl_out_stream.write(val_out);
    }
}






