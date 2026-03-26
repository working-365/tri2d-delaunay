
/**
 * @file test.cpp
 * @brief KV260 Optimized Delaunay Triangulation Testbench
 *
 * Changes: all vertex[i].x -> vtx_x[i], face[i].xxx -> face_topo[i].xxx, etc.
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <algorithm>
#include <ap_fixed.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include "tri2d.hpp"

using namespace std;

// Super triangle vertices coordinates
const fixed_t SUPERTRIANGLE_BOUNDS[3][2] = {
    {fixed_t(-2.5), fixed_t(-1.0)},
    {fixed_t(2.5), fixed_t(-1.0)},
    {fixed_t(0.0), fixed_t(2.5)}
};

int main(int argc, char** argv) {
    FILE* dma_file = nullptr;
    fixed_t x_min = fixed_t(1e30f), x_max = fixed_t(-1e30f);
    fixed_t y_min = fixed_t(1e30f), y_max = fixed_t(-1e30f);
    fixed_t z_min = fixed_t(1e30f), z_max = fixed_t(-1e30f);
    int cnt = 0;

    string pcd_buffer;
    const string end_of_header = "DATA ascii";
    fixed_t x_tmp, y_tmp, z_tmp;
    fixed_t x_coor[MAX_NO_POINTS];
    fixed_t y_coor[MAX_NO_POINTS];
    fixed_t z_coor[MAX_NO_POINTS];

    cout << "=== Delaunay Triangulation Testbench (KV260 Optimized) ===" << endl;

    // ============================================
    // Initialize DCEL data structures (using SoA)
    // ============================================
    for (int i = 0; i < MAX_NO_TRIANGLES; i++) {
        face_topo[i].edge = INVALID_HALFEDGE;
        face_topo[i].prev = INVALID_INDEX;
        face_topo[i].next = INVALID_INDEX;
        face_topo[i].used = 0;
    }

    for (int i = 0; i < MAX_NO_POINTS; i++) {
        vtx_x[i] = fixed_t(0);
        vtx_y[i] = fixed_t(0);
        vtx_hilbert[i] = 0;
        vtx_edge[i] = INVALID_HALFEDGE;
        vtx_prev[i] = INVALID_INDEX;
        vtx_next[i] = INVALID_INDEX;
        vtx_used[i] = 0;
        vtx_processed[i] = 0;
    }

    for (int i = 0; i < MAX_NO_HALFEDGES; i++) {
        he_tail[i] = 0;
        he_twin[i] = INVALID_HALFEDGE;
        he_prev[i] = INVALID_HALFEDGE;
        he_next[i] = INVALID_HALFEDGE;
        he_face[i] = INVALID_INDEX;
        he_used[i] = 0;
    }

    // ============================================
    // Read PCD file
    // ============================================
    string inputFilePath = "D:/src1/bunny/bunny_1024_uniform.pcd";
    if (argc > 1) {
        inputFilePath = argv[1];
    }
    cout << "Opening input file: " << inputFilePath << endl;

    ifstream file(inputFilePath);
    if (!file.is_open()) {
        cerr << "Error: Failed to open file: " << inputFilePath << endl;
        return -1;
    }

    // Skip PCD header
    do {
        getline(file, pcd_buffer);
        if (file.fail()) {
            cerr << "Error: Failed to read PCD header." << endl;
            return -2;
        }
        if (pcd_buffer.find(end_of_header) != string::npos) break;
        if (pcd_buffer.empty()) continue;
    } while (true);

    // Parse PCD data
    cout << "Parsing PCD data points..." << endl;
    while (getline(file, pcd_buffer)) {
        stringstream slicer(pcd_buffer);
        float x_float, y_float, z_float;
        if (!(slicer >> x_float >> y_float >> z_float)) {
            cerr << "Invalid input format in line: " << pcd_buffer << endl;
            continue;
        }

        x_tmp = fixed_t(x_float);
        y_tmp = fixed_t(y_float);
        z_tmp = fixed_t(z_float);

        x_coor[cnt] = x_tmp;
        y_coor[cnt] = y_tmp;
        z_coor[cnt] = z_tmp;

        x_min = min(x_min, x_tmp);
        y_min = min(y_min, y_tmp);
        z_min = min(z_min, z_tmp);
        x_max = max(x_max, x_tmp);
        y_max = max(y_max, y_tmp);
        z_max = max(z_max, z_tmp);

        ++cnt;
    }
    file.close();
    cout << "Parsing completed. Points count: " << cnt << endl;

    // ============================================
    // Set super-triangle vertices (using SoA)
    // ============================================
    for (int i = 0; i < 3; ++i) {
        vtx_x[i] = SUPERTRIANGLE_BOUNDS[i][0];
        vtx_y[i] = SUPERTRIANGLE_BOUNDS[i][1];
        vtx_used[i] = 1;
    }

    // ============================================
    // Normalize point coordinates
    // ============================================
    cout << "Normalizing point coordinates..." << endl;
    fixed_t range_x = x_max - x_min;
    fixed_t range_y = y_max - y_min;

    if (range_x.to_float() == 0.0f || range_y.to_float() == 0.0f) {
        cerr << "Error: Invalid range for normalization." << endl;
        return -5;
    }

    for (int i = 0; i < cnt; ++i) {
        fixed_t norm_x = (x_coor[i] - x_min) / range_x;
        fixed_t norm_y = (y_coor[i] - y_min) / range_y;
        vtx_x[i + 3] = norm_x;
        vtx_y[i + 3] = norm_y;
        vtx_used[i + 3] = 1;
    }

    // ============================================
    // Global Hilbert sort + stride-jump insertion order
    // ============================================
    cout << "Computing global Hilbert sorting with stride sampling..." << endl;
    fixed_t h_x_min = fixed_t(0), h_y_min = fixed_t(0);
    fixed_t h_x_max = fixed_t(1), h_y_max = fixed_t(1);
    const int H_ORDER = 8;

    // Step 1: compute Hilbert codes globally
    int       global_idx[MAX_NO_POINTS];
    hilbert_t global_code[MAX_NO_POINTS];
    for (int i = 0; i < cnt; i++) {
        global_idx[i] = i;
        global_code[i] = getHilbertCode(
            vtx_x[i + 3], vtx_y[i + 3],
            h_x_min, h_y_min, h_x_max, h_y_max, H_ORDER);
    }

    // Step 2: sort indices by Hilbert code
    std::sort(global_idx, global_idx + cnt, [&](int a, int b) {
        return global_code[a] < global_code[b];
        });

    // Step 3: stride-jump reordering
    // First spread points sparsely across the domain, then fill gaps in subsequent
    // rounds to avoid local grid overflow
    int stride = (int)sqrt((double)cnt);  // ~70 for 5000 points
    fixed_t sorted_x[MAX_NO_POINTS];
    fixed_t sorted_y[MAX_NO_POINTS];
    int out_pos = 0;
    for (int start = 0; start < stride; start++) {
        for (int i = start; i < cnt; i += stride) {
            sorted_x[out_pos] = vtx_x[global_idx[i] + 3];
            sorted_y[out_pos] = vtx_y[global_idx[i] + 3];
            out_pos++;
        }
    }

    // Step 4: write back to vertex SoA
    for (int i = 0; i < cnt; i++) {
        vtx_x[i + 3] = sorted_x[i];
        vtx_y[i + 3] = sorted_y[i];
        vtx_used[i + 3] = 1;
    }
    cout << "Stride sampling completed. "
        << "cnt=" << cnt << " stride=" << stride
        << " rounds=" << stride << endl;

    // ============================================
    // Prepare input/output streams
    // ============================================
    hls::stream<axi_i_t> x_in_stream("x_in_stream");
    hls::stream<axi_o_t> cl_out_stream("cl_out_stream");
    axi_i_t in_val;

    string dmaFilePath = "dma_data.bin";
    if (argc > 2) {
        dmaFilePath = argv[2];
    }

    dma_file = fopen(dmaFilePath.c_str(), "wb+");
    if (dma_file == NULL) {
        cerr << "Warning: Unable to open DMA output file. Continuing..." << endl;
    }

    // Write all points to the input stream
    cout << "Writing " << (cnt + 3) << " vertices to input stream..." << endl;
    for (int i = 0; i < cnt + 3; i++) {
        in_val.last = 0;
        float x_val = vtx_x[i].to_float();
        float y_val = vtx_y[i].to_float();

        uint32_t x_bits, y_bits;
        memcpy(&x_bits, &x_val, sizeof(float));
        memcpy(&y_bits, &y_val, sizeof(float));

        in_val.data.range(31, 0) = x_bits;
        in_val.data.range(63, 32) = y_bits;
        in_val.keep = 1;
        in_val.strb = 1;
        in_val.user = 0;
        in_val.last = (i == cnt + 3 - 1);
        in_val.id = 0;
        in_val.dest = 0;

        x_in_stream.write(in_val);

        if (dma_file != NULL) {
            fwrite(&in_val, sizeof(in_val), 1, dma_file);
        }
    }

    if (dma_file != NULL) fclose(dma_file);

    cout << "\n=== Starting Triangulation Process ===" << endl;

    clock_t begin = clock();
volatile int debug_count = 0;
tri2d(x_in_stream, cl_out_stream, debug_count);
printf("debug_count = %d\n", (int)debug_count);
    clock_t end = clock();

    double execution_time_ms = ((double)(end - begin) / CLOCKS_PER_SEC) * 1000;
    cout << "Execution time: " << execution_time_ms << " ms" << endl;

    // ============================================
    // Read triangulation results
    // ============================================
    cout << "\n=== Reading Triangulation Results ===" << endl;
    axi_o_t val_out;
    ap_uint<MNPB * 3> con_list[MAX_FACES];
    int output_count = 0;

    FILE* result_file = fopen("triangulation_result.txt", "w");
    if (result_file == NULL) {
        cerr << "Error: Unable to open output file for writing." << endl;
        return -6;
    }

    fprintf(result_file, "# Delaunay Triangulation Result\n");
    fprintf(result_file, "# VERTICES %d\n", cnt);

    // Write vertex coordinates (using SoA)
    for (int i = 0; i < cnt; i++) {
        fprintf(result_file, "v %f %f\n",
            vtx_x[i + 3].to_float(), vtx_y[i + 3].to_float());
    }

    fprintf(result_file, "# TRIANGLES\n");

    cout << "Reading triangles from output stream..." << endl;
    for (int i = 0; i < MAX_FACES; i++) {
        if (!cl_out_stream.empty()) {
            val_out = cl_out_stream.read();
            con_list[output_count] = val_out.data;

            ap_uint<MNPB> v1 = con_list[output_count].range(MNPB - 1, 0);
            ap_uint<MNPB> v2 = con_list[output_count].range(2 * MNPB - 1, MNPB);
            ap_uint<MNPB> v3 = con_list[output_count].range(3 * MNPB - 1, 2 * MNPB);

            int ov1 = (int)v1 - 3;
            int ov2 = (int)v2 - 3;
            int ov3 = (int)v3 - 3;
            fprintf(result_file, "t %d %d %d\n", ov1, ov2, ov3);

            output_count++;
            if (val_out.last) break;
        }
        else if (i > 0 && i % 1000 == 0) {
            if (cl_out_stream.empty() && output_count > 0) break;
        }
    }

    fclose(result_file);

    cout << "\n=== Triangulation Summary ===" << endl;
    cout << "Total triangles generated: " << output_count << endl;
    cout << "Triangulation completed" << endl;

    return 0;
}
