/**
 * @file tri2d_tb.cpp
 * @brief Simplified 2D Delaunay Triangulation (Tri2D) IP core testbench
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
    {fixed_t(-1.2), fixed_t(-0.1)},
    {fixed_t(2.2), fixed_t(-0.1)},
    {fixed_t(0.5), fixed_t(1.6)}
};

int main(int argc, char** argv) {
    // Variables
    FILE* dma_file = nullptr;
    fixed_t x_min = fixed_t(1e30f), x_max = fixed_t(-1e30f);
    fixed_t y_min = fixed_t(1e30f), y_max = fixed_t(-1e30f);
    fixed_t z_min = fixed_t(1e30f), z_max = fixed_t(-1e30f);
    int cnt = 0;

    string pcd_buffer;
    const string end_of_header = "DATA ascii";
    fixed_t x_tmp = fixed_t(0), y_tmp = fixed_t(0), z_tmp = fixed_t(0);
    fixed_t x_coor[MAX_NO_POINTS];
    fixed_t y_coor[MAX_NO_POINTS];
    fixed_t z_coor[MAX_NO_POINTS];

    cout << "=== Delaunay Triangulation Testbench ===" << endl;

    // Initialize DCEL data structures
    for (int i = 0; i < MAX_NO_TRIANGLES; i++) {
        face[i].edge = INVALID_HALFEDGE;
        face[i].prev = INVALID_INDEX;
        face[i].next = INVALID_INDEX;
        face[i].used = 0;
    }

    for (int i = 0; i < MAX_NO_POINTS; i++) {
        vertex[i].x = fixed_t(0);
        vertex[i].y = fixed_t(0);
        vertex[i].hilbert_code = 0;
        vertex[i].edge = INVALID_HALFEDGE;
        vertex[i].prev = INVALID_INDEX;
        vertex[i].next = INVALID_INDEX;
        vertex[i].used = 0;
        vertex[i].processed = false;
    }

    for (int i = 0; i < MAX_NO_HALFEDGES; i++) {
        he[i].tail = 0;
        he[i].twin = INVALID_HALFEDGE;
        he[i].prev = INVALID_HALFEDGE;
        he[i].next = INVALID_HALFEDGE;
        he[i].face = INVALID_INDEX;
        he[i].used = 0;
    }

    // Open input file
    string inputFilePath = "D:/src1/bunny_sparse_fps_5000.pcd";
    if (argc > 1) {
        inputFilePath = argv[1];
    }
    cout << "Opening input file: " << inputFilePath << endl;

    ifstream file(inputFilePath);
    if (!file.is_open()) {
        cerr << "Error: Failed to open file: " << inputFilePath << endl;
        return -1;
    }

    // Skip PCD file header
    do {
        getline(file, pcd_buffer);
        if (file.fail()) {
            cerr << "Error: Failed to read PCD header." << endl;
            return -2;
        }
        if (pcd_buffer.find(end_of_header) != string::npos) {
            break;
        }
        if (pcd_buffer.empty()) continue;
    } while (true);

    // Parse PCD data section
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

    // Set super-triangle vertices
    for (int i = 0; i < 3; ++i) {
        vertex[i].x = SUPERTRIANGLE_BOUNDS[i][0];
        vertex[i].y = SUPERTRIANGLE_BOUNDS[i][1];
        vertex[i].used = 1;
    }

    // Calculate normalization parameters (excluding super-triangle)
    cout << "Normalizing point coordinates..." << endl;
    fixed_t range_x = x_max - x_min;
    fixed_t range_y = y_max - y_min;

    // Avoid division by zero
    if (range_x.to_float() == 0.0f || range_y.to_float() == 0.0f) {
        cerr << "Error: Invalid range for normalization." << endl;
        return -5;
    }

    // Store normalized points to vertex array
    for (int i = 0; i < cnt; ++i) {
        // Normalize coordinates to [0,1] range
        fixed_t norm_x = (x_coor[i] - x_min) / range_x;
        fixed_t norm_y = (y_coor[i] - y_min) / range_y;

        vertex[i + 3].x = norm_x;
        vertex[i + 3].y = norm_y;
        vertex[i + 3].used = 1;
    }

    // Prepare input/output streams
    hls::stream<axi_i_t> x_in_stream("x_in_stream");
    hls::stream<axi_o_t> cl_out_stream("cl_out_stream");
    axi_i_t in_val;

    // Open DMA file - use command line argument or default path
    string dmaFilePath = "dma_data.bin";
    if (argc > 2) {
        dmaFilePath = argv[2];
    }

    dma_file = fopen(dmaFilePath.c_str(), "wb+");
    if (dma_file == NULL) {
        cerr << "Warning: Unable to open DMA output file. Continuing without DMA file..." << endl;
    }

    // Write all point data to stream and DMA file
    cout << "Writing " << (cnt + 3) << " vertices to input stream..." << endl;
    for (int i = 0; i < cnt + 3; i++) {
        in_val.last = 0;
        float x_val = vertex[i].x.to_float();
        float y_val = vertex[i].y.to_float();

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

        // Write to stream
        x_in_stream.write(in_val);

        // If DMA file successfully opened, also write to file
        if (dma_file != NULL) {
            fwrite(&in_val, sizeof(in_val), 1, dma_file);
        }
    }

    // Close DMA file if successfully opened
    if (dma_file != NULL) {
        fclose(dma_file);
    }

    cout << "\n=== Starting Triangulation Process ===" << endl;

    // Execute triangulation
    clock_t begin = clock();
    tri2d(x_in_stream, cl_out_stream);
    clock_t end = clock();

    double execution_time_ms = ((double)(end - begin) / CLOCKS_PER_SEC) * 1000;
    cout << "Execution time: " << execution_time_ms << " ms" << endl;

    // Read triangle data from output stream
    cout << "\n=== Reading Triangulation Results ===" << endl;
    axi_o_t val_out;
    ap_uint<MNPB * 3> con_list[MAX_FACES];
    int output_count = 0;

    // Open output file
    FILE* result_file = fopen("triangulation_result.txt", "w");
    if (result_file == NULL) {
        cerr << "Error: Unable to open output file for writing." << endl;
        return -6;
    }

    // Write header information
    fprintf(result_file, "# Delaunay Triangulation Result\n");
    fprintf(result_file, "# Format: v1 v2 v3 (vertex indices, 0-based)\n");

    // Read triangles with fixed number of iterations instead of infinite wait
    cout << "Reading triangles from output stream..." << endl;

    for (int i = 0; i < MAX_FACES; i++) {
        if (!cl_out_stream.empty()) {
            val_out = cl_out_stream.read();
            con_list[output_count] = val_out.data;

            // Extract vertex indices
            ap_uint<MNPB> v1 = con_list[output_count].range(MNPB - 1, 0);
            ap_uint<MNPB> v2 = con_list[output_count].range(2 * MNPB - 1, MNPB);
            ap_uint<MNPB> v3 = con_list[output_count].range(3 * MNPB - 1, 2 * MNPB);

            // Write data to output file
            fprintf(result_file, "%u %u %u\n",
                (unsigned int)v1, (unsigned int)v2, (unsigned int)v3);

            output_count++;

            // If it's the last triangle, exit the loop
            if (val_out.last) {
                break;
            }
        }
        else if (i > 0 && i % 1000 == 0) {
            // Avoid infinite loop
            if (cl_out_stream.empty() && output_count > 0) {
                break;
            }
        }
    }

    // Close output file
    fclose(result_file);

    cout << "\n=== Triangulation Summary ===" << endl;
    cout << "Total triangles generated: " << output_count << endl;
    cout << "Triangulation completed" << endl;

    return 0;
}



