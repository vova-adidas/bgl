#pragma once

namespace bgl {
    
    /**
    * @class TriangleScanlineIterator
    * @brief Iterates over all scanlines of a triangle, covering every pixel
    *        that is touched by the triangle edges, including partial coverage.
    *
    * This class generates a series of scanlines (horizontal lines) that
    * intersect the triangle. For each scanline, it provides the range
    * of X coordinates that the triangle covers. It guarantees that all
    * pixels that the triangle touches are included.
    *
    * Features:
    *  - Full coverage of edge pixels.
    *  - Handles triangles of arbitrary orientation.
    *
    * Usage:
    * @code
    * TriangleScanlineIterator iter(ax, ay, bx, by, cx, cy);
    * while (iter.next()) {
    *     for (int x = iter.begin(); x < iter.end(); ++x) {
    *         int y = iter.line();
    *         // process pixel (x, y)
    *     }
    * }
    * @endcode
    */
    class TriangleScanlineIterator {
    private:

        inline constexpr static float EPSILON = 0.001f;

        inline constexpr static auto swap(auto& a, auto& b) {
            auto temp = a;
            a = b;
            b = temp;
        };

        inline constexpr static auto ceil(float x) {
            int i = static_cast<int>(x);
            return x > i ? i + 1 : i;
        };

        inline constexpr static auto floor(float x) {
            int i = static_cast<int>(x);
            return i;
        };

        inline constexpr static auto float_eq(auto x, auto y) {
            auto d = x - y;
            if (d < 0)
                d = -d;
            return d < EPSILON;
        };

        float
            left[6],
            right[6],
            delta_left[5]{},
            delta_right[5]{},
            cur_x0,
            cur_x1;

        int 
            n[5], 
            k = 0, 
            i = 0,
            j = 0,
            cur_y,
            cur_begin,
            cur_end,
            cur_line;

    public:

        inline constexpr TriangleScanlineIterator(const TriangleScanlineIterator&) = delete;

        inline constexpr TriangleScanlineIterator& operator=(const TriangleScanlineIterator&) = delete;

        inline constexpr TriangleScanlineIterator(float ax, float ay, float bx, float by, float cx, float cy) {
            
            if (by > cy) {
                swap(bx, cx);
                swap(by, cy);
            }

            if (ay > by) {
                swap(ax, bx);
                swap(ay, by);
            }

            if (by > cy) {
                swap(bx, cx);
                swap(by, cy);
            }

            const auto saturate = [](float x) { return x < 0 ? 0 : (x < 1 ? x : 1); };

            int min_y = floor(ay);
            int center_y = floor(by);
            int max_y = ceil(cy);
            
            if (static_cast<int>(ay) == static_cast<int>(cy)) {
                auto x0 = (ax < bx ? (ax < cx ? ax : cx) : (bx < cx ? bx : cx));
                auto x1 = (ax > bx ? (ax > cx ? ax : cx) : (bx > cx ? bx : cx));

                left[0] = x0;
                right[0] = x1;
                n[0] = 1;
                for (int i = 1; i < 5; ++i)
                    n[i] = 0;
                k++;
            }
            else
            {

                auto mid_y = by;
                auto left_x = bx;
                auto right_x = ax + (cx - ax) * (mid_y - ay) / (cy - ay);

                if (left_x > right_x)
                    swap(left_x, right_x);

                auto left_edge_correction = left_x - ax < 0 ? 1 : 0;
                auto right_edge_correction = right_x - ax > 0 ? 1 : 0;

                
                auto inv_by_ay = 1 / (by - ay);
                if (min_y < center_y) {
                    auto y = min_y;
                    auto t0 = (y + left_edge_correction - ay) * inv_by_ay;
                    auto t1 = (y + right_edge_correction - ay) * inv_by_ay;

                    left[k] = ax + (left_x - ax) * saturate(t0);
                    right[k] = (ax + (right_x - ax) * saturate(t1));
                    n[k] = 1;
                    k++;
                }

                if (center_y - (min_y + 1) > 0) {
                    
                    n[k] = center_y - (min_y + 1);

                    {
                        auto init = ax + (right_x - ax) * ((min_y + 1 - ay) * inv_by_ay);
                        auto end = ax + (right_x - ax) * ((center_y - ay) * inv_by_ay);

                        right[k] = init;
                        delta_right[k] = (end - init) / n[k];
                        right[k] += right_edge_correction * delta_right[k];
                    }

                    {
                        auto init = ax + (left_x - ax) * ((min_y + 1 - ay) * inv_by_ay);
                        auto end = ax + (left_x - ax) * ((center_y - ay) * inv_by_ay);

                        left[k] = init;
                        delta_left[k] = (end - init) / n[k];
                        left[k] += left_edge_correction * delta_left[k];
                    }
                    k++;
                }

                auto inv_cy_mid_y = 1 / (cy - mid_y);


                {
                    auto y = center_y;

                    float x0 = left_x;
                    auto edge = float_eq(by, ay) ? left_x : (ax + (left_x - ax) * saturate((y - ay) * inv_by_ay));
                    x0 = edge < x0 ? edge : x0;
                    edge = float_eq(cy, mid_y) ? left_x : (left_x + (cx - left_x) * saturate((y + 1 - mid_y) * inv_cy_mid_y));
                    x0 = edge < x0 ? edge : x0;

                    float x1 = right_x;
                    edge = float_eq(by, ay) ? right_x : (ax + (right_x - ax) * saturate((y - ay) * inv_by_ay));
                    x1 = edge > x1 ? edge : x1;
                    edge = float_eq(cy, mid_y) ? right_x : (right_x + (cx - right_x) * saturate((y + 1 - mid_y) * inv_cy_mid_y));
                    x1 = edge > x1 ? edge : x1;

                    n[k] = 1;
                    left[k] = x0;
                    right[k] = x1;
                    ++k;
                }

                left_edge_correction = cx - left_x < 0 ? 1 : 0;
                right_edge_correction = cx - right_x > 0 ? 1 : 0;

                if(max_y - 1 - (center_y + 1) > 0) {
                    n[k] = max_y - 1 - (center_y + 1);
                    {
                        auto init = right_x + (cx - right_x) * ((center_y + 1 - mid_y) * inv_cy_mid_y);
                        auto end = right_x + (cx - right_x) * ((max_y - 1 - mid_y) * inv_cy_mid_y);

                        right[k] = init;
                        delta_right[k] = (end - init) / n[k];
                        right[k] += right_edge_correction * delta_right[k];
                    }

                    {
                        auto init = left_x + (cx - left_x) * ((center_y + 1 - mid_y) * inv_cy_mid_y);
                        auto end = left_x + (cx - left_x) * ((max_y - 1 - mid_y) * inv_cy_mid_y);

                        left[k] = init;
                        delta_left[k] = (end - init) / n[k];
                        left[k] += left_edge_correction * delta_left[k];
                    }
                    ++k;
                }

                {
                    auto y = max_y - 1;
                    if (center_y < y) {
                        auto t0 = (y + left_edge_correction - mid_y) * inv_cy_mid_y;
                        auto t1 = (y + right_edge_correction - mid_y) * inv_cy_mid_y;

                        n[k] = 1;
                        left[k] = ((left_x + (cx - left_x) * saturate(t0)));
                        right[k] = ((right_x + (cx - right_x) * saturate(t1)));
                        ++k;
                    }
                }
            }

            cur_y = min_y;
            cur_x0 = left[j];
            cur_x1 = right[j];
        }

        inline constexpr bool next() {
            
            if (j >= k)
                return false;
       
            cur_begin = floor(cur_x0);
            cur_end = ceil(cur_x1);
            cur_line = cur_y;
            
            cur_x0 += delta_left[j];
            cur_x1 += delta_right[j];

            ++cur_y;
            ++i;


            if (i == n[j]) {
                ++j;
                i = 0;

                cur_x0 = left[j];
                cur_x1 = right[j];
            }
           

            return true;
        }

        inline constexpr int begin() const { return cur_begin; }

        inline constexpr int end() const { return cur_end; }

        inline constexpr int line() const { return cur_line; }
    };

/**
 * @brief Rasterizes a CW triangle in screen space using a top-left rule
 *
 * This function performs triangle rasterization in screen space with per-vertex attributes and
 * a user-defined SIMD warp type (`TWarp`).
 *
 * The algorithm iterates over the triangle tile-by-tile (tiles of size `TILE_SIZE`) and invokes
 * user-provided callbacks (`fragment`, `tile_lock`, `tile_unlock`) for fragment shading and
 * optional multithreaded synchronization.
 * 
 * The fragment callback is fully responsible for writing the results to the appropriate
 * render targets or framebuffers. The rasterizer itself does not perform any implicit pixel writes.
 * 
 * ---
 *
 * @tparam TWarp
 *     User-defined SIMD warp type
 * 
 * @tparam TILE_SIZE
 *     Tile size in pixels (e.g., 8 or 16).
 *
 * @tparam NUM_ATTRS
 *     Number of interpolated vertex attributes (e.g., 3 for UV + depth).
 *
 * ---
 *
 * @param ax, ay
 *     Screen-space coordinates of vertex A (in integer pixels).
 *
 * @param bx, by
 *     Screen-space coordinates of vertex B.
 *
 * @param cx, cy
 *     Screen-space coordinates of vertex C.
 *
 * @param a_attributes
 *     Array of `NUM_ATTRS` vertex attributes for vertex A (e.g., z, u, v, normals, etc.).
 *
 * @param b_attributes
 *     Array of `NUM_ATTRS` vertex attributes for vertex B.
 *
 * @param c_attributes
 *     Array of `NUM_ATTRS` vertex attributes for vertex C.
 *
 * @param fragment
 *     Fragment callback function invoked for each active tile region.
 *     Signature:
 *     @code
 *     void fragment(int tx, int ty,
 *                   int column, int row,
 *                   typename Warp::IntMask mask,
 *                   typename Warp::Float* attrs);
 *     @endcode
 *
 *     Parameters:
 *     - `tx`, `ty` � tile coordinates (in tiles, not pixels);
 *     - `column`, `row` � local pixel position inside the tile;
 *     - `mask` � bitmask of active pixels (edge function > 0);
 *     - `attrs` � interpolated vertex attributes for active pixels.
 *
 * @param tile_lock
 *     Function called before processing each tile.
 *     Used for synchronization in multithreaded rendering.
 *     May be a no-op lambda if synchronization is unnecessary.
 *
 * @param tile_unlock
 *     Function called after a tile has been processed.
 *     Used to release locks or finalize tile buffers.
 *
 * ---
 *
 * @note
 * - All vertex coordinates must be in integer screen-space pixels.
 * - The rasterization uses the **top-left coverage rule**.
 * - Can safely be used in a multithreaded context if `tile_lock` and `tile_unlock`
 *
 * ---
 *
 * @example
 * Example usage:
 * @code
 * topleft_triangle_scanline<MyWarp, 8, 3>(
 *     ax, ay, bx, by, cx, cy,
 *     a_attrs, b_attrs, c_attrs,
 *     [&](int tx, int ty, int column, int row, MyWarp::IntMask mask, MyWarp::Float* attrs) {
 *         // Fragment shading logic for pixels in this tile
 *     },
 *     [] (int tx, int ty) {}, // lock
 *     [] (int tx, int ty) {}  // unlock
 * );
 * @endcode
 */
    template<typename TWarp, int TILE_SIZE, int NUM_ATTRS>
    inline constexpr void topleft_triangle_scanline(
        int ax,
        int ay,
        int bx,
        int by,
        int cx,
        int cy,
        float (&a_attributes)[NUM_ATTRS],
        float (&b_attributes)[NUM_ATTRS],
        float (&c_attributes)[NUM_ATTRS],
        auto fragment,
        auto tile_lock,
        auto tile_unlock) {

        static_assert(TILE_SIZE >= 0, "Error: TILE_SIZE == 0");

        static_assert(TILE_SIZE >= TWarp::SIZE, "Error: TILE_SIZE < BATCH_SIZE");

        static_assert(TILE_SIZE % TWarp::SIZE == 0, "Error: TILE_SIZE % BATCH_SIZE != 0");

        const int SUBPIXEL_RES = 16;

        auto edge = [](auto ax, auto ay, auto bx, auto by, auto x, auto y) {
            return (by - ay) * x + (ax - bx) * y + (bx * ay - ax * by);
        };

        auto is_top_left = [](int dx, int dy) {
            return (dy > 0) || (dy == 0 && dx < 0);
        };

        int area = edge(ax, ay, bx, by, cx, cy);
        float inv_area = 1.f / area;

        if (area <= 0)
            return;

        int iax = static_cast<int>(SUBPIXEL_RES) * ax;
        int ibx = static_cast<int>(SUBPIXEL_RES) * bx;
        int icx = static_cast<int>(SUBPIXEL_RES) * cx;
        int iay = static_cast<int>(SUBPIXEL_RES) * ay;
        int iby = static_cast<int>(SUBPIXEL_RES) * by;
        int icy = static_cast<int>(SUBPIXEL_RES) * cy;

        int max_ix = iax > ibx ? (iax > icx ? iax : icx) : (ibx > icx ? ibx : icx);
        int max_iy = iay > iby ? (iay > icy ? iay : icy) : (iby > icy ? iby : icy);

        int scalar_start_edge[3]{
            edge(ibx, iby, icx, icy, SUBPIXEL_RES / 2, SUBPIXEL_RES / 2) + (is_top_left(icx - ibx, icy - iby) ? 0 : SUBPIXEL_RES),
            edge(icx, icy, iax, iay, SUBPIXEL_RES / 2, SUBPIXEL_RES / 2) + (is_top_left(iax - icx, iay - icy) ? 0 : SUBPIXEL_RES),
            edge(iax, iay, ibx, iby, SUBPIXEL_RES / 2, SUBPIXEL_RES / 2) + (is_top_left(ibx - iax, iby - iay) ? 0 : SUBPIXEL_RES)
        };

        int scalar_dx_edge[3]{
            (icy - iby) * SUBPIXEL_RES,
            (iay - icy) * SUBPIXEL_RES,
            (iby - iay) * SUBPIXEL_RES
        };

        int scalar_dy_edge[3]{
            (ibx - icx) * SUBPIXEL_RES,
            (icx - iax) * SUBPIXEL_RES,
            (iax - ibx) * SUBPIXEL_RES
        };

        alignas(TWarp::ALIGNMENT) int temp_idx[TWarp::SIZE];
        alignas(TWarp::ALIGNMENT) float temp_fidx[TWarp::SIZE];

        for (int i = 0; i < TWarp::SIZE; ++i) {
            temp_idx[i] = i;
            temp_fidx[i] = static_cast<float>(i);
        }

        auto idx = TWarp::load(temp_idx);
        auto fidx = TWarp::load(temp_fidx);

        typename TWarp::Int dx_edge[3];
        typename TWarp::Int dy_edge[3];
        typename TWarp::Int start_edge[3];
        typename TWarp::Int dx_wide_edge[3];

        for (int i = 0; i < 3; ++i) {
            dx_edge[i] = TWarp::scalar(scalar_dx_edge[i]);
            dy_edge[i] = TWarp::scalar(scalar_dy_edge[i]);
            start_edge[i] = TWarp::add(TWarp::mul(idx, dx_edge[i]), TWarp::scalar(scalar_start_edge[i]));
            dx_wide_edge[i] = TWarp::scalar(scalar_dx_edge[i] * TWarp::SIZE);
        }

        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_start[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_dx[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_dy[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_dx_wide[NUM_ATTRS];

        auto w0_dx = TWarp::scalar(scalar_dx_edge[0] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w1_dx = TWarp::scalar(scalar_dx_edge[1] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w2_dx = TWarp::scalar(scalar_dx_edge[2] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w0_dy = TWarp::scalar(scalar_dy_edge[0] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w1_dy = TWarp::scalar(scalar_dy_edge[1] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w2_dy = TWarp::scalar(scalar_dy_edge[2] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));

        auto w0_start = TWarp::fma(fidx, w0_dx, TWarp::scalar(inv_area * scalar_start_edge[0] / (SUBPIXEL_RES * SUBPIXEL_RES)));
        auto w1_start = TWarp::fma(fidx, w1_dx, TWarp::scalar(inv_area * scalar_start_edge[1] / (SUBPIXEL_RES * SUBPIXEL_RES)));
        auto w2_start = TWarp::fma(fidx, w2_dx, TWarp::scalar(inv_area * scalar_start_edge[2] / (SUBPIXEL_RES * SUBPIXEL_RES)));

        for (int i = 0; i < NUM_ATTRS; ++i) {
            auto attr_a = TWarp::scalar(a_attributes[i]);
            auto attr_b = TWarp::scalar(b_attributes[i]);
            auto attr_c = TWarp::scalar(c_attributes[i]);

            a_start[i] = TWarp::fma(attr_a, w0_start, TWarp::fma(attr_b, w1_start, TWarp::mul(attr_c, w2_start)));
            a_dx[i] = TWarp::fma(attr_a, w0_dx, TWarp::fma(attr_b, w1_dx, TWarp::mul(attr_c, w2_dx)));
            a_dy[i] = TWarp::fma(attr_a, w0_dy, TWarp::fma(attr_b, w1_dy, TWarp::mul(attr_c, w2_dy)));
            a_dx_wide[i] = TWarp::mul(a_dx[i], TWarp::scalar(TWarp::SIZE * 1.f));
        }

        TriangleScanlineIterator scanline(1.0f * ax / TILE_SIZE, 1.0f * ay / TILE_SIZE, 1.0f * bx / TILE_SIZE, 1.0f * by / TILE_SIZE, 1.0f * cx / TILE_SIZE, 1.0f * cy / TILE_SIZE);

        int line_batch_left[64];
        int line_batch_right[64];
        int line_batch_i = 0;

        bool has_lines = true;

        // Accumulating work in a batch for potential optimization
        while (has_lines) {

            int line_batch_start_ty = 0;
            line_batch_i = 0;

            while ((has_lines = scanline.next())) {
                if (line_batch_i == 0)
                    line_batch_start_ty = scanline.line();

                auto left = scanline.begin();
                auto right = scanline.end();
                auto rightmost_tile_start_x = (right - 1) * TILE_SIZE * SUBPIXEL_RES;

                //Skip empty tile
                if (max_ix <= rightmost_tile_start_x)
                    right--;

                line_batch_left[line_batch_i] = left;
                line_batch_right[line_batch_i] = right;


                line_batch_i++;
                if (line_batch_i == 64)
                    break;
            }
            int line_batch_size = line_batch_i;

            //Skip empty tile
            int topmost_tile_start_y = (line_batch_start_ty + line_batch_size - 1) * TILE_SIZE * SUBPIXEL_RES;
            if (max_iy <= topmost_tile_start_y)
                --line_batch_size;

            alignas(TWarp::ALIGNMENT) typename TWarp::Float attrs[NUM_ATTRS];
            alignas(TWarp::ALIGNMENT) typename TWarp::Int edges[3];

            for (int ty = line_batch_start_ty, line = 0; ty < line_batch_start_ty + line_batch_size; ++ty, ++line) {

                int column_step = TILE_SIZE / TWarp::SIZE;
         
                int tx_begin = line_batch_left[line];
                int tx_end = line_batch_right[line];
                int sx = tx_begin * TILE_SIZE;
                int sy = ty * TILE_SIZE;

                for (int i = 0; i < NUM_ATTRS; ++i)
                    attrs[i] = TWarp::fma(TWarp::scalar(sy * 1.f), a_dy[i], TWarp::fma(TWarp::scalar(sx * 1.f), a_dx[i], a_start[i]));

                for (int i = 0; i < 3; ++i)
                    edges[i] = TWarp::add(TWarp::mul(TWarp::scalar(sy), dy_edge[i]), TWarp::add(TWarp::mul(TWarp::scalar(sx), dx_edge[i]), start_edge[i]));

                for (int tx = tx_begin; tx < tx_end; ++tx) {
                    tile_lock(tx, ty);

                    for (int row = 0; row < TILE_SIZE; ++row) {

                        for (int column = 0; column < TILE_SIZE / TWarp::SIZE; ++column) {

                            auto edge_mask = TWarp::and_(TWarp::and_(
                                TWarp::greater(edges[0], TWarp::scalar(0)),
                                TWarp::greater(edges[1], TWarp::scalar(0))
                            ), TWarp::greater(edges[2], TWarp::scalar(0)));

                            fragment(tx, ty, column * TWarp::SIZE, row, edge_mask, attrs);

                            for (int i = 0; i < NUM_ATTRS; ++i)
                                attrs[i] = TWarp::add(attrs[i], a_dx_wide[i]);
                                
                            for (int i = 0; i < 3; ++i)
                                edges[i] = TWarp::add(edges[i], dx_wide_edge[i]);
                        }

                        for (int i = 0; i < NUM_ATTRS; ++i)
                            attrs[i] = TWarp::add(attrs[i], TWarp::fma(a_dx_wide[i], TWarp::scalar(-1.f * column_step), a_dy[i]));    
                            
                        for (int i = 0; i < 3; ++i)
                            edges[i] = TWarp::add(edges[i], TWarp::add(dy_edge[i], TWarp::mul(dx_wide_edge[i], TWarp::scalar(-column_step))));
                    }

                    for (int i = 0; i < NUM_ATTRS; ++i)
                        attrs[i] = TWarp::add(attrs[i], TWarp::mul(TWarp::fma(TWarp::scalar(-1.f), a_dy[i], a_dx[i]), TWarp::scalar(1.f * TILE_SIZE)));
                        
                    for (int i = 0; i < 3; ++i)
                        edges[i] = TWarp::add(edges[i], TWarp::mul(TWarp::add(dx_edge[i], TWarp::mul(dy_edge[i], TWarp::scalar(-1))), TWarp::scalar(TILE_SIZE)));

                    tile_unlock(tx, ty);
                }            
            }
        }
    }

    /**
 * @brief Rasterizes a CW triangle in screen space using a top-left rule
 *
 * This function performs triangle rasterization in screen space with per-vertex attributes and
 * a user-defined SIMD warp type (`TWarp`).
 *
 * The algorithm iterates over the triangle tile-by-tile (tiles of size `TILE_SIZE`) and invokes
 * user-provided callbacks (`fragment`, `tile_lock`, `tile_unlock`) for fragment shading and
 * optional multithreaded synchronization.
 *
 * The fragment callback is fully responsible for writing the results to the appropriate
 * render targets or framebuffers. The rasterizer itself does not perform any implicit pixel writes.
 *
 * ---
 *
 * @tparam TWarp
 *     User-defined SIMD warp type
 *
 * @tparam TILE_SIZE
 *     Tile size in pixels (e.g., 8 or 16).
 *
 * @tparam NUM_ATTRS
 *     Number of interpolated vertex attributes (e.g., 3 for UV + depth).
 *
 * @tparam EMPTY_TILE_CULLING
 *     Enables skipping of tiles fully outside the triangle.
 *     Can be disabled for small triangles to avoid extra culling overhead.
 *
 * ---
 *
 * @param ax, ay
 *     Screen-space coordinates of vertex A (in integer pixels).
 *
 * @param bx, by
 *     Screen-space coordinates of vertex B.
 *
 * @param cx, cy
 *     Screen-space coordinates of vertex C.
 *
 * @param a_attributes
 *     Array of `NUM_ATTRS` vertex attributes for vertex A (e.g., z, u, v, normals, etc.).
 *
 * @param b_attributes
 *     Array of `NUM_ATTRS` vertex attributes for vertex B.
 *
 * @param c_attributes
 *     Array of `NUM_ATTRS` vertex attributes for vertex C.
 *
 * @param fragment
 *     Fragment callback function invoked for each active tile region.
 *     Signature:
 *     @code
 *     void fragment(int tx, int ty,
 *                   int column, int row,
 *                   typename Warp::IntMask mask,
 *                   typename Warp::Float* attrs);
 *     @endcode
 *
 *     Parameters:
 *     - `tx`, `ty` � tile coordinates (in tiles, not pixels);
 *     - `column`, `row` � local pixel position inside the tile;
 *     - `mask` � bitmask of active pixels (edge function > 0);
 *     - `attrs` � interpolated vertex attributes for active pixels.
 *
 * @param tile_lock
 *     Function called before processing each tile.
 *     Used for synchronization in multithreaded rendering.
 *     May be a no-op lambda if synchronization is unnecessary.
 *
 * @param tile_unlock
 *     Function called after a tile has been processed.
 *     Used to release locks or finalize tile buffers.
 *
 * ---
 *
 * @note
 * - All vertex coordinates must be in integer screen-space pixels.
 * - The rasterization uses the **top-left coverage rule**.
 * - Can safely be used in a multithreaded context if `tile_lock` and `tile_unlock`
 *
 * ---
 *
 * @example
 * Example usage:
 * @code
 * topleft_triangle_scanline<MyWarp, 8, 3>(
 *     ax, ay, bx, by, cx, cy,
 *     a_attrs, b_attrs, c_attrs,
 *     [&](int tx, int ty, int column, int row, MyWarp::IntMask mask, MyWarp::Float* attrs) {
 *         // Fragment shading logic for pixels in this tile
 *     },
 *     [] (int tx, int ty) {}, // lock
 *     [] (int tx, int ty) {}  // unlock
 * );
 * @endcode
 */
    template<typename TWarp, int TILE_SIZE, int NUM_ATTRS, bool EMPTY_TILE_CULLING>
    inline constexpr void topleft_triangle_bounding_box(
        int ax,
        int ay,
        int bx,
        int by,
        int cx,
        int cy,
        float(&a_attributes)[NUM_ATTRS],
        float(&b_attributes)[NUM_ATTRS],
        float(&c_attributes)[NUM_ATTRS],
        auto fragment,
        auto tile_lock,
        auto tile_unlock) {

        static_assert(TILE_SIZE >= 0, "Error: TILE_SIZE == 0");

        static_assert(TILE_SIZE >= TWarp::SIZE, "Error: TILE_SIZE < BATCH_SIZE");

        static_assert(TILE_SIZE % TWarp::SIZE == 0, "Error: TILE_SIZE % BATCH_SIZE != 0");

        const int SUBPIXEL_RES = 16;

        auto edge = [](auto ax, auto ay, auto bx, auto by, auto x, auto y) {
            return (by - ay) * x + (ax - bx) * y + (bx * ay - ax * by);
            };

        auto is_top_left = [](int dx, int dy) {
            return (dy > 0) || (dy == 0 && dx < 0);
            };

        float area = edge(ax, ay, bx, by, cx, cy);
        float inv_area = 1 / area;

        if (area <= 0)
            return;

        int iax = static_cast<int>(SUBPIXEL_RES) * ax;
        int ibx = static_cast<int>(SUBPIXEL_RES) * bx;
        int icx = static_cast<int>(SUBPIXEL_RES) * cx;
        int iay = static_cast<int>(SUBPIXEL_RES) * ay;
        int iby = static_cast<int>(SUBPIXEL_RES) * by;
        int icy = static_cast<int>(SUBPIXEL_RES) * cy;

        int scalar_start_edge[3]{
            edge(ibx, iby, icx, icy, SUBPIXEL_RES / 2, SUBPIXEL_RES / 2) + (is_top_left(icx - ibx, icy - iby) ? 0 : SUBPIXEL_RES),
            edge(icx, icy, iax, iay, SUBPIXEL_RES / 2, SUBPIXEL_RES / 2) + (is_top_left(iax - icx, iay - icy) ? 0 : SUBPIXEL_RES),
            edge(iax, iay, ibx, iby, SUBPIXEL_RES / 2, SUBPIXEL_RES / 2) + (is_top_left(ibx - iax, iby - iay) ? 0 : SUBPIXEL_RES)
        };

        int scalar_dx_edge[3]{
            (icy - iby) * SUBPIXEL_RES,
            (iay - icy) * SUBPIXEL_RES,
            (iby - iay) * SUBPIXEL_RES
        };

        int scalar_dy_edge[3]{
            (ibx - icx) * SUBPIXEL_RES,
            (icx - iax) * SUBPIXEL_RES,
            (iax - ibx) * SUBPIXEL_RES
        };

        alignas(TWarp::ALIGNMENT) int temp_idx[TWarp::SIZE];
        alignas(TWarp::ALIGNMENT) float temp_fidx[TWarp::SIZE];

        for (int i = 0; i < TWarp::SIZE; ++i) {
            temp_idx[i] = i;
            temp_fidx[i] = i;
        }

        auto idx = TWarp::load(temp_idx);
        auto fidx = TWarp::load(temp_fidx);

        typename TWarp::Int dx_edge[3];
        typename TWarp::Int dy_edge[3];
        typename TWarp::Int start_edge[3];
        typename TWarp::Int dx_wide_edge[3];

        for (int i = 0; i < 3; ++i) {
            dx_edge[i] = TWarp::scalar(scalar_dx_edge[i]);
            dy_edge[i] = TWarp::scalar(scalar_dy_edge[i]);
            start_edge[i] = TWarp::add(TWarp::mul(idx, dx_edge[i]), TWarp::scalar(scalar_start_edge[i]));
            dx_wide_edge[i] = TWarp::scalar(scalar_dx_edge[i] * TWarp::SIZE);
        }

        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_start[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_dx[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_dy[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_dx_wide[NUM_ATTRS];

        auto w0_dx = TWarp::scalar(scalar_dx_edge[0] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w1_dx = TWarp::scalar(scalar_dx_edge[1] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w2_dx = TWarp::scalar(scalar_dx_edge[2] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w0_dy = TWarp::scalar(scalar_dy_edge[0] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w1_dy = TWarp::scalar(scalar_dy_edge[1] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));
        auto w2_dy = TWarp::scalar(scalar_dy_edge[2] * inv_area / (SUBPIXEL_RES * SUBPIXEL_RES));

        auto w0_start = TWarp::fma(fidx, w0_dx, TWarp::scalar(inv_area * scalar_start_edge[0] / (SUBPIXEL_RES * SUBPIXEL_RES)));
        auto w1_start = TWarp::fma(fidx, w1_dx, TWarp::scalar(inv_area * scalar_start_edge[1] / (SUBPIXEL_RES * SUBPIXEL_RES)));
        auto w2_start = TWarp::fma(fidx, w2_dx, TWarp::scalar(inv_area * scalar_start_edge[2] / (SUBPIXEL_RES * SUBPIXEL_RES)));

        for (int i = 0; i < NUM_ATTRS; ++i) {
            auto attr_a = TWarp::scalar(a_attributes[i]);
            auto attr_b = TWarp::scalar(b_attributes[i]);
            auto attr_c = TWarp::scalar(c_attributes[i]);

            a_start[i] = TWarp::fma(attr_a, w0_start, TWarp::fma(attr_b, w1_start, TWarp::mul(attr_c, w2_start)));
            a_dx[i] = TWarp::fma(attr_a, w0_dx, TWarp::fma(attr_b, w1_dx, TWarp::mul(attr_c, w2_dx)));
            a_dy[i] = TWarp::fma(attr_a, w0_dy, TWarp::fma(attr_b, w1_dy, TWarp::mul(attr_c, w2_dy)));
            a_dx_wide[i] = TWarp::mul(a_dx[i], TWarp::scalar(TWarp::SIZE * 1.f));
        }

        int max_ix = ax > bx ? (ax > cx ? ax : cx) : (bx > cx ? bx : cx);
        int min_ix = ax < bx ? (ax < cx ? ax : cx) : (bx < cx ? bx : cx);
        int max_iy = ay > by ? (ay > cy ? ay : cy) : (by > cy ? by : cy);
        int min_iy = ay < by ? (ay < cy ? ay : cy) : (by < cy ? by : cy);

        auto ceil = [](float x) {
            int i = static_cast<int>(x);
            return x > i ? i + 1 : i;
        };

        auto floor = [](float x) {
            int i = static_cast<int>(x);
            return i;
        };

        int tile_start_x = floor(1.f * min_ix / TILE_SIZE);
        int tile_start_y = floor(1.f * min_iy / TILE_SIZE);
        int tile_end_x = ceil(1.f * max_ix / TILE_SIZE);
        int tile_end_y = ceil(1.f * max_iy / TILE_SIZE);

        if ((tile_end_x - 1) * TILE_SIZE >= max_ix) --tile_end_x;
        if ((tile_end_y - 1) * TILE_SIZE >= max_iy) --tile_end_y;

        alignas(TWarp::ALIGNMENT) typename TWarp::Float attrs[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Int edges[3];

        for (int ty = tile_start_y; ty < tile_end_y; ++ty) {

            int tile_start_x_culled = tile_start_x;
            int tile_end_x_culled = tile_end_x;

            if constexpr (EMPTY_TILE_CULLING) {
                const auto max3 = [&](int a, int b, int c) { return a > b ? (a > c ? a : c) : (b > c ? b : c); };
                const auto e0 = [&](int x, int y) { return scalar_start_edge[0] + scalar_dx_edge[0] * x + scalar_dy_edge[0] * y; };
                const auto e1 = [&](int x, int y) { return scalar_start_edge[1] + scalar_dx_edge[1] * x + scalar_dy_edge[1] * y; };
                const auto e2 = [&](int x, int y) { return scalar_start_edge[2] + scalar_dx_edge[2] * x + scalar_dy_edge[2] * y; };

                const auto is_tile_outside = [&](int tx, int ty) {
                    return
                        max3({ e0(tx * TILE_SIZE, ty * TILE_SIZE),
                        e0(tx * TILE_SIZE + TILE_SIZE - 1, ty * TILE_SIZE),
                        e0(tx * TILE_SIZE, ty * TILE_SIZE + TILE_SIZE - 1),
                        e0(tx * TILE_SIZE + TILE_SIZE - 1, ty * TILE_SIZE + TILE_SIZE - 1) }) < 0
                        ||
                        max3({ e1(tx * TILE_SIZE, ty * TILE_SIZE),
                        e1(tx * TILE_SIZE + TILE_SIZE - 1, ty * TILE_SIZE),
                        e1(tx * TILE_SIZE, ty * TILE_SIZE + TILE_SIZE - 1),
                        e1(tx * TILE_SIZE + TILE_SIZE - 1, ty * TILE_SIZE + TILE_SIZE - 1) }) < 0
                        ||
                        max3({ e2(tx * TILE_SIZE, ty * TILE_SIZE),
                        e2(tx * TILE_SIZE + TILE_SIZE - 1, ty * TILE_SIZE),
                        e2(tx * TILE_SIZE, ty * TILE_SIZE + TILE_SIZE - 1),
                        e2(tx * TILE_SIZE + TILE_SIZE - 1, ty * TILE_SIZE + TILE_SIZE - 1) }) < 0;
                    };

                for (int tx = tile_start_x; tx < tile_end_x; ++tx) {
                    if (is_tile_outside(tx, ty))
                        tile_start_x_culled++;
                    else
                        break;
                }

                for (int tx = tile_end_x - 1; tx > tile_start_x_culled; --tx) {

                    if (is_tile_outside(tx, ty))
                        tile_end_x_culled--;
                    else
                        break;
                }
            }

            int column_step = TILE_SIZE / TWarp::SIZE;

            int sy = ty * TILE_SIZE;
            int sx = tile_start_x_culled * TILE_SIZE;

            for (int i = 0; i < NUM_ATTRS; ++i)
                attrs[i] = TWarp::fma(TWarp::scalar(sy * 1.f), a_dy[i], TWarp::fma(TWarp::scalar(sx * 1.f), a_dx[i], a_start[i]));

            for (int i = 0; i < 3; ++i)
                edges[i] = TWarp::add(TWarp::mul(TWarp::scalar(sy), dy_edge[i]), TWarp::add(TWarp::mul(TWarp::scalar(sx), dx_edge[i]), start_edge[i]));

            for (int tx = tile_start_x_culled; tx < tile_end_x_culled; ++tx) {

                tile_lock(tx, ty);

                for (int row = 0; row < TILE_SIZE; ++row) {

                    for (int column = 0; column < TILE_SIZE / TWarp::SIZE; ++column) {

                        auto edge_mask = TWarp::and_(TWarp::and_(
                            TWarp::greater(edges[0], TWarp::scalar(0)),
                            TWarp::greater(edges[1], TWarp::scalar(0))
                        ), TWarp::greater(edges[2], TWarp::scalar(0)));

                        fragment(tx, ty, column * TWarp::SIZE, row, edge_mask, attrs);

                        for (int i = 0; i < NUM_ATTRS; ++i)
                            attrs[i] = TWarp::add(attrs[i], a_dx_wide[i]);

                        for (int i = 0; i < 3; ++i)
                            edges[i] = TWarp::add(edges[i], dx_wide_edge[i]);
                    }

                    for (int i = 0; i < NUM_ATTRS; ++i)
                        attrs[i] = TWarp::add(attrs[i], TWarp::fma(a_dx_wide[i], TWarp::scalar(-1.f * column_step), a_dy[i]));

                    for (int i = 0; i < 3; ++i)
                        edges[i] = TWarp::add(edges[i], TWarp::add(dy_edge[i], TWarp::mul(dx_wide_edge[i], TWarp::scalar(-column_step))));
                }

                for (int i = 0; i < NUM_ATTRS; ++i)
                    attrs[i] = TWarp::add(attrs[i], TWarp::mul(TWarp::fma(TWarp::scalar(-1.f), a_dy[i], a_dx[i]), TWarp::scalar(1.f * TILE_SIZE)));

                for (int i = 0; i < 3; ++i)
                    edges[i] = TWarp::add(edges[i], TWarp::mul(TWarp::add(dx_edge[i], TWarp::mul(dy_edge[i], TWarp::scalar(-1))), TWarp::scalar(TILE_SIZE)));

                tile_unlock(tx, ty);
            }
        }
    }
}