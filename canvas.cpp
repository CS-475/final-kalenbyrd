#include "include/GCanvas.h"
#include "include/GBitmap.h"
#include "include/GColor.h"
#include "include/GRect.h"
#include "include/GPaint.h"
#include "include/GShader.h"
#include "include/GPath.h"
#include <ostream>
#include <iostream>
#include <stack>
#include <vector>
#include <algorithm>


static GPixel premultiplyColor(const GColor& color) {
    if (color.a == 0.0f) {
        return GPixel_PackARGB(0, 0, 0, 0);
    }

    auto clampAndRound = [](float value) -> unsigned {
        return static_cast<unsigned>(std::round(GPinToUnit(value) * 255));
    };

    unsigned a = clampAndRound(color.a);
    unsigned r = clampAndRound(color.r * color.a);
    unsigned g = clampAndRound(color.g * color.a);
    unsigned b = clampAndRound(color.b * color.a);

    r = std::min(r, a);
    g = std::min(g, a);
    b = std::min(b, a);

    return GPixel_PackARGB(a, r, g, b);
}


static inline unsigned div255(unsigned x) {
    return (x + 128) * 257 >> 16;
}

static GPixel blendPixel(const GPixel& src, const GPixel& dst, GBlendMode mode) {
    unsigned src_a = GPixel_GetA(src);
    unsigned dst_a = GPixel_GetA(dst);

    unsigned inv_src_a = 255 - src_a;
    unsigned inv_dst_a = 255 - dst_a;

    unsigned src_r = GPixel_GetR(src);
    unsigned src_g = GPixel_GetG(src);
    unsigned src_b = GPixel_GetB(src);

    unsigned dst_r = GPixel_GetR(dst);
    unsigned dst_g = GPixel_GetG(dst);
    unsigned dst_b = GPixel_GetB(dst);

    switch (mode) {
        case GBlendMode::kClear:
            return GPixel_PackARGB(0, 0, 0, 0); 

        case GBlendMode::kSrc:
            return src; 

        case GBlendMode::kDst:
            return dst; 

        case GBlendMode::kSrcOver:
            if (src_a == 0) return dst; 
            if (src_a == 255) return src; 
            return GPixel_PackARGB(
                src_a + div255(dst_a * inv_src_a),
                src_r + div255(dst_r * inv_src_a),
                src_g + div255(dst_g * inv_src_a),
                src_b + div255(dst_b * inv_src_a)
            );

        case GBlendMode::kDstOver:
            if (dst_a == 0) return src; 
            if (dst_a == 255) return dst;
            return GPixel_PackARGB(
                dst_a + div255(src_a * inv_dst_a),
                dst_r + div255(src_r * inv_dst_a),
                dst_g + div255(src_g * inv_dst_a),
                dst_b + div255(src_b * inv_dst_a)
            );

        case GBlendMode::kSrcIn:
            return GPixel_PackARGB(
                div255(dst_a * src_a),
                div255(dst_a * src_r),
                div255(dst_a * src_g),
                div255(dst_a * src_b)
            );

        case GBlendMode::kDstIn:
            return GPixel_PackARGB(
                div255(src_a * dst_a),
                div255(src_a * dst_r),
                div255(src_a * dst_g),
                div255(src_a * dst_b)
            );

        case GBlendMode::kSrcOut:
            return GPixel_PackARGB(
                div255(inv_dst_a * src_a),
                div255(inv_dst_a * src_r),
                div255(inv_dst_a * src_g),
                div255(inv_dst_a * src_b)
            );

        case GBlendMode::kDstOut:
            return GPixel_PackARGB(
                div255(inv_src_a * dst_a),
                div255(inv_src_a * dst_r),
                div255(inv_src_a * dst_g),
                div255(inv_src_a * dst_b)
            );

        case GBlendMode::kSrcATop:
            return GPixel_PackARGB(
                div255(dst_a * src_a) + div255(inv_src_a * dst_a),
                div255(dst_a * src_r) + div255(inv_src_a * dst_r),
                div255(dst_a * src_g) + div255(inv_src_a * dst_g),
                div255(dst_a * src_b) + div255(inv_src_a * dst_b)
            );

        case GBlendMode::kDstATop:
            return GPixel_PackARGB(
                div255(src_a * dst_a) + div255(inv_dst_a * src_a),
                div255(src_a * dst_r) + div255(inv_dst_a * src_r),
                div255(src_a * dst_g) + div255(inv_dst_a * src_g),
                div255(src_a * dst_b) + div255(inv_dst_a * src_b)
            );

        case GBlendMode::kXor:
            return GPixel_PackARGB(
                div255((inv_src_a * dst_a) + (inv_dst_a * src_a)),
                div255((inv_src_a * dst_r) + (inv_dst_a * src_r)),
                div255((inv_src_a * dst_g) + (inv_dst_a * src_g)),
                div255((inv_src_a * dst_b) + (inv_dst_a * src_b))
            );

        default:
            return src; 
    }
}

struct Edge {
    float x;    
    float slope;    
    int maxY;       
};


bool edgeCompare(const Edge& e1, const Edge& e2) {
    return e1.x < e2.x;
}

class MyGCanvas : public GCanvas {
public:
    MyGCanvas(const GBitmap& bitmap) : bitmap_(bitmap) {
        ctm = GMatrix();
        savedMatrices.push(ctm);
    }

    void save() override {
        savedMatrices.push(ctm);
    }

    void restore() override {
        if (!savedMatrices.empty()) {
            ctm = savedMatrices.top();
            savedMatrices.pop();
        } else {
            fprintf(stderr, "Error: Cannot restore beyond initial state!\n");
        }
    }

    void concat(const GMatrix& matrix) override {

        ctm = GMatrix::Concat(ctm, matrix);
    }

    void clear(const GColor& color) override {
        GPixel clearValue = premultiplyColor(color);
        int width = bitmap_.width();
        int height = bitmap_.height();

        for (int y = 0; y < height; ++y) {
            GPixel* row_addr = bitmap_.getAddr(0, y);
            for (int x = 0; x < width; ++x) {
                row_addr[x] = clearValue;
            }
        }
    }

        struct Edge {
            int yMin;
            int yMax;
            float x;
            float slope;
            int wind;

            bool isValid(int y) const {
                return y >= yMin && y < yMax;
            }

            float computeX(float y) const {
                return x + (y - yMin) * slope;
            }
        };

void drawPath(const GPath& path, const GPaint& paint) override {
    std::vector<Edge> edges;
    GPath::Edger edger(path);

    GPoint pts[GPath::kMaxNextPoints];
    while (auto verb = edger.next(pts)) {
        switch (verb.value()) {
            case GPathVerb::kLine:
                ctm.mapPoints(pts, pts, 2);
                addLineToEdges(pts[0], pts[1], edges);
                break;
            case GPathVerb::kQuad:
                ctm.mapPoints(pts, pts, 3);
                flattenQuad(pts, 0.25, edges);
                break;
            case GPathVerb::kCubic:
                ctm.mapPoints(pts, pts, 4);
                flattenCubic(pts, 0.25, edges);
                break;
            default:
                break;
        }
    }

    if (edges.empty()) return;

    std::sort(edges.begin(), edges.end(), [](const Edge& e1, const Edge& e2) {
        if (e1.yMin != e2.yMin) {
            return e1.yMin < e2.yMin;
        }
        return e1.x < e2.x;
    });

    int y = edges.front().yMin;
    int yMax = edges.back().yMax;

    std::vector<Edge*> activeEdges;

    GShader* shader = paint.peekShader();
    GBlendMode mode = paint.getBlendMode();
    GPixel src_pixel = premultiplyColor(paint.getColor());

    if (shader && !shader->setContext(ctm)) {
        return;
    }

    size_t edgeIndex = 0;

    for (; y < yMax; ++y) {
        if (y < 0 || y >= bitmap_.height()) {
            continue;
        }

        while (edgeIndex < edges.size() && edges[edgeIndex].yMin == y) {
            activeEdges.push_back(&edges[edgeIndex]);
            ++edgeIndex;
        }

        activeEdges.erase(std::remove_if(activeEdges.begin(), activeEdges.end(),
                                         [y](Edge* e) { return y >= e->yMax; }),
                          activeEdges.end());

        std::sort(activeEdges.begin(), activeEdges.end(), [y](Edge* e1, Edge* e2) {
            float x1 = e1->computeX(y + 0.5f);
            float x2 = e2->computeX(y + 0.5f);
            return x1 < x2;
        });

        int w = 0;

        for (size_t i = 0; i < activeEdges.size();) {
            Edge* e = activeEdges[i];
            float x = e->computeX(y + 0.5f);
            int xStart = GRoundToInt(x);

            w += e->wind;

            ++i;

            while (i < activeEdges.size()) {
                Edge* eNext = activeEdges[i];
                float xNext = eNext->computeX(y + 0.5f);
                int xEnd = GRoundToInt(xNext);

                if (w != 0) {
                    int left = std::max(0, xStart);
                    int right = std::min(bitmap_.width(), xEnd);

                    if (left < right) {
                        GPixel* row = bitmap_.getAddr(left, y);
                        int count = right - left;

                        if (shader) {
                            std::vector<GPixel> shaderRow(count);
                            shader->shadeRow(left, y, count, shaderRow.data());

                            for (int k = 0; k < count; ++k) {
                                row[k] = blendPixel(shaderRow[k], row[k], mode);
                            }
                        } else {
                            for (int k = 0; k < count; ++k) {
                                row[k] = blendPixel(src_pixel, row[k], mode);
                            }
                        }
                    }
                }

                w += eNext->wind;
                xStart = xEnd;
                ++i;
            }
        }

        if (w != 0) {
            std::cerr << "Error: Winding count should be zero, but got w = " << w << " at y = " << y << std::endl;
        }
    }
}


    void drawRect(const GRect& rect, const GPaint& paint) override {
        GShader* shader = paint.peekShader();
        GBlendMode mode = paint.getBlendMode();

        if (shader) {
            GPoint points[4] = {
                { rect.left, rect.top },
                { rect.right, rect.top },
                { rect.right, rect.bottom },
                { rect.left, rect.bottom }
            };

            GPoint transformedPoints[4];
            ctm.mapPoints(transformedPoints, points, 4);

            drawLegacyConvexPolygon(transformedPoints, 4, paint);
        } else {
            GPixel src_pixel = premultiplyColor(paint.getColor());

            int clippedStartX = std::max(0, static_cast<int>(std::round(rect.left)));
            int clippedStartY = std::max(0, static_cast<int>(std::round(rect.top)));
            int clippedEndX = std::min(bitmap_.width(), static_cast<int>(std::round(rect.right)));
            int clippedEndY = std::min(bitmap_.height(), static_cast<int>(std::round(rect.bottom)));

            if (clippedStartX >= clippedEndX || clippedStartY >= clippedEndY) {
                return;
            }

            for (int y = clippedStartY; y < clippedEndY; ++y) {
                GPixel* row_addr = bitmap_.getAddr(clippedStartX, y);
                for (int x = 0; x < (clippedEndX - clippedStartX); ++x) {
                    row_addr[x] = blendPixel(src_pixel, row_addr[x], mode);
                }
            }
        }
    }

    void drawLegacyConvexPolygon(const GPoint pts[], int count, const GPaint& paint) {
        
        if (count < 3) {
            return; 
        }

        GShader* shader = paint.peekShader();
        GBlendMode mode = paint.getBlendMode();
        GPixel src_pixel = premultiplyColor(paint.getColor());

        int bitmapHeight = bitmap_.height();
        int bitmapWidth = bitmap_.width();

        struct Edge {
            float x; 
            float slope;
            int maxY; 
        };

        std::vector<std::vector<Edge>> edgeTable(bitmapHeight);

        for (int i = 0; i < count; ++i) {
            GPoint p0 = pts[i];
            GPoint p1 = pts[(i + 1) % count];

            if (p0.y == p1.y) continue;

            if (p0.y > p1.y) std::swap(p0, p1);

            if (p1.y < 0 || p0.y >= bitmapHeight) continue;

            float slope = (p1.x - p0.x) / (p1.y - p0.y);

            int minY = std::max(0, static_cast<int>(std::round(p0.y)));
            int maxY = std::min(bitmapHeight - 1, static_cast<int>(std::round(p1.y)) - 1);

            if (minY > maxY) continue;

            Edge edge = { p0.x + slope * (minY - p0.y), slope, maxY };
            edgeTable[minY].push_back(edge);
        }

        std::vector<Edge> activeEdges;

        if (shader && !shader->setContext(ctm)) {
            return; 
        }

        for (int y = 0; y < bitmapHeight; ++y) {
            for (const Edge& e : edgeTable[y]) {
                activeEdges.push_back(e);
            }

            activeEdges.erase(std::remove_if(activeEdges.begin(), activeEdges.end(),
                                            [y](const Edge& e) { return e.maxY < y; }),
                            activeEdges.end());

            std::sort(activeEdges.begin(), activeEdges.end(), [](const Edge& e1, const Edge& e2) {
                return e1.x < e2.x;
            });

            for (size_t i = 0; i < activeEdges.size(); i += 2) {
                if (i + 1 >= activeEdges.size()) break;

                int xStart = std::max(0, static_cast<int>(std::round(activeEdges[i].x)));
                int xEnd = std::min(bitmapWidth, static_cast<int>(std::round(activeEdges[i + 1].x)));

                if (xStart >= xEnd) continue;

                GPixel* row_addr = bitmap_.getAddr(xStart, y);

                if (shader) {
                    GPixel row[xEnd - xStart];
                    shader->shadeRow(xStart, y, xEnd - xStart, row);

                    for (int x = 0; x < (xEnd - xStart); ++x) {
                        row_addr[x] = blendPixel(row[x], row_addr[x], mode);
                    }
                } else {
                    for (int x = xStart; x < xEnd; ++x) {
                        row_addr[x - xStart] = blendPixel(src_pixel, row_addr[x - xStart], mode);
                    }
                }
            }

            for (Edge& e : activeEdges) {
                e.x += e.slope;
            }
        }
    }



    void drawConvexPolygon(const GPoint pts[], int count, const GPaint& paint) override {
        if (count < 3) {
            return; 
        }

        GPoint transformedPoints[count];
        ctm.mapPoints(transformedPoints, pts, count);

        GShader* shader = paint.peekShader();
        GBlendMode mode = paint.getBlendMode();
        GPixel src_pixel = premultiplyColor(paint.getColor());

        int bitmapHeight = bitmap_.height();
        int bitmapWidth = bitmap_.width();

        struct Edge {
            float x; 
            float slope;
            int maxY; 
        };

        std::vector<std::vector<Edge>> edgeTable(bitmapHeight);

        for (int i = 0; i < count; ++i) {
            GPoint p0 = transformedPoints[i];
            GPoint p1 = transformedPoints[(i + 1) % count];

            if (p0.y == p1.y) continue;

            if (p0.y > p1.y) std::swap(p0, p1);

            if (p1.y < 0 || p0.y >= bitmapHeight) continue;

            float slope = (p1.x - p0.x) / (p1.y - p0.y);

            int minY = std::max(0, static_cast<int>(std::round(p0.y)));
            int maxY = std::min(bitmapHeight - 1, static_cast<int>(std::round(p1.y)) - 1);

            if (minY > maxY) continue;

            Edge edge = { p0.x + slope * (minY - p0.y), slope, maxY };
            edgeTable[minY].push_back(edge);
        }

        std::vector<Edge> activeEdges;

        if (shader && !shader->setContext(ctm)) {
            return; 
        }

        for (int y = 0; y < bitmapHeight; ++y) {
            for (const Edge& e : edgeTable[y]) {
                activeEdges.push_back(e);
            }

            activeEdges.erase(std::remove_if(activeEdges.begin(), activeEdges.end(),
                                            [y](const Edge& e) { return e.maxY < y; }),
                            activeEdges.end());

            std::sort(activeEdges.begin(), activeEdges.end(), [](const Edge& e1, const Edge& e2) {
                return e1.x < e2.x;
            });

            for (size_t i = 0; i < activeEdges.size(); i += 2) {
                if (i + 1 >= activeEdges.size()) break;

                int xStart = std::max(0, static_cast<int>(std::round(activeEdges[i].x)));
                int xEnd = std::min(bitmapWidth, static_cast<int>(std::round(activeEdges[i + 1].x)));

                if (xStart >= xEnd) continue;

                GPixel* row_addr = bitmap_.getAddr(xStart, y);

                if (shader) {
                    GPixel row[xEnd - xStart];
                    shader->shadeRow(xStart, y, xEnd - xStart, row);

                    for (int x = 0; x < (xEnd - xStart); ++x) {
                        row_addr[x] = blendPixel(row[x], row_addr[x], mode);
                    }
                } else {
                    for (int x = xStart; x < xEnd; ++x) {
                        row_addr[x - xStart] = blendPixel(src_pixel, row_addr[x - xStart], mode);
                    }
                }
            }

            for (Edge& e : activeEdges) {
                e.x += e.slope;
            }
        }
    }

    void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count, const int indices[], const GPaint& paint) override {
        for (int i = 0; i < count; ++i) {
            int n = i * 3;
            GPoint p0 = verts[indices[n]];
            GPoint p1 = verts[indices[n + 1]];
            GPoint p2 = verts[indices[n + 2]];

            GColor c0 = colors ? colors[indices[n]] : GColor::RGBA(1.0f, 1.0f, 1.0f, 1.0f);
            GColor c1 = colors ? colors[indices[n + 1]] : GColor::RGBA(1.0f, 1.0f, 1.0f, 1.0f);
            GColor c2 = colors ? colors[indices[n + 2]] : GColor::RGBA(1.0f, 1.0f, 1.0f, 1.0f);

            GPoint t0 = texs ? texs[indices[n]] : GPoint{0, 0};
            GPoint t1 = texs ? texs[indices[n + 1]] : GPoint{0, 0};
            GPoint t2 = texs ? texs[indices[n + 2]] : GPoint{0, 0};

            drawTriangle(p0, p1, p2, c0, c1, c2, t0, t1, t2, paint);
        }
    }

    void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level, const GPaint& paint) {
        int quadCount = (level + 1) * (level + 1);
        std::vector<GPoint> tessVerts(quadCount);
        std::vector<GColor> tessColors(quadCount);
        std::vector<GPoint> tessTexs(quadCount);
        std::vector<int> indices;

        for (int y = 0; y <= level; ++y) {
            float ty = static_cast<float>(y) / level;
            for (int x = 0; x <= level; ++x) {
                float tx = static_cast<float>(x) / level;

                tessVerts[y * (level + 1) + x] = 
                    lerp(lerp(verts[0], verts[1], tx), lerp(verts[3], verts[2], tx), ty);

                if (colors) {
                    tessColors[y * (level + 1) + x] = 
                        lerp(lerp(colors[0], colors[1], tx), lerp(colors[3], colors[2], tx), ty);
                } else {
                    tessColors[y * (level + 1) + x] = GColor::RGBA(1.0f, 1.0f, 1.0f, 1.0f); 
                }

                if (texs) {
                    tessTexs[y * (level + 1) + x] = 
                        lerp(lerp(texs[0], texs[1], tx), lerp(texs[3], texs[2], tx), ty);
                } else {
                    tessTexs[y * (level + 1) + x] = GPoint{0, 0}; 
                }
            }
        }

        for (int y = 0; y < level; ++y) {
            for (int x = 0; x < level; ++x) {
                int idx = y * (level + 1) + x;


                indices.push_back(idx);
                indices.push_back(idx + 1);
                indices.push_back(idx + level + 2);

                indices.push_back(idx);
                indices.push_back(idx + level + 2);
                indices.push_back(idx + level + 1);
            }
        }

        drawMesh(tessVerts.data(), colors ? tessColors.data() : nullptr, texs ? tessTexs.data() : nullptr, indices.size() / 3, indices.data(), paint);
    }





private:
    GMatrix ctm;
    std::stack<GMatrix> savedMatrices;
    GBitmap bitmap_;

    void addLineToEdges(const GPoint& p0, const GPoint& p1, std::vector<Edge>& edges) {
        GPoint start = p0;
        GPoint end = p1;

        int wind = 1;
        if (start.y > end.y) {
            std::swap(start, end);
            wind = -1;
        }

        if (GRoundToInt(start.y) == GRoundToInt(end.y)) {
            return;
        }

        float dy = end.y - start.y;
        float dx = end.x - start.x;

        float slope = (std::abs(dy) < 1e-6) ? 0 : dx / dy;

        int yMin = GRoundToInt(start.y);
        int yMax = GRoundToInt(end.y);
        float x = start.x + (yMin + 0.5f - start.y) * slope;

        edges.push_back({yMin, yMax, x, slope, wind});
    }


    void flattenQuad(const GPoint pts[3], float tolerance, std::vector<Edge>& edges) {
        GPoint dst[5];
        GPath::ChopQuadAt(pts, dst, 0.5f);

        if (needsSubdivision(dst, tolerance)) {
            flattenQuad(dst, tolerance, edges);
            flattenQuad(dst + 2, tolerance, edges);
        } else {
            if (GRoundToInt(pts[0].y) != GRoundToInt(pts[2].y)) {
                addLineToEdges(pts[0], pts[2], edges);
            }
        }
    }

    void flattenCubic(const GPoint pts[4], float tolerance, std::vector<Edge>& edges) {
        GPoint dst[7];
        GPath::ChopCubicAt(pts, dst, 0.5f);

        if (needsSubdivision(dst, tolerance)) {
            flattenCubic(dst, tolerance, edges);
            flattenCubic(dst + 3, tolerance, edges);
        } else {
            if (GRoundToInt(pts[0].y) != GRoundToInt(pts[3].y)) {
                addLineToEdges(pts[0], pts[3], edges);
            }
        }
    }

    bool needsSubdivision(const GPoint pts[], float tolerance) {
        GPoint mid1 = (pts[0] + pts[1]) * 0.5f;
        GPoint mid2 = (pts[1] + pts[2]) * 0.5f;
        GPoint overallMid = (mid1 + mid2) * 0.5f;

        float dist1 = (mid1 - pts[1]).length();
        float dist2 = (mid2 - pts[1]).length();
        return dist1 > tolerance || dist2 > tolerance;
    }

    void drawTriangle(const GPoint& p0, const GPoint& p1, const GPoint& p2, const GColor& c0, const GColor& c1, const GColor& c2,
        const GPoint& t0, const GPoint& t1, const GPoint& t2, const GPaint& paint) {
        GPoint points[3] = { p0, p1, p2 };
        GPoint transformedPoints[3];
        ctm.mapPoints(transformedPoints, points, 3);

        float minX = std::min({transformedPoints[0].x, transformedPoints[1].x, transformedPoints[2].x});
        float minY = std::min({transformedPoints[0].y, transformedPoints[1].y, transformedPoints[2].y});
        float maxX = std::max({transformedPoints[0].x, transformedPoints[1].x, transformedPoints[2].x});
        float maxY = std::max({transformedPoints[0].y, transformedPoints[1].y, transformedPoints[2].y});

        minX = std::max(0.0f, std::floor(minX));
        minY = std::max(0.0f, std::floor(minY));
        maxX = std::min(static_cast<float>(bitmap_.width() - 1), std::ceil(maxX));
        maxY = std::min(static_cast<float>(bitmap_.height() - 1), std::ceil(maxY));

        GMatrix M(
            transformedPoints[1].x - transformedPoints[0].x, transformedPoints[2].x - transformedPoints[0].x, transformedPoints[0].x,
            transformedPoints[1].y - transformedPoints[0].y, transformedPoints[2].y - transformedPoints[0].y, transformedPoints[0].y
        );

        std::optional<GMatrix> optInvM = M.invert();
        if (!optInvM.has_value()) {
            return;
        }
        GMatrix invM = optInvM.value();

        GShader* shader = paint.peekShader();
        if (shader) {
            GMatrix texMatrix(
                t1.x - t0.x, t2.x - t0.x, t0.x,
                t1.y - t0.y, t2.y - t0.y, t0.y
            );

            GMatrix combinedMatrix = texMatrix * invM;
            if (!shader->setContext(combinedMatrix)) {
                return;
            }
        }

        GBlendMode blendMode = paint.getBlendMode();

        for (int y = static_cast<int>(minY); y <= static_cast<int>(maxY); ++y) {
            for (int x = static_cast<int>(minX); x <= static_cast<int>(maxX); ++x) {
                GPoint pixelPoint = {x + 0.5f, y + 0.5f};
                GPoint barycentric;
                invM.mapPoints(&barycentric, &pixelPoint, 1);
                float u = barycentric.x;
                float v = barycentric.y;
                float w = 1 - u - v;

                if (u >= 0 && v >= 0 && w >= 0) {
                    GColor interpolatedColor = {
                        w * c0.r + u * c1.r + v * c2.r,
                        w * c0.g + u * c1.g + v * c2.g,
                        w * c0.b + u * c1.b + v * c2.b,
                        w * c0.a + u * c1.a + v * c2.a
                    };

                    GPixel blendedPixel;
                    if (shader) {
                        GPixel shaderPixel;
                        shader->shadeRow(x, y, 1, &shaderPixel);

                        blendedPixel = premultiplyColor(interpolatedColor);
                        blendedPixel = blendPixel(shaderPixel, blendedPixel, blendMode);
                    } else {
                        blendedPixel = premultiplyColor(interpolatedColor);
                    }

                    GPixel* destPixel = bitmap_.getAddr(x, y);
                    *destPixel = blendPixel(blendedPixel, *destPixel, blendMode);
                }
            }
        }
    }

    GPoint lerp(const GPoint& a, const GPoint& b, float t) const {
        return {a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
    }


    GColor lerp(const GColor& a, const GColor& b, float t) const {
        return GColor::RGBA(
            a.r + t * (b.r - a.r),
            a.g + t * (b.g - a.g),
            a.b + t * (b.b - a.b),
            a.a + t * (b.a - a.a)
        );
    }
};



std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& bitmap) {
    return std::make_unique<MyGCanvas>(bitmap);
}