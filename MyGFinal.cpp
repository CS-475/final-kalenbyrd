#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GFinal.h"
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>

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


class MyGFinal : public GFinal {
    public:

    /**
     * The voronoi shader is defined by an array of points, each with an associated color.
     * The color at any (x,y) is the color of the closest point from the array.
     */
    std::shared_ptr<GShader> createVoronoiShader(const GPoint points[],
                                                 const GColor colors[],
                                                 int count) override {
        if (count < 1 || points == nullptr || colors == nullptr) {
            return nullptr;
        }

        class VoronoiShader : public GShader {
        public:
            VoronoiShader(const GPoint points[], const GColor colors[], int count)
                : fPoints(points, points + count), fColors(colors, colors + count) {
            }

            bool isOpaque() override {
                for (const auto& color : fColors) {
                    if (color.a < 1.0f) {
                        return false;
                    }
                }
                return true;
            }

            bool setContext(const GMatrix& ctm) override {
                if (auto inverted = ctm.invert()) {
                    fInverse = *inverted;
                    return true;
                }
                return false;
            }

            void shadeRow(int x, int y, int count, GPixel row[]) override {

                std::vector<GPoint> points(count);
                for (int i = 0; i < count; ++i) {
                    points[i] = {x + i + 0.5f, y + 0.5f};
                }
                
                std::vector<GPoint> localPoints(count);
                fInverse.mapPoints(localPoints.data(), points.data(), count);

                for (int i = 0; i < count; ++i) {
                    const GPoint& src = localPoints[i];
                    int closestIndex = 0;
                    float minD = std::numeric_limits<float>::max();

                    for (int j = 0; j < fPoints.size(); ++j) {
                        float dx = src.x - fPoints[j].x;
                        float dy = src.y - fPoints[j].y;
                        float dist = dx * dx + dy * dy;
                        if (dist < minD) {
                            minD = dist;
                            closestIndex = j;
                        }
                    }

                    row[i] = premultiplyColor(fColors[closestIndex]);
                }
            }

        private:
            std::vector<GPoint> fPoints;
            std::vector<GColor> fColors;
            GMatrix fInverse;
        };

        return std::make_shared<VoronoiShader>(points, colors, count);
    }

        /**
     *  Return a sweep-gradient shader, centered at 'center', starting wiht color[0] at  startRadians
     *  and ending with colors[count-1] at startRadians+2pi. The colors are distributed evenly around the sweep.
     */
    std::shared_ptr<GShader> createSweepGradient(GPoint center, float startRadians,
                                                         const GColor colors[], int count) {
        return nullptr;
    }

        /*
     *  Returns a new type of linear gradient. In this variant, the "count" colors are
     *  positioned along the line p0...p1 not "evenly", but according to the pos[] array.
     *
     *  pos[] holds "count" values, each 0...1, which specify the percentage along the
     *  line where the each color lies.
     *
     *  e.g. pos[] = {0, 0.25, 1} would mean 3 colors positioned as follows:
     *
     *      color[0] ..... color[1] ..... ..... ..... color[2]
     *
     *  color[i] is positioned by computing (1 - pos[i])*p0 + pos[i]*p1
     *
     *  For this API, pos[] will always be monotonic, with p[0] == 0 and p[count-1] == 1.0
     *
     *  For simplicity, assume that we're using "clamp" tiling.
     */
    std::shared_ptr<GShader> createLinearPosGradient(GPoint p0, GPoint p1,
                                                             const GColor colors[],
                                                             const float pos[],
                                                             int count) {
        if (count < 2 || !colors || !pos) {
            return nullptr;
        }

        class LinearPosGradient : public GShader {
        public:
            LinearPosGradient(GPoint p0, GPoint p1, const GColor colors[], const float pos[], int count)
                : p0_(p0), p1_(p1), colors_(colors, colors + count), positions_(pos, pos + count) {}

            bool isOpaque() override {
                for (const auto& color : colors_) {
                    if (color.a < 1.0f) {
                        return false;
                    }
                }
                return true;
            }

            bool setContext(const GMatrix& ctm) override {
                float dx = p1_.x - p0_.x;
                float dy = p1_.y - p0_.y;
                float D = std::sqrt(dx * dx + dy * dy);

                if (D == 0) {
                    return false;
                }

                GMatrix gradientMatrix = GMatrix(
                    dx / D, -dy / D, p0_.x,
                    dy / D,  dx / D, p0_.y
                );

                GMatrix combinedMatrix = GMatrix::Concat(ctm, gradientMatrix);
                if (auto inv = combinedMatrix.invert()) {
                    inverse_ = *inv;
                    return true;
                }
                return false;
            }

            void shadeRow(int x, int y, int count, GPixel row[]) override {
                std::vector<GPoint> points(count);
                for (int i = 0; i < count; ++i) {
                    points[i] = {x + i + 0.5f, y + 0.5f};
                }

                std::vector<GPoint> localPoints(count);
                inverse_.mapPoints(localPoints.data(), points.data(), count);

                for (int i = 0; i < count; ++i) {
                    float t = (localPoints[i].x - p0_.x) / (p1_.x - p0_.x);
                    t = std::max(0.0f, std::min(1.0f, t));

                    size_t j;
                    for (j = 0; j < positions_.size() - 1; ++j) {
                        if (t <= positions_[j + 1]) {
                            break;
                        }
                    }

                    float segment_t = (t - positions_[j]) / (positions_[j + 1] - positions_[j]);
                    GColor inter = {
                        (1 - segment_t) * colors_[j].r + segment_t * colors_[j + 1].r,
                        (1 - segment_t) * colors_[j].g + segment_t * colors_[j + 1].g,
                        (1 - segment_t) * colors_[j].b + segment_t * colors_[j + 1].b,
                        (1 - segment_t) * colors_[j].a + segment_t * colors_[j + 1].a
                    };

                    row[i] = premultiplyColor(inter);
                }
            }

        private:
            GPoint p0_, p1_;
            std::vector<GColor> colors_;
            std::vector<float> positions_;
            GMatrix inverse_;
        };

        return std::make_shared<LinearPosGradient>(p0, p1, colors, pos, count);
    }


        /*
     *  Returns an instance to a shader that will proxy to a "realShader", and transform
     *  its output using the GColorMatrix provided.
     *
     *  Note: the GColorMatrix is defined to operate on unpremul GColors
     *
     *  Note: the resulting colors (after applying the colormatrix) may be out of bounds
     *        for color componets. If this happens they should be clamped to legal values.
     */
    std::shared_ptr<GShader> createColorMatrixShader(const GColorMatrix&,
                                                             GShader* realShader) {
        return nullptr;
    }

        /**
     *  Construct a path that, when drawn, will look like a stroke of the specified polygon.
     *  - count is the number of points in the polygon (it will be >= 2)
     *  - width is the thickness of the stroke that should be centered on the polygon edges
     *  - isClosed specifies if the polygon should appear closed (true) or open (false).
     *
     *  Any caps or joins needed should be round (circular).
     */
    std::shared_ptr<GPath> strokePolygon(const GPoint[], int count, float width, bool isClosed) {
        return nullptr;
    }


        /*
     *  Draw the corresponding mesh constructed from a quad with each side defined by a
     *  quadratic bezier, evaluating them to produce "level" interior lines (same convention
     *  as drawQuad().
     *
     *  pts[0]    pts[1]    pts[2]
     *  pts[7]              pts[3]
     *  pts[6]    pts[5]    pts[4]
     *
     *  Evaluate values within the mesh using the Coons Patch formulation:
     *
     *  value(u,v) = TB(u,v) + LR(u,v) - Corners(u,v)
     *
     *     Top is quadratic bezier: pts[0], pts[1], pts[2]
     *  Bottom is quadratic bezier: pts[6], pts[5], pts[4]
     *    Left is quadratic bezier: pts[0], pts[7], pts[6]
     *   Right is quadratic bezier: pts[2], pts[3], pts[4]
     *
     *  Where
     *      TB is computed by first evaluating the Top and Bottom curves at (u), and then
     *      linearly interpolating those points by (v)
     *
     *      LR is computed by first evaluating the Left and Right curves at (v), and then
     *      linearly interpolating those points by (u)
     *
     *      Corners is computed by our standard "drawQuad" evaluation using the 4 corners 0,2,4,6
     */
    void drawQuadraticCoons(GCanvas*, const GPoint pts[8], const GPoint tex[4],
                                    int level, const GPaint&) {

                                    }
    
};




std::unique_ptr<GFinal> GCreateFinal() {
    return std::make_unique<MyGFinal>();
}

