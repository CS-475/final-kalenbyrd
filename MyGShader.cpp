#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
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

class MyBitmapShader : public GShader {
public:
    MyBitmapShader(const GBitmap& bitmap, const GMatrix& localMatrix, GTileMode tileMode)
        : bitmap_(bitmap), localMatrix_(localMatrix), tileMode_(tileMode) {}

    bool isOpaque() override {
        return bitmap_.isOpaque();
    }

    bool setContext(const GMatrix& ctm) override {
        GMatrix combinedMatrix = GMatrix::Concat(ctm, localMatrix_);
        if (auto inv = combinedMatrix.invert()) {
            inverse_ = *inv;
            return true;
        } else {
            return false;
        }
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {
        GPoint srcPoints[count];
        for (int i = 0; i < count; ++i) {
            srcPoints[i] = { x + i + 0.5f, y + 0.5f };
        }

        inverse_.mapPoints(srcPoints, srcPoints, count);

        int bmWidth = bitmap_.width();
        int bmHeight = bitmap_.height();

        for (int i = 0; i < count; ++i) {
            float sx = srcPoints[i].x;
            float sy = srcPoints[i].y;

            switch (tileMode_) {
                case GTileMode::kClamp:
                    sx = std::clamp(sx, 0.0f, bmWidth - 1.0f);
                    sy = std::clamp(sy, 0.0f, bmHeight - 1.0f);
                    break;
                case GTileMode::kRepeat:
                    sx = sx - bmWidth * std::floor(sx / bmWidth);
                    sy = sy - bmHeight * std::floor(sy / bmHeight);
                    break;
                case GTileMode::kMirror:
                    sx = std::abs(sx) - bmWidth * std::floor(std::abs(sx) / bmWidth);
                    sy = std::abs(sy) - bmHeight * std::floor(std::abs(sy) / bmHeight);
                    if (static_cast<int>(std::floor(sx / bmWidth)) % 2 == 1) {
                        sx = bmWidth - sx;
                    }
                    if (static_cast<int>(std::floor(sy / bmHeight)) % 2 == 1) {
                        sy = bmHeight - sy;
                    }
                    break;
            }

            int ix = static_cast<int>(std::round(sx));
            int iy = static_cast<int>(std::round(sy));

            ix = std::max(0, std::min(ix, bmWidth - 1)); 
            iy = std::max(0, std::min(iy, bmHeight - 1));

            row[i] = *bitmap_.getAddr(ix, iy);
        }
    }

private:
    GBitmap bitmap_;
    GMatrix localMatrix_;
    GMatrix inverse_; 
    GTileMode tileMode_;
};

std::shared_ptr<GShader> GCreateBitmapShader(const GBitmap& bitmap, const GMatrix& localMatrix, GTileMode tileMode) {
    if (bitmap.width() <= 0 || bitmap.height() <= 0) {
        return nullptr;
    }
    return std::make_shared<MyBitmapShader>(bitmap, localMatrix, tileMode);
}

float clamp(float val, float minVal, float maxVal) {
    return std::max(minVal, std::min(val, maxVal));
}

class GradientGShader : public GShader {
public:
    GradientGShader(GPoint p0, GPoint p1, const GColor* colors, int count, GTileMode tileMode)
        : p0_(p0), p1_(p1), colors_(colors), count_(count), tileMode_(tileMode) {}

    bool isOpaque() override {
        for (int i = 0; i < count_; ++i) {
            if (colors_[i].a < 1.0f) {
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
        } else {
            return false;
        }
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {
        GPoint srcPoints[count];
        for (int i = 0; i < count; ++i) {
            srcPoints[i] = { x + i + 0.5f, y + 0.5f };
        }

        inverse_.mapPoints(srcPoints, srcPoints, count);

        for (int i = 0; i < count; ++i) {
            float t = (srcPoints[i].x - p0_.x) / (p1_.x - p0_.x);

            switch (tileMode_) {
                case GTileMode::kClamp:
                    t = clamp(t, 0.0f, 1.0f);
                    break;
                case GTileMode::kRepeat:
                    t = t - std::floor(t);
                    break;
                case GTileMode::kMirror:
                    t = std::abs(t - 2 * std::floor(t / 2));
                    if (static_cast<int>(std::floor(t)) % 2 == 1) {
                        t = 1.0f - t;
                    }
                    break;
            }

            int idx = static_cast<int>(t * (count_ - 1));
            GColor blendedColor = lerpColor(colors_[idx], colors_[idx + 1], t * (count_ - 1) - idx);
            row[i] = premultiplyColor(blendedColor);
        }
    }

private:
    GColor lerpColor(const GColor& c0, const GColor& c1, float t) {
        return {
            (1 - t) * c0.r + t * c1.r,
            (1 - t) * c0.g + t * c1.g,
            (1 - t) * c0.b + t * c1.b,
            (1 - t) * c0.a + t * c1.a
        };
    }

    GPoint p0_, p1_;
    const GColor* colors_;
    int count_;
    GMatrix inverse_;
    GTileMode tileMode_;
};

std::shared_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode tileMode) {
    if (count < 1) {
        return nullptr;
    }
    return std::make_shared<GradientGShader>(p0, p1, colors, count, tileMode);
}

class ComposShader : public GShader {
public:
    ComposShader(std::shared_ptr<GShader> shaderA, std::shared_ptr<GShader> shaderB)
        : fShaderA(std::move(shaderA)), fShaderB(std::move(shaderB)) {}

    bool isOpaque() override {
        if (!fShaderA || !fShaderB) {
            std::cerr << "Error: One of the shaders is null in isOpaque." << std::endl;
            return false;
        }
        return fShaderA->isOpaque() && fShaderB->isOpaque();
    }

    bool setContext(const GMatrix& ctm) override {
        if (!fShaderA || !fShaderB) {
            std::cerr << "Error: One of the shaders is null in setContext." << std::endl;
            return false;
        }
        bool contextA = fShaderA->setContext(ctm);
        bool contextB = fShaderB->setContext(ctm);
        if (!contextA || !contextB) {
            std::cerr << "Error: Failed to set context for one of the shaders." << std::endl;
        }
        return contextA && contextB;
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {
        if (!fShaderA || !fShaderB) {
            std::cerr << "Error: One of the shaders is null in shadeRow." << std::endl;
            std::fill(row, row + count, GPixel_PackARGB(0, 0, 0, 0));
            return;
        }

        if (count <= 0) {
            std::cerr << "Error: Invalid count value in shadeRow: " << count << std::endl;
            return;
        }

        std::vector<GPixel> rowA(count);
        std::vector<GPixel> rowB(count);

        fShaderA->shadeRow(x, y, count, rowA.data());
        fShaderB->shadeRow(x, y, count, rowB.data());

        for (int i = 0; i < count; ++i) {
            row[i] = blendPixel(rowA[i], rowB[i]);
        }
    }

private:
    std::shared_ptr<GShader> fShaderA;
    std::shared_ptr<GShader> fShaderB;

    GPixel blendPixel(const GPixel& a, const GPixel& b) {
        int aR = GPixel_GetR(a);
        int aG = GPixel_GetG(a);
        int aB = GPixel_GetB(a);
        int aA = GPixel_GetA(a);

        int bR = GPixel_GetR(b);
        int bG = GPixel_GetG(b);
        int bB = GPixel_GetB(b);
        int bA = GPixel_GetA(b);

        // Blend components using source-over blending
        int blendedA = aA + (bA * (255 - aA)) / 255;
        int blendedR = (aR * aA + bR * bA * (255 - aA) / 255) / 255;
        int blendedG = (aG * aA + bG * bA * (255 - aA) / 255) / 255;
        int blendedB = (aB * aA + bB * bA * (255 - aA) / 255) / 255;

        return GPixel_PackARGB(std::clamp(blendedA, 0, 255), std::clamp(blendedR, 0, 255), std::clamp(blendedG, 0, 255), std::clamp(blendedB, 0, 255));
    }
};