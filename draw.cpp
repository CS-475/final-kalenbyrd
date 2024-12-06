#include "include/GCanvas.h"
#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GRandom.h"
#include "include/GShader.h"
#include "include/GMatrix.h"
#include "include/GPath.h"
#include <cmath>

std::string GDrawSomething(GCanvas* canvas, GISize dimension) {

    canvas->clear({1, 1, 1, 1});


    GPoint p0 = {0, 0};
    GPoint p1 = {static_cast<float>(dimension.width), static_cast<float>(dimension.height)};

    GColor colors[] = {
        {1, 1, 0, 0},  
        {1, 0, 1, 0},  
        {1, 0, 0, 1},
        {1, 1, 1, 0},  
    };
    const int colorCount = sizeof(colors) / sizeof(GColor);

    auto shader = GCreateLinearGradient(p0, p1, colors, colorCount);


    GPaint paint(shader);

    GRect rect = GRect::WH(dimension.width, dimension.height);
    canvas->drawRect(rect, paint);

    GPoint center = {dimension.width / 2.0f, dimension.height / 2.0f};
    float maxRadius = std::min(dimension.width, dimension.height) / 2.0f;

    const int numShapes = 30;
    GRandom rand;

    for (int i = 0; i < numShapes; ++i) {
        float angle = i * 2 * M_PI / numShapes;
        float nextAngle = (i + 1) * 2 * M_PI / numShapes;

        float innerRadius = maxRadius * rand.nextF() * 0.5f;

        GPoint p0 = {
            center.x + innerRadius * std::cos(angle),
            center.y + innerRadius * std::sin(angle)
        };
        GPoint p1 = {
            center.x + maxRadius * std::cos(angle),
            center.y + maxRadius * std::sin(angle)
        };
        GPoint p2 = {
            center.x + maxRadius * std::cos(nextAngle),
            center.y + maxRadius * std::sin(nextAngle)
        };
        GPoint p3 = {
            center.x + innerRadius * std::cos(nextAngle),
            center.y + innerRadius * std::sin(nextAngle)
        };

        GPoint quad[4] = {p0, p1, p2, p3};

        GColor color = {
            1.0f,          
            rand.nextF(), 
            rand.nextF(),  
            rand.nextF()  
        };

        GPaint shapePaint(color);

        canvas->drawConvexPolygon(quad, 4, shapePaint);
    }

    return "Gradient Background with Circular Pattern";
}