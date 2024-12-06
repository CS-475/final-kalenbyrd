#include "include/GPath.h"
#include "include/GPathBuilder.h"
#include "include/GRect.h"
#include <algorithm>
#include <cassert>


/*
GRect GPath::bounds() const {
    if (fPts.empty()) {
        return GRect::LTRB(0, 0, 0, 0);
    }

    float minX = fPts[0].x;
    float maxX = fPts[0].x;
    float minY = fPts[0].y;
    float maxY = fPts[0].y;

    GPath::Iter iter(*this);
    GPoint pts[GPath::kMaxNextPoints];

    while (auto verb = iter.next(pts)) {
        switch (verb.value()) {
            case GPathVerb::kMove:
                // Move doesn't contribute to bounds
                break;
            case GPathVerb::kLine:
                for (int i = 0; i < 2; ++i) {
                    minX = std::min(minX, pts[i].x);
                    maxX = std::max(maxX, pts[i].x);
                    minY = std::min(minY, pts[i].y);
                    maxY = std::max(maxY, pts[i].y);
                }
                break;
            case GPathVerb::kQuad: {
                // Endpoint bounds
                for (int i = 0; i < 3; ++i) {
                    minX = std::min(minX, pts[i].x);
                    maxX = std::max(maxX, pts[i].x);
                    minY = std::min(minY, pts[i].y);
                    maxY = std::max(maxY, pts[i].y);
                }

                // Calculate tight bounds using extremum points
                float aX = pts[0].x - 2 * pts[1].x + pts[2].x;
                float bX = 2 * (pts[1].x - pts[0].x);
                if (fabs(aX) > 1e-6) {
                    float t = -bX / (2 * aX);
                    if (t > 0 && t < 1) {
                        float x = (1 - t) * (1 - t) * pts[0].x + 2 * (1 - t) * t * pts[1].x + t * t * pts[2].x;
                        minX = std::min(minX, x);
                        maxX = std::max(maxX, x);
                    }
                }

                float aY = pts[0].y - 2 * pts[1].y + pts[2].y;
                float bY = 2 * (pts[1].y - pts[0].y);
                if (fabs(aY) > 1e-6) {
                    float t = -bY / (2 * aY);
                    if (t > 0 && t < 1) {
                        float y = (1 - t) * (1 - t) * pts[0].y + 2 * (1 - t) * t * pts[1].y + t * t * pts[2].y;
                        minY = std::min(minY, y);
                        maxY = std::max(maxY, y);
                    }
                }
                break;
            }
            case GPathVerb::kCubic: {
                // Endpoint bounds
                for (int i = 0; i < 4; ++i) {
                    minX = std::min(minX, pts[i].x);
                    maxX = std::max(maxX, pts[i].x);
                    minY = std::min(minY, pts[i].y);
                    maxY = std::max(maxY, pts[i].y);
                }

                // Calculate tight bounds using extremum points
                float aX = -pts[0].x + 3 * pts[1].x - 3 * pts[2].x + pts[3].x;
                float bX = 2 * (pts[0].x - 2 * pts[1].x + pts[2].x);
                float cX = -pts[0].x + pts[1].x;
                float discriminantX = bX * bX - 4 * aX * cX;
                if (aX != 0 && discriminantX >= 0) {
                    float sqrtDisc = sqrt(discriminantX);
                    float t1 = (-bX + sqrtDisc) / (2 * aX);
                    float t2 = (-bX - sqrtDisc) / (2 * aX);
                    for (float t : {t1, t2}) {
                        if (t > 0 && t < 1) {
                            float x = powf(1 - t, 3) * pts[0].x + 3 * powf(1 - t, 2) * t * pts[1].x + 3 * (1 - t) * t * t * pts[2].x + powf(t, 3) * pts[3].x;
                            minX = std::min(minX, x);
                            maxX = std::max(maxX, x);
                        }
                    }
                }

                float aY = -pts[0].y + 3 * pts[1].y - 3 * pts[2].y + pts[3].y;
                float bY = 2 * (pts[0].y - 2 * pts[1].y + pts[2].y);
                float cY = -pts[0].y + pts[1].y;
                float discriminantY = bY * bY - 4 * aY * cY;
                if (aY != 0 && discriminantY >= 0) {
                    float sqrtDisc = sqrt(discriminantY);
                    float t1 = (-bY + sqrtDisc) / (2 * aY);
                    float t2 = (-bY - sqrtDisc) / (2 * aY);
                    for (float t : {t1, t2}) {
                        if (t > 0 && t < 1) {
                            float y = powf(1 - t, 3) * pts[0].y + 3 * powf(1 - t, 2) * t * pts[1].y + 3 * (1 - t) * t * t * pts[2].y + powf(t, 3) * pts[3].y;
                            minY = std::min(minY, y);
                            maxY = std::max(maxY, y);
                        }
                    }
                }
                break;
            }
        }
    }

    return GRect::LTRB(minX, minY, maxX, maxY);
}
*/

GRect GPath::bounds() const {
    if (fPts.empty()) {
        return GRect::LTRB(0, 0, 0, 0);
    }

    float minX = fPts[0].x;
    float maxX = fPts[0].x;
    float minY = fPts[0].y;
    float maxY = fPts[0].y;

    GPath::Iter iter(*this);
    GPoint pts[GPath::kMaxNextPoints];

    while (auto verb = iter.next(pts)) {
        switch (verb.value()) {
            case GPathVerb::kMove:
                break;
            case GPathVerb::kLine:
                for (int i = 0; i < 2; ++i) {
                    minX = std::min(minX, pts[i].x);
                    maxX = std::max(maxX, pts[i].x);
                    minY = std::min(minY, pts[i].y);
                    maxY = std::max(maxY, pts[i].y);
                }
                break;
            case GPathVerb::kQuad: {
                for (int i = 0; i < 3; ++i) {
                    minX = std::min(minX, pts[i].x);
                    maxX = std::max(maxX, pts[i].x);
                    minY = std::min(minY, pts[i].y);
                    maxY = std::max(maxY, pts[i].y);
                }
                break;
            }
            case GPathVerb::kCubic: {
                for (int i = 0; i < 4; ++i) {
                    minX = std::min(minX, pts[i].x);
                    maxX = std::max(maxX, pts[i].x);
                    minY = std::min(minY, pts[i].y);
                    maxY = std::max(maxY, pts[i].y);
                }
                break;
            }
        }
    }

    return GRect::LTRB(minX, minY, maxX, maxY);
} 


void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t) {
    GPoint ab = { (1 - t) * src[0].x + t * src[1].x, (1 - t) * src[0].y + t * src[1].y };
    GPoint bc = { (1 - t) * src[1].x + t * src[2].x, (1 - t) * src[1].y + t * src[2].y };
    GPoint abc = { (1 - t) * ab.x + t * bc.x, (1 - t) * ab.y + t * bc.y };

    dst[0] = src[0];
    dst[1] = ab;
    dst[2] = abc;
    dst[3] = bc;
    dst[4] = src[2];
}


void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t) {
    GPoint ab = { (1 - t) * src[0].x + t * src[1].x, (1 - t) * src[0].y + t * src[1].y };
    GPoint bc = { (1 - t) * src[1].x + t * src[2].x, (1 - t) * src[1].y + t * src[2].y };
    GPoint cd = { (1 - t) * src[2].x + t * src[3].x, (1 - t) * src[2].y + t * src[3].y };

    GPoint abc = { (1 - t) * ab.x + t * bc.x, (1 - t) * ab.y + t * bc.y };
    GPoint bcd = { (1 - t) * bc.x + t * cd.x, (1 - t) * bc.y + t * cd.y };

    GPoint abcd = { (1 - t) * abc.x + t * bcd.x, (1 - t) * abc.y + t * bcd.y };

    dst[0] = src[0];
    dst[1] = ab;
    dst[2] = abc;
    dst[3] = abcd;
    dst[4] = bcd;
    dst[5] = cd;
    dst[6] = src[3];
}


void GPathBuilder::addRect(const GRect& rect, GPathDirection dir) {
    GPoint corners[4] = {
        {rect.left, rect.top}, 
        {rect.right, rect.top},  
        {rect.right, rect.bottom}, 
        {rect.left, rect.bottom}  
    };

    int indices[4];
    if (dir == GPathDirection::kCW) {
        indices[0] = 0; indices[1] = 1; indices[2] = 2; indices[3] = 3;
    } else {
        indices[0] = 0; indices[1] = 3; indices[2] = 2; indices[3] = 1;
    }

    this->moveTo(corners[indices[0]]);
    for (int i = 1; i < 4; ++i) {
        this->lineTo(corners[indices[i]]);
    }
}

void GPathBuilder::addPolygon(const GPoint pts[], int count) {
    if (count < 1) return;

    this->moveTo(pts[0]);
    for (int i = 1; i < count; ++i) {
        this->lineTo(pts[i]);
    }
}

void GPathBuilder::addCircle(GPoint center, float radius, GPathDirection dir) {
    const float kappa = 0.551915f; 

    GPoint p0 = { center.x + radius, center.y };
    GPoint p1 = { center.x, center.y + radius }; 
    GPoint p2 = { center.x - radius, center.y };    
    GPoint p3 = { center.x, center.y - radius };          


    GPoint c0 = { center.x + radius, center.y + kappa * radius };     
    GPoint c1 = { center.x + kappa * radius, center.y + radius };     
    GPoint c2 = { center.x - kappa * radius, center.y + radius };     
    GPoint c3 = { center.x - radius, center.y + kappa * radius };     
    GPoint c4 = { center.x - radius, center.y - kappa * radius };      
    GPoint c5 = { center.x - kappa * radius, center.y - radius };     
    GPoint c6 = { center.x + kappa * radius, center.y - radius };     
    GPoint c7 = { center.x + radius, center.y - kappa * radius };   

    if (dir == GPathDirection::kCW) {
        this->moveTo(p0);
        this->cubicTo(c0, c1, p1);
        this->cubicTo(c2, c3, p2);
        this->cubicTo(c4, c5, p3);
        this->cubicTo(c6, c7, p0);
    } else {
        this->moveTo(p0);
        this->cubicTo(c7, c6, p3);
        this->cubicTo(c5, c4, p2);
        this->cubicTo(c3, c2, p1);
        this->cubicTo(c1, c0, p0);
    }
}