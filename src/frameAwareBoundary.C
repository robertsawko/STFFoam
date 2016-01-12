#include "frameAwareBoundary.H"

frameAwareBoundary::frameAwareBoundary() : velocity_(0, 0, 0) {}

frameAwareBoundary::frameAwareBoundary(const frameAwareBoundary &b)
    : velocity_(b.velocity_) {}

vector frameAwareBoundary::currentVelocity() const { return velocity_; }

void frameAwareBoundary::correct(const vector &V) { velocity_ = V; }
