#include "vector.H"

using Foam::vector;

class frameAwareBoundary {
  private:
    vector velocity_;

  public:

    frameAwareBoundary();
    frameAwareBoundary(const frameAwareBoundary&);

    void correct(const vector &);
    vector currentVelocity() const;
};
