#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <vector>

class Constraints
{
public:
    // Constructor
    Constraints() = default;

    // Destructor
    ~Constraints() = default;

    // Add an anchor index
    void addAnchorIndex(int index);

    // Add a handle index
    void addHandleIndex(int index);

    // Get all anchor indices
    const std::vector<int> &getAnchorIndices() const;

    // Get all handle indices
    const std::vector<int> &getHandleIndices() const;

    int getAnchorIndex(int index) const;

    int getHandleIndex(int index) const;

    int getAnchorSize() const;

    int getHandleSize() const;

    // Clear all anchor indices
    void clearAnchorIndices();

    // Clear all handle indices
    void clearHandleIndices();

    void importConstraints(char *filePath);

    void exportConstraints(char *filePath);

private:
    std::vector<int> anchorIndices;
    std::vector<int> handleIndices;
};

#endif // CONSTRAINTS_H
