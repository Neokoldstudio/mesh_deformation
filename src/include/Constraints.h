#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <vector>
#include <Eigen>

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

    // Get anchor at index
    int getAnchorIndex(int index) const;

    // Get handle at index
    int getHandleIndex(int index) const;

    int getAnchorSize() const;

    int getHandleSize() const;

    // Clear all anchor indices
    void clearAnchorIndices();

    // Clear all handle indices
    void clearHandleIndices();

    // Set the transformation for all handles
    void setHandleTransformation(const Eigen::Affine3d &transformation);

    // Get the transformation for all handles
    Eigen::Affine3d getHandleTransformation() const;

    void importConstraints(char *filePath);

    void exportConstraints(char *filePath);

private:
    std::vector<int> anchorIndices;
    std::vector<int> handleIndices;
    Eigen::Affine3d handleTransformation;
};

#endif // CONSTRAINTS_H
