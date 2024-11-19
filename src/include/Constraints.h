#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <vector>
#include <Eigen/Dense>

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

    void importConstraints(const char *filePath);

    void exportConstraints(const char *filePath);

    void updateTransformation(const Eigen::Matrix4f &T);

    Eigen::Matrix4f getTransformation() const;

private:
    std::vector<int> anchorIndices;
    std::vector<int> handleIndices;
    Eigen::Matrix4f transformation = Eigen::Matrix4f::Identity(); // Store a single transformation matrix
};

#endif // CONSTRAINTS_H
