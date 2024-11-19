#include "include/Constraints.h"
#include <nlohmann/json.hpp>
#include <cstdio>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdexcept>

void Constraints::addAnchorIndex(int index)
{
    anchorIndices.push_back(index);
}

void Constraints::addHandleIndex(int index)
{
    handleIndices.push_back(index);
}

const std::vector<int> &Constraints::getAnchorIndices() const
{
    return anchorIndices;
}

const std::vector<int> &Constraints::getHandleIndices() const
{
    return handleIndices;
}

int Constraints::getAnchorIndex(int index) const
{
    if (index >= 0 && index < anchorIndices.size())
    {
        return anchorIndices[index];
    }
    throw std::out_of_range("Index out of range");
}

int Constraints::getHandleIndex(int index) const
{
    if (index >= 0 && index < handleIndices.size())
    {
        return handleIndices[index];
    }
    throw std::out_of_range("Index out of range");
}

int Constraints::getAnchorSize() const
{
    return anchorIndices.size();
}

int Constraints::getHandleSize() const
{
    return handleIndices.size();
}

void Constraints::clearAnchorIndices()
{
    anchorIndices.clear();
}

void Constraints::clearHandleIndices()
{
    handleIndices.clear();
}

void Constraints::updateTransformation(const Eigen::Matrix4f &T)
{
    transformation = T;
}

Eigen::Matrix4f Constraints::getTransformation() const
{
    return transformation;
}

void Constraints::importConstraints(const char *filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open())
    {
        printf("Error opening file");
        return;
    }

    nlohmann::json j;
    file >> j;

    anchorIndices = j["anchors"].get<std::vector<int>>();
    handleIndices = j["handles"].get<std::vector<int>>();
    std::vector<float> transform = j["transformation"].get<std::vector<float>>();
    transformation = Eigen::Map<Eigen::Matrix4f>(transform.data());

    file.close();
}

void Constraints::exportConstraints(const char *filePath)
{
    std::ofstream file(filePath, std::ios::out | std::ios::trunc);
    if (!file.is_open())
    {
        printf("Error opening file");
        return;
    }

    nlohmann::json j;
    j["anchors"] = anchorIndices;
    j["handles"] = handleIndices;
    std::vector<float> transform(transformation.data(), transformation.data() + transformation.size());
    j["transformation"] = transform;

    file << j.dump(4);
    file.close();
}