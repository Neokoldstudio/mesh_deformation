#include "include/Constraints.h"
#include <cstdio>
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

void Constraints::importConstraints(char *filePath)
{
    FILE *file = fopen(filePath, "r");
    if (file == NULL)
    {
        printf("Error opening file\n");
        return;
    }

    char line[128];
    while (fgets(line, sizeof(line), file))
    {
        char *token = strtok(line, " ");
        if (strcmp(token, "a") == 0)
        {
            token = strtok(NULL, " ");
            addAnchorIndex(atoi(token));
        }
        else if (strcmp(token, "h") == 0)
        {
            token = strtok(NULL, " ");
            addHandleIndex(atoi(token));
        }
    }

    fclose(file);
}

void Constraints::exportConstraints(char *filePath)
{
    FILE *file = fopen(filePath, "w");
    if (file == NULL)
    {
        file = fopen(filePath, "w+");
        if (file == NULL)
        {
            printf("Error opening file\n");
            return;
        }
    }

    for (int i = 0; i < anchorIndices.size(); ++i)
    {
        fprintf(file, "a %d\n", anchorIndices[i]);
    }

    for (int i = 0; i < handleIndices.size(); ++i)
    {
        fprintf(file, "h %d\n", handleIndices[i]);
    }

    fclose(file);
}

void Constraints::setHandleTransformation(const Eigen::Affine3d &transformation)
{
    handleTransformation = transformation;
}

Eigen::Affine3d Constraints::getHandleTransformation() const
{
    return handleTransformation;
}