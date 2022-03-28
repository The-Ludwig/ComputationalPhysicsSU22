#pragma once
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <sstream>

class NumpySaver
{
    typedef double T;
    const char *file_name;
    bool warn = false;
    std::vector<std::vector<T> *> nums;
    std::stringstream header;
    int size = -1;
    unsigned int precision = 5;

    void save_it()
    {

        assert(nums.size() != 0);
        std::ofstream file;
        file.precision(5);
        file << std::fixed;
        file.open(file_name, std::ofstream::out | std::ofstream::trunc);

        file << "#" << header.str() << "\n";

        for (unsigned int i = 0; i < size; i++)
        {
            for (int j = 0; j < nums.size(); j++)
            {

                file << (*nums[j])[i] << " ";
            }

            file << "\n";
        }

        file.close();
    }

public:
    enum signal
    {
        save,
        warning_off,
        warning_on
    };

    NumpySaver(const char *file_name)
    {
        this->file_name = file_name;
    }

    ~NumpySaver()
    {
        if (nums.size() > 0)
            save_it();
        for (int i = 0; i < nums.size(); i++)
            delete nums[i];
    }

    void setPrecision(unsigned int precision)
    {
        this->precision = precision;
    }

    NumpySaver &operator<<(signal sig)
    {
        switch (sig)
        {
        case signal::save:
            if (nums.size() > 0)
            {
                save_it();
            }
            else
            {
                std::cerr << "Trying to save, but no data was given (NumpySaver)";
            }
            break;
        case signal::warning_off:
            warn = false;
        case signal::warning_on:
            warn = true;
        default:
            break;
        }
        return *this;
    }

    template <typename Derived>
    NumpySaver &operator<<(const Eigen::MatrixBase<Derived> &mat)
    {
        if (size == -1)
            size = mat.size();
        if (size != mat.size())
            std::cerr << "Error saving numpy comform file: all nums/vectors must have the same length. Weird errors will occur!";

        std::vector<T> *num = new std::vector<T>(mat.size());
        Eigen::VectorXd::Map(&((*num)[0]), mat.size()) = mat.reshaped();
        nums.push_back(num);
        return *this;
    }

    NumpySaver &operator<<(const std::vector<T> &mat)
    {
        if (size == -1)
            size = mat.size();
        if (size != mat.size())
            std::cerr << "Error saving numpy comform file: all nums/vectors must have the same length. Weird errors will occur!";
        std::vector<T> *num = new std::vector<T>(mat);
        nums.push_back(num);
        return *this;
    }

    NumpySaver &operator<<(const T &val)
    {
        if (size == -1)
            size = 1;
        if (size != 1)
            std::cerr << "Error saving numpy comform file: all nums/vectors must have the same length. Weird errors will occur!";
        std::vector<T> *num = new std::vector<T>(1);
        num->push_back(val);
        nums.push_back(num);
        return *this;
    }

    NumpySaver &operator<<(const char *name)
    {
        header << name << " ";
        return *this;
    }

    NumpySaver &operator<<(std::string &name)
    {
        header << name << " ";
        return *this;
    }
};