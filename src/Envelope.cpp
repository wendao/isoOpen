#include "Envelope.h"
#include <stdexcept>

/**
 * @brief Construct a new Envelope with a given charge state.
 * @param charge The charge of the envelope (e.g., 1, 2, ..., 6).
 */
Envelope::Envelope(int charge)
    : charge_(charge)
{
    // 初始化时不需要额外操作
}

/**
 * @brief Add a peak (m/z, intensity) to this envelope.
 * @param mz The m/z value of the peak.
 * @param intensity The intensity of the peak.
 */
void Envelope::AddPeak(double mz, double intensity) {
    mzs_.push_back(mz);
    intensities_.push_back(intensity);
}

/**
 * @brief Get the charge state of this envelope.
 * @return int The charge state.
 */
int Envelope::GetCharge() const {
    return charge_;
}

/**
 * @brief Get the number of peaks in this envelope.
 * @return std::size_t The number of peaks.
 */
std::size_t Envelope::GetLength() const {
    return mzs_.size();
}

/**
 * @brief Get the m/z value of the last (highest m/z) peak.
 * @return double The m/z of the last peak.
 */
double Envelope::GetLastMz() const {
    if (mzs_.empty()) {
        throw std::out_of_range("Envelope contains no peaks");
    }
    return mzs_.back();
}

/**
 * @brief Get the m/z value of a peak at a specific index.
 * @param index The zero-based index of the peak.
 * @return double The m/z value.
 */
double Envelope::GetMz(std::size_t index) const {
    if (index >= mzs_.size()) {
        throw std::out_of_range("Index out of range in GetMz");
    }
    return mzs_[index];
}

/**
 * @brief Get the intensity of a peak at a specific index.
 * @param index The zero-based index of the peak.
 * @return double The intensity value.
 */
double Envelope::GetIntensity(std::size_t index) const {
    if (index >= intensities_.size()) {
        throw std::out_of_range("Index out of range in GetIntensity");
    }
    return intensities_[index];
}

