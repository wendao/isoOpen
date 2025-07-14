#ifndef ENVELOPE_H
#define ENVELOPE_H

#include <vector>
#include <cstddef>

/**
 * @brief Represents an isotope envelope detected in an MS1 spectrum.
 * Stores the charge state and a list of m/z and intensity values for the peaks.
 */
class Envelope {
public:
    /**
     * @brief Construct a new Envelope with a given charge state.
     * @param charge The charge of the envelope (e.g., 1, 2, ..., 6).
     */
    explicit Envelope(int charge);

    /**
     * @brief Add a peak (m/z, intensity) to this envelope.
     * @param mz The m/z value of the peak.
     * @param intensity The intensity of the peak.
     */
    void AddPeak(double mz, double intensity);

    /**
     * @brief Get the charge state of this envelope.
     * @return int The charge state.
     */
    int GetCharge() const;

    /**
     * @brief Get the number of peaks in this envelope.
     * @return std::size_t The number of peaks.
     */
    std::size_t GetLength() const;

    /**
     * @brief Get the m/z value of the last (highest m/z) peak.
     * @return double The m/z of the last peak.
     */
    double GetLastMz() const;

    /**
     * @brief Get the m/z value of a peak at a specific index.
     * @param index The zero-based index of the peak.
     * @return double The m/z value.
     */
    double GetMz(std::size_t index) const;

    /**
     * @brief Get the intensity of a peak at a specific index.
     * @param index The zero-based index of the peak.
     * @return double The intensity value.
     */
    double GetIntensity(std::size_t index) const;

private:
    int charge_;                               ///< Envelope charge state
    std::vector<double> mzs_;                  ///< m/z values of detected peaks
    std::vector<double> intensities_;          ///< Intensities of detected peaks
};

#endif // ENVELOPE_H

