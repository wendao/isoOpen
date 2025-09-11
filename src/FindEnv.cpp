#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include "Envelope.h"

/**
 * 从一张 MS1 质谱图中提取 isotope envelope 模式。
 *
 * @param envelopes  用于存放识别出的所有 Envelope 对象指针，函数开始时会清空原有内容
 * @param mz          mz 值数组（升序）
 * @param intensity   对应的强度数组
 * @param lenMin      envelope 最小长度（含），默认 6
 * @param lenMax      envelope 最大长度（含），默认 14
 */
void FindEnvelope(std::vector<Envelope*>& envelopes,
                  const std::vector<double>& mz,
                  const std::vector<double>& intensity,
                  int lenMin = 6,
                  int lenMax = 14)
{
    envelopes.clear();
    size_t n = mz.size();
    if (n == 0 || intensity.size() != n) return;

    // 13C–12C 质量差常数（Da）
    const double C13C12 = 1.003354835;
    // 匹配容差，可根据仪器分辨率调节
    const double tol = 0.01;

    // 标记峰是否已归入某个 envelope
    std::vector<bool> used(n, false);
    // 电荷优先级：高电荷先扫描
    const std::vector<int> charges = {7, 6, 5, 4, 3, 2, 1};

    for (int charge : charges) {
        double deltaM = C13C12 / charge;

        for (size_t i = 0; i < n; ++i) {
            if (used[i]) continue;

            // 如果从这里出发连最小长度都匹配不到，就跳过
            if (mz[i] + (lenMin - 1) * deltaM > mz.back() + tol)
                break;

            std::vector<size_t> idxs;
            idxs.reserve(lenMax);
            idxs.push_back(i);

            double target = mz[i] + deltaM;
            while ((int)idxs.size() < lenMax) {
                // 在 mz[idxs.back()+1..end) 中找第一个 >= target - tol
                auto lb = std::lower_bound(mz.begin() + idxs.back() + 1, mz.end(), target - tol);
                if (lb == mz.end()) break;

                size_t j = lb - mz.begin();
                if (!used[j] && std::abs(mz[j] - target) <= tol) {
                    idxs.push_back(j);
                    target += deltaM;
                } else {
                    break;
                }
            }

            if ((int)idxs.size() >= lenMin) {
                // 构造一个新 envelope，并标记这些峰已用
                Envelope* e = new Envelope(charge);
                for (size_t idx : idxs) {
                    e->AddPeak(mz[idx], intensity[idx]);
                    used[idx] = true;
                }
                envelopes.push_back(e);
            }
        }
    }
}

int main(int argc, char* argv[]) {
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::ifstream;
    using std::istringstream;
    using std::vector;
    using std::string;

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <ms1_data.txt>\n";
        return EXIT_FAILURE;
    }

    ifstream infile(argv[1]);
    if (!infile) {
        cerr << "Error: cannot open file " << argv[1] << "\n";
        return EXIT_FAILURE;
    }

    vector<double> mzs;
    vector<double> ints;
    string line;
    int spectrum_idx = 0;
    vector<Envelope*> envelopes;

    auto process_spectrum = [&]() {
        if (mzs.empty()) return;
        ++spectrum_idx;

        // 释放上一谱的 Envelope 对象
        for (auto *e : envelopes) delete e;
        envelopes.clear();

        cout << "Spectrum " << spectrum_idx 
             << ": points=" << mzs.size() << endl;

        FindEnvelope(envelopes, mzs, ints, 6, 15);

        // 这里可以遍历 envelopes，打印或者处理每个 envelope
        for (auto *e : envelopes) {
            cout << "  Charge=" << e->GetCharge()
                 << "  Length=" << e->GetLength()
                 << "  Last m/z=" << e->GetLastMz()
                 << endl;
        }

        mzs.clear();
        ints.clear();
    };

    while (getline(infile, line)) {
        if (line.empty()) {
            continue;
        }
        // 如果是新谱开始标志，先处理上一谱
        if (line[0] == 'S') {
            process_spectrum();
            continue;
        }
        // 跳过所有以 'I' 开头的注释行
        if (line[0] == 'I') {
            continue;
        }
        // 否则应为数据行：四列数值，取前两列 mz 和 intensity
        istringstream iss(line);
        double mz, intensity;
        // 读两列，忽略后面两列
        if (iss >> mz >> intensity) {
            mzs.push_back(mz);
            ints.push_back(intensity);
        }
    }
    // 文件末尾，处理最后一谱
    process_spectrum();

    return EXIT_SUCCESS;
}

