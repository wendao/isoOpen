# Envelope 同位素峰包络检测算法

## 概述

`Envelope` 类和相关算法用于从质谱 MS1 数据中检测同位素峰包络 (Isotope Envelope)。该算法基于同位素峰之间的固定质量间隔来识别具有特定电荷态的峰系列。

## 核心原理

### 同位素峰间隔

对于带电荷 $z$ 的肽段离子，相邻同位素峰（${}^{13}\mathrm{C}$ 替换 ${}^{12}\mathrm{C}$）之间的 m/z 间隔为：

$$\Delta m/z = \frac{m_{^{13}\mathrm{C}} - m_{^{12}\mathrm{C}}}{z} = \frac{1.003354835}{z}$$

其中：
- ${}^{13}\mathrm{C}$ 与 ${}^{12}\mathrm{C}$ 的质量差 = 1.003354835 Da
- $z$ 为电荷态 (1, 2, 3, ..., 7)

## 数据结构

### Envelope 类 (Envelope.h/cpp)

```cpp
class Envelope {
private:
    int charge_;                    // 电荷态
    std::vector<double> mzs_;      // m/z 值列表
    std::vector<double> intensities_; // 峰强度列表

public:
    explicit Envelope(int charge);
    // Move semantics
    Envelope(Envelope&&) noexcept = default;
    Envelope& operator=(Envelope&&) noexcept = default;

    void AddPeak(double mz, double intensity);
    void Reserve(std::size_t capacity);
    [[nodiscard]] int GetCharge() const;
    [[nodiscard]] std::size_t GetLength() const;
    [[nodiscard]] double GetLastMz() const;
    [[nodiscard]] double GetMz(std::size_t index) const;
    [[nodiscard]] double GetIntensity(std::size_t index) const;
    [[nodiscard]] const std::vector<double>& GetMzs() const;
    [[nodiscard]] const std::vector<double>& GetIntensities() const;
};
```

`Envelope` 是一个轻量级数据结构，使用 `std::unique_ptr` 管理生命周期。

## 检测算法

### FindEnv.cpp 主要函数

#### 1. CalculateSpacingError()

计算峰包络的间隔误差（单位：ppm），使用**中位数**（抗异常值优于均值）。

```cpp
double CalculateSpacingError(const std::vector<double>& mz, int charge)
```

- **输入**: m/z 列表、电荷态
- **输出**: 中位数 ppm 误差
- **公式**:
  $$\text{error}_{ppm} = \frac{|\Delta m/z_{actual} - \Delta m/z_{expected}|}{\Delta m/z_{expected}} \times 10^6$$

其中 $\Delta m/z_{expected} = 1.003354835 / charge$

#### 2. SplitAndValidate()

递归拆分超长包络。在最低强度内部峰处拆分，左右子包络递归验证。

```cpp
static void SplitAndValidate(const std::vector<size_t>& idxs,
                             const std::vector<double>& intensity,
                             int lenMin, int lenMax,
                             std::vector<std::vector<size_t>>& result)
```

- 若长度 < lenMin：丢弃
- 若长度 <= lenMax：保留
- 若长度 > lenMax：在最低强度内部峰处拆分，递归处理子包络

#### 3. FindAllCandidates()

枚举所有可能的候选包络，不做贪心选择。

```cpp
void FindAllCandidates(std::vector<CandidateEnvelope>& candidates,
                       const std::vector<double>& mz,
                       const std::vector<double>& intensity,
                       int lenMin, int lenMax,
                       double minIntensity,
                       double minTotalIntensity,
                       double maxSpacingPpm)
```

**参数**:
- `charges = {6, 5, 4, 3, 2, 1}` - 默认尝试的电荷态
- `lenMin = 6` - 最小峰数
- `lenMax = 15` - 最大峰数（超长包络会被拆分）
- `ppmTol = 12.0` - 通用默认动态 ppm 容差
- `minIntensity = 0.0` - 可选的强度预过滤阈值
- `minTotalIntensity` / `maxSpacingPpm` - 进入冲突图前的候选质量门槛

UPS2 的 `--profile ups2-gt` 使用 `ppmTol=6`、`lenMin=3`、
`minTotalIntensity=3500`、`maxSpacingPpm=15000`。绝对强度门槛与仪器量纲相关，
不应直接搬到其他数据集。

**算法流程**:

```
若 minIntensity > 0: 过滤低强度峰，建立原始索引映射

对于每个电荷态 z:
    deltaM = 1.003354835 / z

    对于每个起始峰 i:
        1. 初始化包络，包含种子峰 i

        Phase A: 向后扩展（低 m/z 方向）
            target = mz[i] - deltaM
            从 i-1 开始向低 m/z 扫描
            在动态 ppm 容差内选最佳匹配峰（最小 m/z 偏差，tie-break 强度）
            继续扩展直到无匹配峰

        Phase B: 向前扩展（高 m/z 方向）
            target = mz[i] + deltaM
            从 i+1 开始向高 m/z 扫描
            在动态 ppm 容差内选最佳匹配峰（最小 m/z 偏差，tie-break 强度）
            继续扩展直到无匹配峰

        2. 若使用预过滤，将过滤索引映射回原始索引

        3. SplitAndValidate: 若长度 > lenMax，递归拆分为有效子包络

        4. 计算总强度和中位间距误差；不满足质量门槛的候选在冲突选择前丢弃

        5. 对保留候选计算评分，添加到候选列表
```

**评分公式**:

$$\text{score} = \log(\sum I + 1) \times 2.0 - \text{spacingError} \times 0.001 + \log_2(n) \times 0.5$$

其中：
- $\sum I$ = 所有峰强度之和
- spacingError = 中位数间隔误差 (ppm)
- $n$ = 峰数量

#### 4. ImprovedSelection()

冲突图 + 贪心 + 局部搜索优化的非贪婪选择算法。

```cpp
void ImprovedSelection(std::vector<std::unique_ptr<Envelope>>& envelopes,
                      std::vector<CandidateEnvelope>& candidates)
```

**执行步骤**:

1. **按评分降序排列**所有候选
2. **构建冲突图**: 重叠 ≥50% 的候选视为冲突（m/z 范围预检查加速）
3. **贪心初选**: 按分数降序依次选择不冲突的候选
4. **局部搜索优化** (最多 500 次迭代):
   - 遍历未选中候选，查找与之冲突且分数更低的已选候选
   - 尝试替换，若总得分提升则接受
5. **收集结果**: 选中的 unique_ptr 移动到输出向量，未选中的自动析构

#### 5. FindEnvelope()

主入口函数，执行完整的包络检测流程。

```cpp
void FindEnvelope(std::vector<std::unique_ptr<Envelope>>& envelopes,
                  const std::vector<double>& mz,
                  const std::vector<double>& intensity,
                  int lenMin = 6,
                  int lenMax = 15,
                  double minIntensity = 0.0,
                  double minTotalIntensity = 0.0,
                  double maxSpacingPpm = infinity)
```

## 输入格式

程序读取 MS1 文本格式文件:

```
S                  # 谱图开始标记
I rt=1.23          # 信息行（跳过）
H ...              # Header 行（跳过）
mz  intensity  ... # 数据行（至少两列）
mz  intensity  ...
S                  # 下一个谱图开始
...
```

## 输出示例

```
Spectrum 1: points=500
  Charge=3  Length=8  Last m/z=523.41
  Charge=2  Length=10  Last m/z=785.23
```

## 算法特点

### 优点

1. **双向扩展**: 从种子峰同时向低 m/z 和高 m/z 方向生长，避免遗漏单同位素峰
2. **非贪心候选查找**: 先找出所有可能的包络，再做选择，避免局部最优
3. **多电荷态支持**: 默认同时尝试 6 种电荷态 (6→1)
4. **择优匹配**: 在容差窗口内选 m/z 偏差最小的峰，偏差相同时选强度最高者
5. **动态 ppm 容差**: 容差随 m/z 缩放（通用默认 12 ppm），比固定 Da 更合理
6. **超长包络拆分**: 超过 lenMax 的包络递归拆分为有效子包络，而非直接丢弃
7. **强度预过滤**: 可选的 minIntensity 参数过滤噪声峰
8. **评分机制**: 综合考虑强度、中位数间隔精度、长度三个因素
9. **冲突图 + 局部搜索**: 比纯贪心更接近全局最优
10. **候选前置过滤**: 强度和间距门槛在冲突图选择前应用，避免低质量候选阻塞更好的重叠候选
11. **内存安全**: 使用 `std::unique_ptr` 管理 Envelope 生命周期，无内存泄漏

### 与 Python 版的对比

| 特性 | C++ 版 | Python 版 |
|------|--------|-----------|
| 扩展方向 | 双向 | 双向 |
| 选择策略 | 冲突图 + 贪心 + 局部搜索 | 贪婪（强度排序）+ 去重 |
| 峰匹配 | 择优（最小偏差 + 强度 tie-break） | 择优（最小偏差 + 强度 tie-break） |
| 容差 | 12 ppm 动态（运行时可改） | max(0.001 Da, 10 ppm) |
| 超长包络 | 递归拆分 | 递归拆分 |
| 强度预过滤 | minIntensity 参数 | min_intensity 参数 |
| 电荷顺序 | {6,5,4,3,2,1} | {7,6,5,4,3,2} |
| 间距误差 | 中位数 | 中位数 |
| 内存管理 | unique_ptr (RAII) | GC |
| 数据格式 | .ms1 文本 | .mzML/.mzXML |

## 编译与运行

```bash
# 编译（Linux GCC / macOS Homebrew libomp 自动适配）
bash comp.sh

# 通用默认运行
./find_envelope ms1_data.txt

# UPS2 GT 优化配置
./find_envelope ms1_data.txt --profile ups2-gt
```

## 相关文件

- `Envelope.h` - Envelope 类头文件
- `Envelope.cpp` - Envelope 类实现
- `FindEnv.cpp` - 包络检测算法实现
- `test.ms1`, `fast_filter.ms1` - 测试数据
