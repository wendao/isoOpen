from typing import Iterator, Tuple
import numpy as np
from pyteomics import mzml, mzxml

def ms1_spectra(path: str) -> Iterator[Tuple[np.ndarray, np.ndarray]]:
    """
    逐谱返回 (mz, intensity)；仅包含 MS1。
    支持 .mzML(.gz) 与 .mzXML(.gz)
    """
    # 根据扩展名选择读取器
    lower = path.lower()
    if lower.endswith((".mzml", ".mzml.gz")):
        opener = mzml.MzML
        level_key = "ms level"
    elif lower.endswith((".mzxml", ".mzxml.gz")):
        opener = mzxml.MzXML
        level_key = "msLevel"
    else:
        raise ValueError("只支持 mzML/mzXML（可.gz）")

    # 流式读取，内存友好
    with opener(path) as reader:
        for spec in reader:
            # 仅保留 MS1
            if int(spec.get(level_key, 0)) != 1:
                continue

            # 提取并转成连续的 float64 数组
            mz = np.asarray(spec["m/z array"], dtype=np.float64)
            intensity = np.asarray(spec["intensity array"], dtype=np.float64)

            if mz.size == 0 or intensity.size == 0:
                continue

            yield mz, intensity

