from pathlib import Path

import pandas as pd

from gtf.restore_split_bed import (
    main,
)  # 请将 `your_script` 替换为你的 Python 脚本文件名（不带 `.py`）


def test_restore_bed(tmp_path):
    # 构造输入文件路径
    bed_file = tmp_path / "test.bed"
    split_file = tmp_path / "test.split.bed"
    output_file = tmp_path / "test.restore.bed"

    # 写入 test.bed 内容（模拟拆分后的子染色体区域）
    bed_content = "split1\t100\t200\nsplit2\t0\t50\n"
    bed_file.write_text(bed_content)

    # 写入 test.split.bed 内容（映射关系）
    split_content = "chr1\t0\t500\tsplit1\nchr1\t500\t1000\tsplit2\n"
    split_file.write_text(split_content)

    # 调用待测试的 main 函数
    main(bed=bed_file, split_bed=split_file, restore_bed=output_file)

    # 读取输出文件内容
    result_df = pd.read_csv(output_file, sep="\t", header=None)
    result_df.columns = ["chrom", "start", "end"]

    # 构造期望结果
    expected = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 500],
            "end": [200, 550],
        }
    )

    # 比较结果
    pd.testing.assert_frame_equal(result_df, expected)

    print("✅ 测试通过，输出正确")
