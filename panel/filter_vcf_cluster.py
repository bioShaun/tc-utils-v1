import gzip
import logging
import shutil
import subprocess
from collections import deque
from pathlib import Path

import typer

# 配置日志
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)

app = typer.Typer()


def find_bgzip(bgzip_path: Path = None) -> Path:
    """
    查找 bgzip 可执行文件，优先使用指定路径，否则在系统 PATH 中搜索。
    """
    if bgzip_path:
        exe = bgzip_path / "bgzip" if bgzip_path.is_dir() else bgzip_path
        if exe.exists():
            return exe
        raise FileNotFoundError(f"未在指定位置找到 bgzip: {exe}")
    system_path = shutil.which("bgzip")
    if system_path:
        return Path(system_path)
    raise FileNotFoundError("未在系统 PATH 中找到 bgzip 可执行文件")


def open_vcf(path: Path):
    """
    打开 VCF 文件，自动处理 .vcf 和 .vcf.gz 格式。
    """
    if path.suffix in [".gz", ".bgz"]:
        return gzip.open(path, "rt")
    return open(path, "r")


def write_header(infile, outfile_handle):
    """
    将所有以 '#' 开头的头部行复制到输出文件，并返回第一条非头部变异记录。
    跳过空行，确保正确保留头部信息。
    如果没有找到数据行，则返回 None。
    """
    for line in infile:
        # 跳过空行
        if not line.strip():
            continue
        # 写入所有头部行
        if line.startswith("#"):
            outfile_handle.write(line)
            continue
        # 返回第一条变异记录，用于初始化缓冲区
        return line
    return None


def filter_stream(
    infile_path: Path, outfile_path: Path, window_size: int, snp_count: int
):
    """
    流式处理 VCF 文件，使用滑动窗口方法删除聚簇 SNP。

    infile_path: 输入 VCF 文件路径（.vcf 或 .vcf.gz）
    outfile_path: 输出过滤后 VCF 文件路径（未压缩）
    window_size: 聚簇检测窗口大小（bp）
    snp_count: 聚簇中 SNP 数量减 1，用于定义缓冲区大小
    """
    buffer = deque(maxlen=snp_count + 1)
    with open_vcf(infile_path) as infile, open(outfile_path, "w") as out:
        # 写入头部并获取第一条变异记录
        first_line = write_header(infile, out)
        if not first_line:
            logger.warning("输入 VCF 文件中未找到变异记录。")
            return

        # 辅助函数：解析变异行，只提取 CHROM 和 POS
        def parse_variant(line):
            if not line.strip():
                return None
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                return None
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                return None
            return {"line": line, "chrom": chrom, "pos": pos, "remove": False}

        # 初始化缓冲区
        first_var = parse_variant(first_line)
        if first_var:
            buffer.append(first_var)

        # 逐行处理
        for line in infile:
            # 跳过头部或空行
            if not line.strip() or line.startswith("#"):
                continue
            var = parse_variant(line)
            if not var:
                continue
            buffer.append(var)
            # 如果缓冲区填满，则检测并输出最旧记录
            if len(buffer) == buffer.maxlen:
                oldest, newest = buffer[0], buffer[-1]
                if oldest["chrom"] == newest["chrom"] and (
                    newest["pos"] - oldest["pos"] <= window_size
                ):
                    # 标记所有缓冲区记录为待删除
                    for entry in buffer:
                        entry["remove"] = True
                # 弹出最旧记录，未标记则写入
                to_write = buffer.popleft()
                if not to_write["remove"]:
                    out.write(to_write["line"])
        # 清空剩余缓冲区
        while buffer:
            entry = buffer.popleft()
            if not entry["remove"]:
                out.write(entry["line"])


@app.command()
def main(
    vcf_file: Path = typer.Argument(
        ..., exists=True, help="输入 VCF 文件路径（.vcf 或 .vcf.gz）"
    ),
    output_file: Path = typer.Argument(
        ..., help="输出过滤后 VCF 文件路径（未压缩，后续会使用 bgzip 压缩）"
    ),
    window_size: int = typer.Option(60, help="聚簇检测窗口大小（bp）"),
    snp_count: int = typer.Option(
        2,
        min=1,
        help="定义簇的大小。一个簇由 (snp_count + 1) 个 SNP 组成。"
        "如果这些 SNP 在 window_size 内，则所有 (snp_count + 1) 个 SNP 都将被移除。"
        "例如，默认值 2 表示 3 个 SNP 构成一个簇。",
    ),
    bgzip_path: Path = typer.Option(None, help="bgzip 可执行文件或其所在目录路径"),
):
    """
    过滤 VCF 文件中的聚簇 SNP，并使用 bgzip 压缩输出结果。
    """
    logger.info(f"开始过滤 VCF: {vcf_file} -> {output_file}")
    try:
        filter_stream(vcf_file, output_file, window_size, snp_count)
        bgzip_exe = find_bgzip(bgzip_path)
        subprocess.run([str(bgzip_exe), "-f", str(output_file)], check=True)
        logger.info(f"成功压缩输出文件为 {output_file}.gz")
    except Exception as e:
        logger.error(f"过滤或压缩过程中出现错误: {e}")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
