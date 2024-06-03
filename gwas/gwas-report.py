from pathlib import Path
from typing import Optional, Tuple

from reportlab.graphics.charts.barcharts import VerticalBarChart  # 图表类
from reportlab.graphics.charts.legends import Legend  # 图例类
from reportlab.graphics.shapes import Drawing  # 绘图工具
from reportlab.lib import colors  # 颜色模块
from reportlab.lib.pagesizes import A4  # 页面的标志尺寸(8.5*inch, 11*inch)
from reportlab.lib.styles import getSampleStyleSheet  # 文本样式
from reportlab.lib.units import cm  # 单位：cm
from reportlab.pdfbase import pdfmetrics  # 注册字体
from reportlab.pdfbase.ttfonts import TTFont  # 字体类
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,  # 报告内容相关类
    SimpleDocTemplate,
    Spacer,
    Table,
)
from reportlab.platypus.frames import Frame

# 注册字体(提前准备好字体文件, 如果同一个文件需要多种字体可以注册多个)
DATA_PATH = Path(__file__).parent / "data"
FONT_PATH = DATA_PATH / "NotoSansSC-Regular.ttf"
LOGO_PATH = DATA_PATH / "logo.png"
pdfmetrics.registerFont(TTFont("NotoSansSC", str(FONT_PATH)))

from functools import partial

import PIL


class Graphs:
    # 绘制标题
    @staticmethod
    def draw_title(title: str):
        # 获取所有样式表
        style = getSampleStyleSheet()
        # 拿到标题样式
        ct = style["Heading1"]
        # 单独设置样式相关属性
        ct.fontName = "NotoSansSC"  # 字体名
        ct.fontSize = 18  # 字体大小
        ct.leading = 50  # 行间距
        ct.alignment = 1  # 居中
        ct.bold = True
        # 创建标题对应的段落，并且返回
        return Paragraph(title, ct)

    # 绘制小标题
    @staticmethod
    def draw_little_title(title: str):
        # 获取所有样式表
        style = getSampleStyleSheet()
        # 拿到标题样式
        ct = style["Normal"]
        # 单独设置样式相关属性
        ct.fontName = "NotoSansSC"  # 字体名
        ct.fontSize = 12  # 字体大小
        ct.leading = 30  # 行间距
        # 创建标题对应的段落，并且返回
        return Paragraph(title, ct)

    # 绘制普通段落内容
    @staticmethod
    def draw_text(text: str):
        # 获取所有样式表
        style = getSampleStyleSheet()
        # 获取普通样式
        ct = style["Normal"]
        ct.fontName = "NotoSansSC"
        ct.fontSize = 10
        ct.wordWrap = "CJK"  # 设置自动换行
        ct.alignment = 0  # 左对齐
        ct.firstLineIndent = 16  # 第一行开头空格
        ct.leading = 25
        return Paragraph(text, ct)

    # 绘制图片
    @staticmethod
    def draw_img(path, width):
        img = Image(path)  # 读取指定路径下的图片
        height = image_auto_height(path, width)
        img.drawWidth = width * cm
        img.drawHeight = height * cm  # 设置图片的宽度
        return img

    @staticmethod
    def draw_table(*args):
        # 列宽度
        col_width = 120
        style = [
            ("FONTNAME", (0, 0), (-1, -1), "NotoSansSC"),  # 字体
            ("FONTSIZE", (0, 0), (-1, 0), 10),  # 第一行的字体大小
            ("FONTSIZE", (0, 1), (-1, -1), 10),  # 第二行到最后一行的字体大小
            ("BACKGROUND", (0, 0), (-1, 0), "#d5dae6"),  # 设置第一行背景颜色
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),  # 第一行水平居中
            ("ALIGN", (0, 1), (-1, -1), "LEFT"),  # 第二行到最后一行左右左对齐
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),  # 所有表格上下居中对齐
            ("TEXTCOLOR", (0, 0), (-1, -1), colors.darkslategray),  # 设置表格内文字颜色
            (
                "GRID",
                (0, 0),
                (-1, -1),
                0.5,
                colors.grey,
            ),  # 设置表格框线为grey色，线宽为0.5
            # ('SPAN', (0, 1), (0, 2)),  # 合并第一列二三行
            # ('SPAN', (0, 3), (0, 4)),  # 合并第一列三四行
            # ('SPAN', (0, 5), (0, 6)),  # 合并第一列五六行
            # ('SPAN', (0, 7), (0, 8)),  # 合并第一列五六行
        ]
        table = Table(args, colWidths=col_width, style=style)
        return table


def image_auto_height(path: Path, width: int) -> int:
    img = PIL.Image.open(path)
    height = int(width * img.height / img.width)
    return height


def image_paths(data_path: Path) -> Tuple[str, str, str]:
    coverage_path = ""
    coverage_hist_path = ""
    genome_region_path = ""
    for img in data_path.glob("*.png"):
        if "coverage.density" in img.name:
            coverage_path = img
        if "coverage.hist" in img.name:
            coverage_hist_path = img
        if "genome_region.pie" in img.name:
            genome_region_path = img
    return str(coverage_path), str(coverage_hist_path), str(genome_region_path)


def header(canvas, doc, content):
    canvas.saveState()
    w, h = content.wrap(doc.width, doc.topMargin)
    content.drawOn(
        canvas, doc.leftMargin, doc.height + doc.bottomMargin + doc.topMargin - h
    )
    canvas.restoreState()


def report(
    species: str,
    genome: str,
    probe_tag: str,
    data_path: Path,
    report_path: Path,
    maf: Optional[float] = 0.05,
    gc_low: Optional[float] = 0.3,
    gc_high: Optional[float] = 0.6,
    match: Optional[str] = "30,1",
) -> None:
    # 生成pdf文件
    doc = SimpleDocTemplate(
        str(report_path),
        pagesize=A4,
        leftMargin=2.2 * cm,
        rightMargin=2.2 * cm,
        topMargin=1.5 * cm,
        bottomMargin=2.5 * cm,
    )
    banner_height = image_auto_height(LOGO_PATH, 150)
    header_content = Image(str(LOGO_PATH), width=150, height=banner_height)

    # 创建内容对应的空列表
    content = list()

    content.append(Spacer(doc.width, banner_height + 1))

    content.append(Graphs.draw_title("芯片设计报告"))

    # 添加表格
    content.append(Graphs.draw_little_title("基本信息"))
    data = [
        ("物种", "基因组", "探针数"),
        (f"{species}", f"{genome}", f"{probe_tag}"),
    ]
    content.append(Graphs.draw_table(*data))
    content.append(Spacer(doc.width, 30))

    # 设计参数
    content.append(Graphs.draw_little_title("设计参数"))
    if maf:
        content.append(
            Graphs.draw_text(f"探针覆盖的位点在群体数据中最小等位基因频率 >= {maf}")
        )
    if gc_low and gc_high:
        content.append(Graphs.draw_text(f"探针GC含量在 {gc_low} - {gc_high} 之间"))
    if match:
        match_len, match_count = match.split(",")
        content.append(
            Graphs.draw_text(
                f"探针{match_len}bp以上长度比对到基因组不超过{match_count}个位置"
            )
        )
    content.append(PageBreak())

    coverage_path, coverage_hist_path, genome_region_path = image_paths(data_path)
    if coverage_path:
        content.append(Graphs.draw_little_title("芯片在基因组的覆盖"))
        content.append(Graphs.draw_img(coverage_path, 18))
        content.append(Graphs.draw_img(coverage_hist_path, 18))
        content.append(PageBreak())

    if genome_region_path:
        content.append(Graphs.draw_little_title("芯片在基因区域的分布"))
        content.append(Graphs.draw_img(genome_region_path, 10))

    doc.build(
        content,
        onFirstPage=partial(header, content=header_content),
        onLaterPages=partial(header, content=header_content),
    )
