from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import typer
from attrs import define, field
from fpdf import FPDF

DATA_PATH = Path(__file__).parent / "data"
FONT_PATH = DATA_PATH / "NotoSansSC-Regular.ttf"
LOGO_PATH = DATA_PATH / "logo.png"


# 创建PDF文档类
class PDF(FPDF):
    def header(self):
        self.image(str(LOGO_PATH), 10, 8, 33)
        self.set_font("SimHei", "", 12)
        self.set_text_color(128, 128, 128)  # 设置文字颜色为灰色
        # Moving cursor to the right:
        self.cell(80)
        # Printing title:
        self.cell(30, 5, "GWAS 分析说明")
        self.set_text_color(0, 0, 0)  # 还原文字颜色为黑色
        self.ln(20)

    def chapter_title(self, title):
        self.set_font("SimHei", "", 12)
        self.cell(0, 10, title, 0, 1, "L")
        self.ln(5)

    def chapter_body(self, body):
        self.set_font("SimHei", "", 12)
        self.multi_cell(0, 10, body)
        self.ln()

    def add_image(self, image_path):
        self.image(image_path, x=10, h=150)
        self.ln(10)

    def add_wide_image(self, image_path):
        self.image(image_path, x=10, w=150)
        self.ln(10)

    def check_page_break(self, h):
        # 如果剩余空间不足以容纳内容，则添加新页面
        if self.get_y() + h > self.page_break_trigger:
            self.add_page()


@define
class GwasReportItems:

    bestK: int = field(factory=int)
    sections: list = field(init=False)

    def __attrs_post_init__(self):
        self.sections = [
            {
                "title": "1. 系统发育树分析",
                "body": "IQ-TREE 2是一款高效的系统发育树构建软件，支持多种替换模型和快速计算。它通过最大似然法（Maximum Likelihood）来推断系统发育树，并使用快速自举分析来评估树节点的可靠性。",
                "image_dir": "系统发育树分析",
                "images": ["snp.circular.png"],
            },
            {
                "title": "2. 主成分分析",
                "body": "PLINK是一款常用的基因组数据分析工具，支持多种遗传统计分析。主成分分析（PCA）通过降维技术，将高维基因型数据投影到低维空间，以揭示样本间的遗传结构和群体分层。",
                "image_dir": "主成分分析",
                "images": ["pca.png"],
            },
            {
                "title": "3. 群体结构分析",
                "body": "ADMIXTURE是一款快速估计种群结构的软件，通过最大似然法估计每个个体的混合比例，揭示个体在不同种群中的遗传成分。",
                "image_dir": "群体结构分析",
                "images": [
                    "bestK.png",
                    f"plot/snp.{self.bestK}.Q.png",
                ],
                "is_wide_image": [False, True],
            },
            {
                "title": "4. 连锁不平衡分析",
                "body": "PopLDdecay是一款专门用于大规模SNP数据的连锁不平衡（LD）分析工具。它通过计算SNP对之间的LD值（如r²），评估基因组范围内的LD衰减模式。",
                "image_dir": "连锁不平衡分析",
                "images": ["LDdecay.png"],
            },
            {
                "title": "5. 性状关联分析",
                "body": "GEMMA是一款用于基因组宽关联分析（GWAS）和混合线性模型分析的软件。它通过混合线性模型（MLM）控制群体结构和亲缘关系，识别与性状相关的遗传变异。",
                "image_dir": "性状关联分析",
                "images": ["manhatton.png", "qq.png"],
            },
        ]


# 生成PDF报告
def generate_report(gwas_results_path: Path, bestK: int):
    pdf = PDF()
    pdf.add_font("SimHei", "", str(FONT_PATH), uni=True)

    section_height = 60
    image_height = 100

    gwasItems = GwasReportItems(bestK=bestK)

    for section in gwasItems.sections:
        title = section["title"]
        body = section["body"]
        image_dir = gwas_results_path / section["image_dir"]
        if not image_dir.exists() or not image_dir.is_dir():
            continue
        images = [str(image_dir / img) for img in section["images"]]

        pdf.add_page()
        pdf.chapter_title(title)
        ##pdf.check_page_break(section_height)  # 假设标题和内容占用60单位高度
        pdf.chapter_body(body)
        for n, image in enumerate(images):
            pdf.check_page_break(image_height)  # 假设每张图片占用100单位高度
            if "is_wide_image" in section and section["is_wide_image"][n]:
                pdf.add_wide_image(image)
            else:
                pdf.add_image(image)

    output_pdf = gwas_results_path / "analysis_report.pdf"
    pdf.output(f"{output_pdf}")


def main(gwas_results_path: Path, bestk: int):
    if not gwas_results_path.exists() or not gwas_results_path.is_dir():
        typer.echo(f"路径 {gwas_results_path} 不存在或不是一个目录", err=True)
        raise typer.Exit(code=1)
    generate_report(gwas_results_path, bestk)


if __name__ == "__main__":
    typer.run(main)
