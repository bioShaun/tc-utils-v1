from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import typer

from fpdf import FPDF
from attrs import define, field

DATA_PATH = Path(__file__).parent / "data"
FONT_PATH = DATA_PATH / "NotoSansSC-Regular.ttf"
LOGO_PATH = DATA_PATH / "logo.png"




# 创建PDF文档类
class PDF(FPDF):
    def header(self):
        self.image(LOGO_PATH, 10, 8, 33)
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
class GwasReportPlot:

    gwas_result_path: Path = field(factory=Path)
    tree_png: Path = field(init=False)

    def __attrs_post_init__(self):
        self.tree_png = gwas_result_path / 


# 生成系统发育树分析图（示例）
def generate_phylogenetic_tree_plot():
    # 假设这里有生成系统发育树的代码
    plt.figure(figsize=(6, 4))
    plt.title("系统发育树")
    plt.text(0.5, 0.5, "系统发育树图示例", ha="center", va="center")
    plt.savefig("tree_visualization.png")
    plt.close()


# 生成PCA分析图（示例）
def generate_pca_plot():
    # 假设这里有PCA分析的代码
    plt.figure(figsize=(6, 4))
    plt.title("PCA分析结果")
    plt.scatter(np.random.rand(50), np.random.rand(50))
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.savefig("pca_plot.png")
    plt.close()


# 生成ADMIXTURE分析图（示例）
def generate_admixture_plot():
    # 假设这里有ADMIXTURE分析的代码
    plt.figure(figsize=(6, 4))
    plt.title("ADMIXTURE分析结果")
    plt.bar(range(10), np.random.rand(10))
    plt.xlabel("Individuals")
    plt.ylabel("Ancestry Proportion")
    plt.savefig("admixture_plot.png")
    plt.close()


# 生成LD衰减图（示例）
def generate_ld_decay_plot():
    # 假设这里有LD分析的代码
    plt.figure(figsize=(6, 4))
    plt.title("LD衰减曲线")
    plt.plot(np.linspace(0, 100, 100), np.exp(-np.linspace(0, 100, 100) / 20))
    plt.xlabel("Distance (kb)")
    plt.ylabel("LD (r^2)")
    plt.savefig("ld_decay_plot.png")
    plt.close()


# 生成GWAS分析图（示例）
def generate_gwas_plot():
    # 假设这里有GWAS分析的代码
    plt.figure(figsize=(6, 4))
    plt.title("Manhattan Plot")
    plt.scatter(np.arange(100), -np.log10(np.random.rand(100)))
    plt.xlabel("Chromosome Position")
    plt.ylabel("-log10(p-value)")
    plt.savefig("manhattan_plot.png")
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.title("QQ Plot")
    plt.scatter(np.sort(np.random.rand(100)), -np.log10(np.random.rand(100)))
    plt.xlabel("Expected -log10(p-value)")
    plt.ylabel("Observed -log10(p-value)")
    plt.savefig("qq_plot.png")
    plt.close()


# 生成PDF报告
def generate_report():
    pdf = PDF()
    pdf.add_font("SimHei", "", FONT_PATH, uni=True)
    pdf.set_font("SimHei", "", 14)
    pdf.add_page()

    # 系统发育树分析
    pdf.chapter_title("1. 系统发育树分析")
    pdf.chapter_body(
        "IQ-TREE2 是一款高效的系统发育树构建软件，支持多种替换模型和快速计算。它通过最大似然法（Maximum Likelihood）来推断系统发育树，并使用快速自举分析来评估树节点的可靠性。系统发育树显示了各个样本之间的进化关系，自举支持值表明了树节点的可靠性。下图展示了构建的系统发育树："
    )
    # generate_phylogenetic_tree_plot()
    pdf.add_image("test.png")
    pdf.check_page_break(100)

    # PCA分析
    pdf.chapter_title("2. 主成分分析")
    pdf.chapter_body(
        "PLINK 是一款常用的基因组数据分析工具，支持多种遗传统计分析。主成分分析（PCA）通过降维技术，将高维基因型数据投影到低维空间，以揭示样本间的遗传结构和群体分层。主成分分析结果展示了样本在前两个主成分空间中的分布。下图为PCA分析结果的散点图："
    )
    # generate_pca_plot()
    pdf.add_image("test.png")
    pdf.check_page_break(100)

    # ADMIXTURE分析
    pdf.chapter_title("3. 群体结构分析")
    pdf.chapter_body(
        "ADMIXTURE 是一款快速估计种群结构的软件，通过最大似然法估计每个个体的混合比例，揭示个体在不同种群中的遗传成分。群体结构分析结果展示了每个个体的混合比例。下图分别展示了使用Delta K方法选择的最佳K值和最佳K值的群体结构分析结果："
    )
    # generate_admixture_plot()
    pdf.add_image("test.png")
    pdf.check_page_break(100)

    # LD分析
    pdf.chapter_title("4. 连锁不平衡分析")
    pdf.chapter_body(
        "PopLDdecay 是一款专门用于大规模SNP数据的连锁不平衡（LD）分析工具。它通过计算SNP对之间的LD值（如r²），评估基因组范围内的LD衰减模式。连锁不平衡分析结果展示了不同距离范围内的LD值。下图为LD衰减曲线："
    )
    # generate_ld_decay_plot()
    pdf.add_image("test.png")
    pdf.check_page_break(100)

    # GWAS分析
    pdf.chapter_title("5. 性状关联分析")
    pdf.chapter_body(
        "GEMMA 是一款用于基因组宽关联分析（GWAS）和混合线性模型分析的软件。它通过混合线性模型（MLM）控制群体结构和亲缘关系，识别与性状相关的遗传变异。性状关联分析结果展示了与性状显著关联的SNP位点。下图为GWAS结果的Manhattan图和QQ图："
    )
    # generate_gwas_plot()
    pdf.add_image("test.png")
    pdf.add_image("test.png")

    pdf.output("analysis_report.pdf")


if __name__ == "__main__":
    typer.run(generate_report)
