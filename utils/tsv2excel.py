from pathlib import Path

import pandas as pd
import typer

def file2excel(input_path: Path, excel_path: Path):
    """
    将CSV或TSV文件转换为Excel文件
    
    参数:
        input_path: CSV或TSV文件路径
        excel_path: 输出的Excel文件路径
    """
    # 检测文件类型并选择适当的分隔符
    if str(input_path).lower().endswith('.tsv'):
        # TSV文件使用制表符分隔
        df = pd.read_table(input_path)
    elif str(input_path).lower().endswith('.csv'):
        # CSV文件使用逗号分隔
        df = pd.read_csv(input_path)
    else:
        # 尝试自动检测分隔符
        with open(input_path, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()
            
        if '\t' in first_line:
            df = pd.read_table(input_path)
        elif ',' in first_line:
            df = pd.read_csv(input_path)
        else:
            # 如果无法确定，默认使用pandas的自动检测
            df = pd.read_csv(input_path, sep=None, engine='python')
    
    # 转换为Excel格式
    df.to_excel(excel_path, index=False)
    
    return df.shape  # 返回数据框的行数和列数，便于确认转换结果

if __name__ == "__main__":
    typer.run(file2excel)
