import re
import sqlite3
from pathlib import Path

import requests
import typer
from fsspec import config
from loguru import logger

app = typer.Typer()


PROJECT_API = "http://183.221.124.252:8095/api/project/sequencing"

TIMEOUT = 60


@app.command()
def upload_one(
    report_path: Path,
    project_name: str,
    project_code: str,
    ref_genome: str,
    biochip_code: str,
):
    url = "http://183.221.124.252:8190/api/v1/capture"
    headers = {
        "Token": "14CF2E831D6F4B853F37B386D74B1EDA",
    }
    files = {
        "upload[]": open(
            f"{report_path}",
            "rb",
        )
    }
    data = {
        "project_name": f"{project_name}",
        "project_code": f"{project_code}",
        "ref_genome": f"{ref_genome}",
        "biochip_code": f"{biochip_code}",
    }
    response = requests.post(
        url, headers=headers, files=files, data=data, verify=False, timeout=TIMEOUT
    )
    # print(response.status_code, response.text)
    if response.status_code != 200:
        return False
    if (
        response.json()["errMessage"]
        and response.json()["errMessage"] != "业务错误: project already exists"
    ):
        return False
    return True


def chip_code_from_file(report_path: Path) -> str:
    report_name = report_path.name
    chip_code = re.sub(r"\.capture\.report\.tsv$", "", report_name)
    if re.search("[sS][0-9]+$", chip_code):
        chip_code = re.sub(r"[-_][sS][0-9]+$", "", chip_code)
    return chip_code


def fetch_project_name(project_code: str) -> str:
    """获取项目名"""
    params = {"projectId": project_code, "sample": "false"}

    response = requests.get(PROJECT_API, params=params, timeout=TIMEOUT)
    if response.status_code == 200:
        return response.json()["projectName"]
    logger.error(f"{project_code}: 请求失败，状态码：{response.status_code}")
    logger.error(response.text)
    return ""


def fetch_project_genome(project_dir: Path) -> str:
    """获取项目基因组"""
    project_code = project_dir.name
    for sh_path in project_dir.glob("*.sh"):
        pattern = rf"pipe_gts_{re.escape(project_code)}_(.+)\.sh"
        match = re.match(pattern, sh_path.name)
        if match:
            genome_version = match.group(1)
            if re.search("[-_]split$", genome_version):
                return re.sub("[-_]split$", "", genome_version)
            return genome_version
    return ""


def init_db(db_file: Path) -> None:
    db_file.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(db_file) as conn:
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS project_info (
                project_code TEXT PRIMARY KEY,
                chip_code TEXT,
                upload_status TEXT,
                last_modified TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """
        )


def add_or_update(
    db_file: Path, chip_code: str, project_code: str, upload_status: str
) -> None:
    with sqlite3.connect(db_file) as conn:
        cur = conn.cursor()
        cur.execute(
            """
        REPLACE INTO project_info (chip_code, project_code, upload_status, last_modified)
        VALUES (?, ?, ?, CURRENT_TIMESTAMP)
        """,
            (chip_code, project_code, upload_status),
        )
        conn.commit()


def query_upload_status(db_file: Path, project_code: str) -> str:
    with sqlite3.connect(db_file) as conn:
        cur = conn.cursor()
        cur.execute(
            """
        SELECT upload_status FROM project_info
        WHERE project_code = ?
        ORDER BY last_modified DESC
        LIMIT 1
        """,
            (project_code,),
        )
        res = cur.fetchone()
        if res is None:
            return "empty"
        return res[0]


@app.command()
def upload_batch(base_dir: Path, config_file: Path) -> None:
    log_dir = config_file.parent
    log_dir.mkdir(parents=True, exist_ok=True)
    logger.add(
        "log.txt", rotation="10 MB", retention="7 days", encoding="utf-8", level="INFO"
    )
    if config_file.exists():
        typer.echo(f"使用已有配置文件: {config_file}")
    else:
        init_db(config_file)
    for report_path in base_dir.glob("*/results/summary/bamdst/*.capture.report.tsv"):
        chip_code = chip_code_from_file(report_path)
        project_base_dir = report_path.parent.parent.parent.parent
        genome_version = fetch_project_genome(project_base_dir)
        project_code = report_path.parent.parent.parent.parent.name
        logger.info(f"开始处理: {project_code}")
        if not genome_version:
            logger.warning(f"未找到项目基因组版本: {project_base_dir}")
            continue
        project_name = fetch_project_name(project_code)
        if not project_name:
            logger.warning(f"未找到项目名称: {project_code}")
            continue
        upload_status = query_upload_status(
            db_file=config_file, project_code=project_code
        )
        if upload_status == "success":
            logger.info(f"项目 {project_name} 已在数据库中，跳过上传")
            continue
        success = upload_one(
            report_path=report_path,
            project_name=project_name,
            project_code=project_code,
            ref_genome=genome_version,
            biochip_code=chip_code,
        )
        if success:
            logger.success(f"上传成功: {project_code}")
            add_or_update(
                db_file=config_file,
                chip_code=chip_code,
                project_code=project_code,
                upload_status="success",
            )
        else:
            logger.error(f"上传失败: {project_code} | {chip_code}")
            add_or_update(
                db_file=config_file,
                chip_code=chip_code,
                project_code=project_code,
                upload_status="fail",
            )


if __name__ == "__main__":
    app()
