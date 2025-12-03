# pytest 快速指南

## 安装

```bash
pip install pytest
```

## 运行测试

### 基本命令

```bash
# 运行当前目录及子目录下的所有测试
pytest

# 运行特定文件
pytest test_extract_gene_id_map_pytest.py

# 运行特定测试函数
pytest test_extract_gene_id_map_pytest.py::test_with_group_map

# 详细输出
pytest -v

# 显示最慢的10个测试
pytest --durations=10

# 显示详细失败信息
pytest -vv

# 在第一个失败时停止
pytest -x

# 运行所有测试但显示所有失败
pytest --maxfail=5
```

### 标记测试

pytest允许你标记测试，方便选择性地运行：

```bash
# 运行所有标记为"unit"的测试
pytest -m unit

# 运行所有未标记为"slow"的测试
pytest -m "not slow"

# 查看所有标记
pytest --markers
```

### 代码覆盖率

```bash
# 安装pytest-cov
pip install pytest-cov

# 运行测试并生成覆盖率报告
pytest --cov=your_module tests/

# 生成HTML覆盖率报告
pytest --cov=your_module --cov-report=html tests/
```

## 测试结构

### Fixture（夹具）

Fixture用于在测试前后设置和清理资源：

```python
import pytest

@pytest.fixture
def my_fixture():
    # 设置
    resource = setup_resource()
    yield resource
    # 清理
    cleanup_resource(resource)

def test_using_fixture(my_fixture):
    assert my_fixture.something
```

### 参数化测试

```python
@pytest.mark.parametrize("input,expected", [
    ("a", "A"),
    ("b", "B"),
    ("c", "C"),
])
def test_upper(input, expected):
    assert input.upper() == expected
```

## 断言

pytest使用原生Python断言，不需要特殊函数：

```python
def test_assertions():
    assert 1 + 1 == 2
    assert "hello" in "hello world"
    assert [1, 2, 3] == [1, 2, 3]
    assert {"a": 1} == {"a": 1}
```

## 最佳实践

1. **测试文件名**：`test_*.py`
2. **测试函数名**：`test_*()`
3. **使用fixture**管理资源（临时文件、数据库连接等）
4. **使用assert**而不是print语句
5. **每个测试独立**（不依赖其他测试）
6. **测试名称描述性强**（清晰表达测试意图）
7. **使用标记**区分不同类型的测试（unit, integration, slow等）
