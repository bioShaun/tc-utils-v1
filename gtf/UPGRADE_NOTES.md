# 升级到 pytest 测试框架 - v2.0

## 概述

本次升级将测试框架从原生Python测试升级到pytest，为测试提供更强大的功能和更好的可维护性。

## 新增文件

### 测试文件

1. **`test_extract_gene_id_map_pytest.py`**
   - pytest版本的基因映射提取测试
   - 使用fixture和assertion语法
   - 3个测试用例全部通过

2. **`test_filter_gene_map_by_chromosome_pytest.py`**
   - pytest版本的染色体过滤测试
   - 4个测试用例全部通过

### 配置文件

3. **`pytest.ini`**
   - pytest主配置文件
   - 配置测试路径、输出格式、标记等

4. **`conftest.py`**
   - pytest共享配置文件
   - 可用于定义全局fixture和钩子

5. **`run_pytest.py`**
   - pytest测试运行器
   - 一键运行所有pytest测试

6. **`requirements-dev.txt`**
   - 开发依赖列表
   - 包含pytest及相关插件

7. **`pytest_guide.md`**
   - pytest使用指南
   - 包含快速入门和最佳实践

## 主要改进

### 1. 测试代码更简洁

**Before（原生Python）：**
```python
def test_something():
    result = do_something()
    if result != expected:
        print(f"Expected {expected}, got {result}")
        return False
    return True
```

**After（pytest）：**
```python
def test_something():
    assert do_something() == expected
```

### 2. 自动Fixture管理

```python
@pytest.fixture
def temp_dir():
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir)

def test_with_temp_dir(temp_dir):
    # 使用fixture，无需手动创建和清理
    assert True
```

### 3. 更详细的测试报告

pytest自动生成彩色输出和详细的测试信息：
- 测试通过/失败状态
- 测试运行时间
- 失败时的详细错误信息
- 测试覆盖率报告（可选）

### 4. 更灵活的测试运行

```bash
# 运行所有测试
pytest

# 运行特定测试
pytest test_extract_gene_id_map_pytest.py::test_with_group_map

# 详细输出
pytest -v

# 显示最慢的测试
pytest --durations=10
```

## 兼容性

- 保留所有原有的原生Python测试文件
- 两个测试框架可并行运行
- 原有测试脚本 (`run_all_tests.py`) 仍然有效

## 推荐使用方式

**新项目或新测试：** 使用pytest框架
**现有测试：** 可以逐步迁移或保持现有格式

## 安装和使用

```bash
# 安装依赖
pip install -r tests/requirements-dev.txt

# 运行pytest测试
cd tests
python run_pytest.py

# 或使用pytest命令
pytest -v
```

## 测试结果

✅ 所有7个测试用例通过：
- `test_extract_gene_id_map_pytest.py`: 3个测试
- `test_filter_gene_map_by_chromosome_pytest.py`: 4个测试

## 迁移计划

建议逐步将现有测试迁移到pytest框架：

1. ✅ 新增pytest测试文件（已完成）
2. ✅ 保留原生Python测试（已完成）
3. ✅ 验证pytest测试功能（已完成）
4. 🔄 可选：逐步迁移现有测试到pytest格式

## 优势总结

1. **更简洁**：无需手动返回True/False
2. **更清晰**：失败时自动显示错误信息
3. **更灵活**：支持fixture、参数化、标记等
4. **更专业**：工业级测试框架
5. **更高效**：更好的测试组织和发现
6. **更可扩展**：支持插件和自定义扩展

## 升级完成时间

2025年12月3日

## 下一步计划

1. 添加更多测试用例
2. 集成代码覆盖率测试
3. 添加CI/CD配置示例
4. 编写更详细的测试最佳实践指南
